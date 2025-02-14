import json
import os
import imaplib
import email
import re
import smtplib
import time
import random
import hashlib
import requests
from email.mime.text import MIMEText
from email.header import Header

import urllib3
from Bio import Entrez
from tenacity import retry, stop_after_attempt, wait_exponential, wait_fixed
import ssl
import certifi

# ==== 环境变量配置 ====
os.environ.update({
    'GMAIL_ADDRESS': 'jihaibiao012@gmail.com',
    'GMAIL_APP_PASSWORD': 'hbxaosexacavrars',
    'BAIDU_APP_ID': '20250214002273327',  # 替换为你的百度APP ID
    'BAIDU_SECRET_KEY': 'UszzWMkFZFzR8YmpRzPB'  # 替换为你的百度密钥
})

# 验证环境变量
EMAIL = os.getenv('GMAIL_ADDRESS')
PASSWORD = os.getenv('GMAIL_APP_PASSWORD')
BAIDU_APP_ID = os.getenv('BAIDU_APP_ID')
BAIDU_SECRET_KEY = os.getenv('BAIDU_SECRET_KEY')
assert all([EMAIL, PASSWORD, BAIDU_APP_ID, BAIDU_SECRET_KEY]), "环境变量未正确设置"

# ==== 其他配置保持不变 ====
IMAP_SERVER = 'imap.gmail.com'
SMTP_SERVER = 'smtp.gmail.com'
IMAP_TIMEOUT = 60
Entrez.email = EMAIL

ssl._create_default_https_context = ssl._create_unverified_context
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)


class BaiduTranslator:
    """百度翻译API封装类"""

    def __init__(self, app_id, secret_key):
        if not app_id or not secret_key:
            raise ValueError("需要有效的百度翻译API凭证")

        self.app_id = app_id
        self.secret_key = secret_key
        self.base_url = "https://fanyi-api.baidu.com/api/trans/vip/translate"

        # 配置带证书验证的会话
        self.session = requests.Session()
        self.session.verify = certifi.where()
        self.timeout = 15

    def _generate_sign(self, query, salt):
        """生成百度API签名"""
        sign_str = f"{self.app_id}{query}{salt}{self.secret_key}"
        return hashlib.md5(sign_str.encode('utf-8')).hexdigest()

    @retry(stop=stop_after_attempt(3),
           wait=wait_exponential(multiplier=1, min=2, max=10),
           reraise=True)
    def translate(self, text, target_lang='zh', max_length=4500):
        """执行翻译"""
        if not text:
            return ""

        # 百度API单次请求最大6000字节
        text = text[:max_length]

        try:
            salt = str(random.randint(32768, 65536))
            sign = self._generate_sign(text, salt)

            params = {
                'q': text,
                'from': 'en',
                'to': target_lang,
                'appid': self.app_id,
                'salt': salt,
                'sign': sign
            }

            response = self.session.get(
                self.base_url,
                params=params,
                timeout=self.timeout
            )

            response.raise_for_status()
            result = response.json()

            if 'error_code' in result:
                error_msg = result.get('error_msg', '未知错误')
                print(f"百度翻译错误 (代码{result['error_code']}): {error_msg}")
                return text[:200] + "[...]"

            return '\n'.join([item['dst'] for item in result['trans_result']])

        except requests.exceptions.SSLError as e:
            print(f"SSL验证失败: {str(e)}")
            raise
        except requests.exceptions.Timeout:
            print(f"请求超时 ({self.timeout}s)")
            raise
        except requests.exceptions.RequestException as e:
            print(f"请求失败: {str(e)}")
            raise
        except KeyError as e:
            print(f"响应格式错误: {str(e)}")
            return text

    def close(self):
        self.session.close()


@retry(stop=stop_after_attempt(3), wait=wait_fixed(10))
def fetch_stork_emails():
    """获取Stork推送的未读邮件（优化超时和重试）"""
    try:
        print("正在连接IMAP服务器...")
        mail = imaplib.IMAP4_SSL(IMAP_SERVER, timeout=IMAP_TIMEOUT)  # 显式设置超时
        mail.login(EMAIL, PASSWORD)
        mail.select('inbox')
        _, data = mail.search(None, 'UNSEEN', '(FROM "support@storkapp.me")')
        email_ids = data[0].split()
        print(f"IMAP连接成功，找到 {len(email_ids)} 封未读邮件")
        return mail, email_ids
    except Exception as e:
        print(f"IMAP连接失败: {str(e)}")
        raise


@retry(stop=stop_after_attempt(3), wait=wait_exponential(multiplier=1, min=4, max=10))
def send_summary_email(content):
    """使用Gmail发送汇总邮件（添加详细错误日志）"""
    try:
        msg = MIMEText(content, 'html', 'utf-8')
        msg['Subject'] = Header('每日论文摘要推送', 'utf-8')
        msg['From'] = Header(f"论文助手 <{EMAIL}>", 'utf-8')
        msg['To'] = 'jihaibiao012@163.com'

        with smtplib.SMTP(SMTP_SERVER, 587, timeout=30) as server:  # 设置SMTP超时
            server.ehlo()
            server.starttls()
            server.login(EMAIL, PASSWORD)
            server.sendmail(EMAIL, ['jihaibiao012@163.com'], msg.as_string())
        print("邮件发送成功！")
    except smtplib.SMTPException as e:
        print(f"SMTP协议错误: {str(e)}")
        raise
    except Exception as e:
        print(f"邮件发送未知错误: {str(e)}")
        raise


def extract_paper_info(msg):
    body = ""
    if msg.is_multipart():
        for part in msg.walk():
            if part.get_content_type() == "text/plain":
                body = part.get_payload(decode=True).decode('utf-8', errors='ignore')
                break
    else:
        body = msg.get_payload(decode=True).decode('utf-8', errors='ignore')

    print("邮件正文内容：", body)  # 调试输出

    # 改进后的正则表达式
    pattern = r"([A-Z][^\.]+?\.)\s+by\s+([^\(]+)\s+\((\d{4})\)\s+([^\(]+)\s+\(impact factor:\s+([0-9.]+)\)[^P]*PMID:\s+(\d+)\s+doi:\s+([^\s]+)"
    matches = re.findall(pattern, body, re.DOTALL)

    print("提取到的论文信息：", matches)  # 调试输出

    return [{'title': match[0].strip(), 'pmid': match[5], 'doi': match[6]} for match in matches]


@retry(stop=stop_after_attempt(3), wait=wait_exponential(multiplier=1, min=4, max=10))
def get_abstract_from_pubmed(pmid):
    """通过 PubMed API 获取论文摘要"""
    try:
        handle = Entrez.efetch(db="pubmed", id=pmid, rettype="abstract", retmode="text")
        abstract = handle.read()
        handle.close()
        return abstract.strip() if abstract else "未找到摘要"
    except Exception as e:
        print(f"获取摘要失败 (PMID: {pmid}): {str(e)}")
        return "摘要获取失败"


@retry(stop=stop_after_attempt(3), wait=wait_exponential(multiplier=1, min=4, max=10))
def main():
    try:
        print("===== 脚本启动 =====")

        # 1. 连接 IMAP 服务器并获取邮件
        print("尝试连接 IMAP 服务器...")
        mail, email_ids = fetch_stork_emails()
        print(f"连接成功，共找到 {len(email_ids)} 封未读邮件")

        translator = BaiduTranslator(BAIDU_APP_ID, BAIDU_SECRET_KEY)
        all_translations = []

        # 2. 处理每封邮件
        for e_id in email_ids:
            print(f"正在处理邮件 ID: {e_id}")
            _, data = mail.fetch(e_id, '(RFC822)')
            msg = email.message_from_bytes(data[0][1])

            # 3. 提取论文信息
            print("提取论文信息中...")
            papers = extract_paper_info(msg)
            print(f"提取到 {len(papers)} 篇论文")

            # 4. 处理每篇论文
            for paper in papers:
                try:
                    print(f"正在处理论文: {paper['title']} (PMID: {paper['pmid']})")

                    # 获取摘要
                    abstract = get_abstract_from_pubmed(paper['pmid'])
                    print(f"获取到摘要: {abstract[:50]}...")  # 只打印前 50 个字符

                    # 翻译
                    zh_title = translator.translate(paper['title'])
                    print(f"翻译后的标题: {zh_title}")
                    zh_abstract = translator.translate(abstract) if abstract else "无可用摘要"
                    print(f"翻译后的摘要: {zh_abstract[:50]}...")

                    # 添加到邮件内容
                    all_translations.append(f"""
                    <h3>{zh_title}</h3>
                    <p><b>原文标题:</b> {paper['title']}</p>
                    <p><b>摘要:</b> {zh_abstract}</p>
                    <p><b>PMID:</b> {paper['pmid']} | <b>DOI:</b> {paper['doi']}</p>
                    <hr>
                    """)
                except Exception as e:
                    print(f"处理论文失败: {str(e)}")

        # 5. 发送邮件
        if all_translations:
            print("生成摘要邮件内容...")
            html_content = f"""
            <html>
                <body>
                    <h2>今日 Stork 论文推送摘要 ({len(all_translations)} 篇)</h2>
                    {"".join(all_translations)}
                </body>
            </html>
            """
            print("邮件内容准备完成，准备发送邮件...")
            send_summary_email(html_content)
        else:
            print("今日无新论文推送")

    except Exception as e:
        print(f"!!! 主程序异常: {str(e)} !!!")
        raise
    finally:
        if 'mail' in locals():
            try:
                mail.close()
                mail.logout()
                print("IMAP 连接已关闭")
            except Exception as e:
                print(f"关闭 IMAP 连接时出错: {str(e)}")
        print("===== 脚本结束 =====")


if __name__ == "__main__":
    main()