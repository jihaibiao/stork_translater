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

# ==== 核心修复部分 ====

os.environ.update({
    'GMAIL_ADDRESS': 'jihaibiao012@gmail.com',
    'GMAIL_APP_PASSWORD': 'hbxaosexacavrars',
    'TENCENT_APP_ID': 'AKIDGOJ7nJl00oa9VeaMtCF8sPksseBTWLwI',
    'TENCENT_APP_KEY': '93jAMK7lDSSEl25tSntkNaidTHvMIF5j'
})

# 验证变量（确保后续代码正确读取）
EMAIL = os.getenv('GMAIL_ADDRESS')
PASSWORD = os.getenv('GMAIL_APP_PASSWORD')
assert EMAIL and PASSWORD, "环境变量未正确设置"

# ==== 配置信息 ====
IMAP_SERVER = 'imap.gmail.com'
SMTP_SERVER = 'smtp.gmail.com'
IMAP_TIMEOUT = 60  # 增加IMAP超时时间
Entrez.email = EMAIL

# 正确设置 SSL 上下文
ssl._create_default_https_context = ssl._create_unverified_context
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)


class TencentTranslator:
    """腾讯翻译API封装类 (包含SSL证书验证和详细错误处理)"""

    def __init__(self, app_id, app_key):
        if not app_id or not app_key:
            raise ValueError("腾讯翻译API需要有效的app_id和app_key")

        self.app_id = app_id
        self.app_key = app_key
        self.base_url = "https://api.ai.qq.com/fcgi-bin/nlp/nlp_texttranslate"

        # 创建带证书验证的会话
        self.session = requests.Session()
        self.session.verify = certifi.where()  # 使用certifi的CA证书包
        self.timeout = 15  # 请求超时时间

    def _get_sign(self, params):
        """生成请求签名 (修正排序和格式问题)"""
        # 确保参数排序符合腾讯API要求
        sorted_params = sorted(params.items(), key=lambda x: x[0])
        param_str = '&'.join([f"{k}={v}" for k, v in sorted_params])
        sign_str = f"{param_str}&app_key={self.app_key}"

        # 计算MD5时使用UTF-8编码
        return hashlib.md5(sign_str.encode('utf-8')).hexdigest().upper()

    @retry(stop=stop_after_attempt(3),
           wait=wait_exponential(multiplier=1, min=2, max=10),
           reraise=True)
    def translate(self, text, target_lang='zh', max_length=4500):
        """执行翻译 (包含详细错误处理)"""
        if not text:
            return ""

        # 分割过长的文本 (腾讯API限制5000字节)
        text = text[:max_length]

        try:
            # 1. 准备签名参数
            params = {
                'app_id': self.app_id,
                'time_stamp': str(int(time.time())),
                'nonce_str': ''.join(random.choices('abcdefghijklmnopqrstuvwxyz0123456789', k=16)),
                'text': text,
                'source': 'en',
                'target': target_lang
            }
            params['sign'] = self._get_sign(params)

            # 2. 发送带证书验证的请求
            response = self.session.post(
                self.base_url,
                data=params,
                timeout=self.timeout,
                headers={'Content-Type': 'application/x-www-form-urlencoded'}
            )

            # 3. 处理响应
            response.raise_for_status()
            result = response.json()

            if result.get('ret') == 0:
                return result['data']['target_text']
            else:
                error_msg = result.get('msg', '未知错误')
                print(f"腾讯翻译API错误 (代码{result['ret']}): {error_msg}")
                return text[:200] + "[...]"  # 返回部分原文防止信息丢失

        except requests.exceptions.SSLError as e:
            print(f"SSL验证失败: {str(e)}")
            print("尝试更新certifi证书包:pip install --upgrade certifi")
            raise
        except requests.exceptions.Timeout:
            print(f"请求超时 ({self.timeout}s),请检查网络连接")
            raise
        except requests.exceptions.RequestException as e:
            print(f"请求失败: {str(e)}")
            raise
        except KeyError as e:
            print(f"响应数据格式错误,缺少字段: {str(e)}")
            return text
        except json.JSONDecodeError:
            print("响应不是有效JSON格式")
            return text

    def close(self):
        """关闭会话连接"""
        self.session.close()

        # ...处理响应...
# ==== 关键修改1：为IMAP连接添加超时和重试 ====
@retry(stop=stop_after_attempt(3), wait=wait_fixed(10))  # 每10秒重试一次，最多3次
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

# ==== 关键修改2：增强SMTP错误处理 ====
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

        translator = TencentTranslator(os.getenv('TENCENT_APP_ID'), os.getenv('TENCENT_APP_KEY'))
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
