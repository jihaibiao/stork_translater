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
from bs4 import BeautifulSoup

# ==== 环境变量配置 ====
os.environ.update({
    'GMAIL_ADDRESS': 'jihaibiao012@gmail.com',
    'GMAIL_APP_PASSWORD': 'hbxaosexacavrars',
    'BAIDU_APP_ID': '20250214002273327',
    'BAIDU_SECRET_KEY': 'UszzWMkFZFzR8YmpRzPB'
})

# 验证环境变量
EMAIL = os.getenv('GMAIL_ADDRESS')
PASSWORD = os.getenv('GMAIL_APP_PASSWORD')
BAIDU_APP_ID = os.getenv('BAIDU_APP_ID')
BAIDU_SECRET_KEY = os.getenv('BAIDU_SECRET_KEY')
assert all([EMAIL, PASSWORD, BAIDU_APP_ID, BAIDU_SECRET_KEY]), "环境变量未正确设置"

# ==== 配置信息 ====
IMAP_SERVER = 'imap.gmail.com'
SMTP_SERVER = 'smtp.gmail.com'
IMAP_TIMEOUT = 60
Entrez.email = EMAIL

ssl._create_default_https_context = ssl._create_unverified_context
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)


class BaiduTranslator:
    """百度翻译API封装（修复版）"""

    def __init__(self, app_id, secret_key):
        self.app_id = app_id
        self.secret_key = secret_key
        self.base_url = "https://fanyi-api.baidu.com/api/trans/vip/translate"
        self.session = requests.Session()
        self.session.verify = certifi.where()
        self.timeout = 15

    def _generate_sign(self, query, salt):
        sign_str = f"{self.app_id}{query}{salt}{self.secret_key}"
        return hashlib.md5(sign_str.encode('utf-8')).hexdigest()

    @retry(stop=stop_after_attempt(3), wait=wait_exponential(multiplier=1, min=2, max=10))
    def translate(self, text, target_lang='zh', max_length=4500):
        if not text:
            return ""

        text = text[:max_length]
        salt = str(random.randint(32768, 65536))

        try:
            params = {
                'q': text,
                'from': 'en',
                'to': target_lang,
                'appid': self.app_id,
                'salt': salt,
                'sign': self._generate_sign(text, salt)
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
                raise ValueError(f"百度翻译错误 {result['error_code']}: {error_msg}")

            return ' '.join([item['dst'] for item in result['trans_result']])

        except Exception as e:
            print(f"翻译失败: {str(e)}")
            return f"[翻译失败] {text[:200]}"


def extract_paper_info(msg):
    """论文信息提取（增强正则表达式版）"""
    # 解析邮件正文
    body = ""
    if msg.is_multipart():
        for part in msg.walk():
            content_type = part.get_content_type()
            if content_type == "text/html":
                body = part.get_payload(decode=True).decode('utf-8', errors='ignore')
                break
    else:
        body = msg.get_payload(decode=True).decode('utf-8', errors='ignore')

    # 使用BeautifulSoup清理内容
    soup = BeautifulSoup(body, 'html.parser')
    text = ' '.join(soup.stripped_strings)

    # 调试输出
    with open("debug_email.txt", "w", encoding="utf-8") as f:
        f.write(text)

    # 优化后的正则表达式（精确匹配论文条目）
    pattern = re.compile(
        r'(?P<title>[A-Z][^\.]+?\.)\s+by\s+'  # 标题（以大写字母开头，句号结尾）
        r'(?P<authors>.+?)\s+\('  # 作者
        r'(?P<year>\d{4})\)\s+'  # 年份
        r'(?P<journal>.+?)\s+\(impact\s+factor:\s*'
        r'(?P<impact_factor>\d+\.?\d*)\)'  # 影响因子
        r'.*?PMID:\s+(?P<pmid>\d+)\s+'  # PMID
        r'doi:\s+(?P<doi>10\.\S+)',  # DOI（以10.开头）
        re.DOTALL | re.IGNORECASE
    )

    papers = []
    for match in pattern.finditer(text):
        papers.append({
            'title': match.group('title').strip(),
            'pmid': match.group('pmid'),
            'doi': match.group('doi'),
            'authors': match.group('authors'),
            'year': match.group('year'),
            'journal': match.group('journal'),
            'impact_factor': match.group('impact_factor')
        })

    print(f"提取到 {len(papers)} 篇论文")
    return papers


@retry(stop=stop_after_attempt(3), wait=wait_fixed(10))
def fetch_stork_emails():
    try:
        print("🔄 连接IMAP服务器...")
        mail = imaplib.IMAP4_SSL(IMAP_SERVER, timeout=IMAP_TIMEOUT)
        mail.login(EMAIL, PASSWORD)
        mail.select('inbox')
        _, data = mail.search(None, 'UNSEEN', '(FROM "support@storkapp.me")')
        email_ids = data[0].split()
        print(f"✅ 找到 {len(email_ids)} 封未读邮件")
        return mail, email_ids
    except Exception as e:
        print(f"❌ IMAP连接失败: {str(e)}")
        raise


@retry(stop=stop_after_attempt(3), wait=wait_exponential(multiplier=1, min=4, max=10))
def send_summary_email(content):
    try:
        msg = MIMEText(content, 'html', 'utf-8')
        msg['Subject'] = Header('📚 每日论文摘要推送', 'utf-8')
        msg['From'] = Header(f"论文助手 <{EMAIL}>", 'utf-8')
        msg['To'] = 'jihaibiao012@163.com'

        with smtplib.SMTP(SMTP_SERVER, 587, timeout=30) as server:
            server.ehlo()
            server.starttls()
            server.login(EMAIL, PASSWORD)
            server.sendmail(EMAIL, ['jihaibiao012@163.com'], msg.as_string())
        print("📧 邮件发送成功！")
    except Exception as e:
        print(f"❌ 邮件发送失败: {str(e)}")
        raise


@retry(stop=stop_after_attempt(3), wait=wait_exponential(multiplier=1, min=4, max=10))
def get_abstract_from_pubmed(pmid):
    """修复摘要获取的decode错误"""
    try:
        handle = Entrez.efetch(db="pubmed", id=pmid, rettype="abstract", retmode="text")
        abstract = handle.read()  # 直接获取字符串内容
        handle.close()
        return abstract.strip() or "未找到摘要"
    except Exception as e:
        print(f"❌ 获取摘要失败 PMID {pmid}: {str(e)}")
        return "摘要获取失败"


def main():
    try:
        print("\n=== 🚀 论文助手开始运行 ===")
        mail, email_ids = fetch_stork_emails()
        translator = BaiduTranslator(BAIDU_APP_ID, BAIDU_SECRET_KEY)
        all_translations = []

        for e_id in email_ids:
            print(f"\n📨 处理邮件 {e_id.decode()}...")
            _, data = mail.fetch(e_id, '(RFC822)')
            msg = email.message_from_bytes(data[0][1])

            papers = extract_paper_info(msg)
            if not papers:
                print("⚠️ 未发现有效论文信息")
                continue

            for paper in papers:
                print(f"\n🔍 处理论文: {paper['title']}")
                abstract = get_abstract_from_pubmed(paper['pmid'])

                # 翻译处理
                zh_title = translator.translate(paper['title'])
                zh_abstract = translator.translate(abstract) if abstract else "无可用摘要"

                # 构建内容
                all_translations.append(f"""
                <div style="margin-bottom: 2rem; padding: 1rem; border-left: 4px solid #2196F3;">
                    <h3 style="color: #2c3e50; margin-top: 0;">{zh_title}</h3>
                    <p><strong>📖 原文标题:</strong> {paper['title']}</p>
                    <p><strong>👥 作者:</strong> {paper['authors']} ({paper['year']})</p>
                    <p><strong>🏛️ 期刊:</strong> {paper['journal']} (IF: {paper['impact_factor']})</p>
                    <p><strong>📄 摘要:</strong> {zh_abstract}</p>
                    <p style="font-size: 0.9em; color: #666;">
                        <a href="https://pubmed.ncbi.nlm.nih.gov/{paper['pmid']}" target="_blank">PubMed</a> | 
                        <a href="https://doi.org/{paper['doi']}" target="_blank">全文链接</a>
                    </p>
                </div>
                """)

        if all_translations:
            html_content = f"""
            <html>
                <body style="font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; max-width: 800px; margin: auto;">
                    <h1 style="color: #2c3e50; border-bottom: 2px solid #2196F3; padding-bottom: 0.5rem;">
                        📰 今日文献推送 ({len(all_translations)}篇)
                    </h1>
                    {"".join(all_translations)}
                    <footer style="margin-top: 2rem; text-align: center; color: #666; font-size: 0.9em;">
                        🚀 由文献鸟助手自动生成 | 📧 有问题请联系 {EMAIL}
                    </footer>
                </body>
            </html>
            """
            send_summary_email(html_content)
        else:
            print("ℹ️ 今日无新论文需要处理")

    except Exception as e:
        print(f"\n❌ 发生严重错误: {str(e)}")
    finally:
        if 'mail' in locals():
            try:
                mail.close()
                mail.logout()
            except:
                pass
        print("\n=== 🏁 运行结束 ===")


if __name__ == "__main__":
    main()