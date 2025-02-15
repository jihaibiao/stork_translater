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
import xml.etree.ElementTree as ET

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
    """从邮件提取PMID和影响因子"""
    body = ""
    if msg.is_multipart():
        for part in msg.walk():
            content_type = part.get_content_type()
            if content_type == "text/html":
                body = part.get_payload(decode=True).decode('utf-8', errors='ignore')
                break
    else:
        body = msg.get_payload(decode=True).decode('utf-8', errors='ignore')

    soup = BeautifulSoup(body, 'html.parser')
    text = ' '.join(soup.stripped_strings)

    papers = []
    pattern = re.compile(
        r'PMID:\s+(?P<pmid>\d+).*?impact\s+factor:\s*(?P<impact_factor>\d+\.?\d*)',
        re.DOTALL
    )

    for match in pattern.finditer(text):
        papers.append({
            'pmid': match.group('pmid'),
            'impact_factor': match.group('impact_factor')
        })

    print(f"提取到 {len(papers)} 篇论文的PMID和影响因子")
    return papers


@retry(stop=stop_after_attempt(3), wait=wait_fixed(10))
def get_pubmed_details(pmid):
    """从PubMed获取完整文献信息"""
    try:
        handle = Entrez.efetch(db="pubmed", id=pmid, retmode="xml")
        data = handle.read()
        handle.close()

        root = ET.fromstring(data)
        article = root.find('.//PubmedArticle')

        # 提取标题
        title = article.find('.//ArticleTitle').text.strip()

        # 提取作者
        authors = []
        for author in article.findall('.//Author'):
            lastname = author.find('LastName').text if author.find('LastName') is not None else ''
            forename = author.find('ForeName').text if author.find('ForeName') is not None else ''
            if lastname or forename:
                authors.append(f"{forename} {lastname}".strip())
        authors_str = ", ".join(authors[:3]) + (" et al." if len(authors) > 3 else "")

        # 提取期刊信息
        journal = article.find('.//Journal/Title').text
        year = article.find('.//PubDate/Year').text if article.find('.//PubDate/Year') is not None else ''

        # 提取DOI
        doi = ""
        for id in article.findall('.//ArticleId'):
            if id.attrib.get('IdType') == 'doi':
                doi = id.text
                break

        return {
            'title': title,
            'authors': authors_str,
            'journal': journal,
            'year': year,
            'doi': doi
        }
    except Exception as e:
        print(f"❌ 获取PubMed数据失败 PMID {pmid}: {str(e)}")
        return None


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
                print(f"\n🔍 处理PMID: {paper['pmid']}")
                pubmed_data = get_pubmed_details(paper['pmid'])
                if not pubmed_data:
                    continue

                # 获取摘要
                abstract = get_abstract_from_pubmed(paper['pmid'])

                # 翻译处理
                zh_title = translator.translate(pubmed_data['title'])
                zh_abstract = translator.translate(abstract) if abstract else "无可用摘要"

                # 合并数据
                full_data = {
                    **pubmed_data,
                    **paper,
                    'abstract': abstract,
                    'zh_title': zh_title,
                    'zh_abstract': zh_abstract
                }

                # 构建内容
                all_translations.append(f"""
                <div style="margin-bottom: 2rem; padding: 1.5rem; background: #f8faff; border-radius: 8px; box-shadow: 0 2px 12px rgba(28,87,223,0.1); border-left: 4px solid #1a73e8;">
                    <h3 style="color: #1a3d6d; margin: 0 0 0.8rem 0; font-size: 1.1rem; line-height: 1.4;">
                        {full_data['zh_title']}
                    </h3>
                    <div style="color: #4a5568; line-height: 1.6;">
                        <p style="margin: 0.4rem 0;">
                            <span style="font-weight: 600;">📖 原文标题:</span> 
                            <span style="color: #2d3748;">{full_data['title']}</span>
                        </p>
                        <p style="margin: 0.4rem 0;">
                            <span style="font-weight: 600;">👥 作者:</span> 
                            {full_data['authors']}
                        </p>
                        <p style="margin: 0.4rem 0;">
                            <span style="font-weight: 600;">🏛️ 期刊:</span> 
                            {full_data['journal']} ({full_data['year']}, IF: {full_data['impact_factor']})
                        </p>
                        <div style="margin: 1rem 0; padding: 0.8rem; background: white; border-radius: 6px; border: 1px solid #e2e8f0;">
                            <span style="font-weight: 600;">📄 摘要:</span> 
                            <div style="color: #4a5568; margin-top: 0.4rem;">
                                {full_data['zh_abstract']}
                            </div>
                        </div>
                        <div style="margin-top: 1rem;">
                            <a href="https://pubmed.ncbi.nlm.nih.gov/{full_data['pmid']}" 
                               target="_blank"
                               style="display: inline-block; padding: 6px 12px; background: #1a73e8; color: white; border-radius: 4px; text-decoration: none; margin-right: 8px;">
                               PubMed
                            </a>
                            {f'<a href="https://doi.org/{full_data["doi"]}" target="_blank" style="display: inline-block; padding: 6px 12px; background: #38a169; color: white; border-radius: 4px; text-decoration: none;">Full Text</a>' if full_data['doi'] else ''}
                        </div>
                    </div>
                </div>
                """)

        # 剩余部分保持不变...
        if all_translations:
            html_content = f"""
            <html>
                <body style="font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, sans-serif; max-width: 700px; margin: 0 auto; padding: 2rem 1rem; background-color: #f7fafc;">
                    <header style="text-align: center; margin-bottom: 2.5rem;">
                        <h1 style="color: #1a365d; margin: 0 0 0.5rem 0; font-size: 1.8rem; display: flex; align-items: center; gap: 0.8rem; justify-content: center;">
                            <span style="background: #1a73e8; color: white; padding: 6px 12px; border-radius: 6px;">📰 今日文献</span>
                            <span>推送 ({len(all_translations)}篇)</span>
                        </h1>
                    </header>
                    {"".join(all_translations)}
                    <footer style="margin-top: 3rem; text-align: center; color: #718096; font-size: 0.85rem; padding-top: 1.5rem; border-top: 1px solid #e2e8f0;">
                        🚀 由论文助手自动生成 | 📧 反馈请联系 {EMAIL}
                    </footer>
                </body>
            </html>
            """
            send_summary_email(html_content)
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


# 其他辅助函数保持不变...
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

        # 使用更稳定的SMTP_SSL连接
        context = ssl.create_default_context()
        with smtplib.SMTP_SSL(SMTP_SERVER, 465, timeout=60, context=context) as server:
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

if __name__ == "__main__":
    main()