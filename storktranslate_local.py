
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
# 新增在文件顶部（所有import之后）
from dotenv import load_dotenv
load_dotenv()  # 加载本地.env文件
# ==== 环境变量配置 ====


# 验证环境变量
EMAIL = os.getenv('GMAIL_ADDRESS')
PASSWORD = os.getenv('GMAIL_APP_PASSWORD')
BAIDU_APP_ID = os.getenv('BAIDU_APP_ID')
BAIDU_SECRET_KEY = os.getenv('BAIDU_SECRET_KEY')
RECIPIENT_EMAIL = os.getenv('RECIPIENT_EMAIL')
# 验证环境变量（添加RECIPIENT_EMAIL检查）
assert all([EMAIL, PASSWORD, BAIDU_APP_ID, BAIDU_SECRET_KEY, RECIPIENT_EMAIL]), "环境变量未正确设置"

# ==== 配置信息 ====
IMAP_SERVER = 'imap.gmail.com'
SMTP_SERVER = 'smtp.gmail.com'
IMAP_TIMEOUT = 60
Entrez.email = EMAIL

ssl._create_default_https_context = ssl._create_unverified_context
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)


class BaiduTranslator:
    """百度翻译API封装"""

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
    """增强版文献信息提取"""
    body = ""
    if msg.is_multipart():
        for part in msg.walk():
            content_type = part.get_content_type()
            if content_type == "text/html":
                body = part.get_payload(decode=True).decode('utf-8', errors='ignore')
                break
    else:
        body = msg.get_payload(decode=True).decode('utf-8', errors='ignore')

    # 优化匹配模式，支持多文献
    pattern = re.compile(
        r'PMID:\s*(?P<pmid>\d+).*?'
        r'impact\s*factor:\s*(?P<impact_factor>\d+\.?\d*)',
        re.DOTALL
    )

    papers = []
    seen_pmids = set()

    for match in re.finditer(pattern, body.replace('\n', ' ')):
        pmid = match.group('pmid')
        if pmid not in seen_pmids:
            papers.append({
                'pmid': pmid,
                'impact_factor': match.group('impact_factor')
            })
            seen_pmids.add(pmid)

    print(f"提取到 {len(papers)} 篇论文信息")
    return papers


def extract_paper_info(msg):
    """增强版多文献提取（支持多关键词分组）"""
    body = ""
    if msg.is_multipart():
        for part in msg.walk():
            content_type = part.get_content_type()
            if content_type == "text/html":
                body = part.get_payload(decode=True).decode('utf-8', errors='ignore')
                break
    else:
        body = msg.get_payload(decode=True).decode('utf-8', errors='ignore')

    # 清理HTML标签并保留结构
    soup = BeautifulSoup(body, 'html.parser')
    text = soup.get_text(separator=' ', strip=True)

    # 优化匹配模式（支持多关键词分组）
    pattern = re.compile(
        r'PMID:\s*(?P<pmid>\d+).*?'  # 匹配PMID
        r'(?:impact\s*factor:\s*(?P<impact_factor>\d+\.?\d*))?.*?'  # 匹配影响因子（可选）
        r'(?:doi:\s*(?P<doi>10\.\S+))?',  # 匹配DOI（可选）
        re.DOTALL | re.IGNORECASE
    )

    papers = []
    current_pmid = None
    current_data = {}

    # 分段处理文献条目
    for section in re.split(r'(?=\bPMID:\s*\d+)', text):
        # 提取核心信息
        pmid_match = re.search(r'PMID:\s*(\d+)', section)
        if pmid_match:
            if current_pmid:  # 保存上一条记录
                papers.append(current_data)

            current_pmid = pmid_match.group(1)
            current_data = {
                'pmid': current_pmid,
                'impact_factor': 'N/A',
                'doi': 'N/A'
            }

            # 提取影响因子
            if_match = re.search(r'impact\s*factor:\s*(\d+\.?\d*)', section)
            if if_match:
                current_data['impact_factor'] = if_match.group(1)

            # 提取DOI
            doi_match = re.search(r'doi:\s*(10\.\S+)', section)
            if doi_match:
                current_data['doi'] = doi_match.group(1)

    # 添加最后一条记录
    if current_pmid:
        papers.append(current_data)

    # 去重逻辑优化
    seen = set()
    unique_papers = []
    for p in papers:
        if p['pmid'] not in seen:
            seen.add(p['pmid'])
            unique_papers.append(p)

    print(f"提取到 {len(unique_papers)} 篇论文信息")
    return unique_papers


@retry(stop=stop_after_attempt(3), wait=wait_fixed(10))
def get_pubmed_details(pmid):
    """从PubMed获取完整元数据"""
    try:
        handle = Entrez.efetch(db="pubmed", id=pmid, retmode="xml")
        data = handle.read()
        handle.close()

        root = ET.fromstring(data)
        article = root.find('.//PubmedArticle')

        # 提取元数据
        title = article.find('.//ArticleTitle').text.strip()
        journal = article.find('.//Journal/Title').text
        year = article.find('.//PubDate/Year').text if article.find('.//PubDate/Year') else ''
        doi = next((id.text for id in article.findall('.//ArticleId') if id.get('IdType') == 'doi'), '')

        # 处理作者信息
        authors = []
        for author in article.findall('.//Author'):
            last = author.findtext('LastName', '')
            fore = author.findtext('ForeName', '')
            if last or fore:
                authors.append(f"{fore} {last}".strip())
        author_str = ", ".join(authors[:3]) + (" et al." if len(authors) > 3 else "")

        return {
            'title': title,
            'authors': author_str,
            'journal': journal,
            'year': year,
            'doi': doi
        }
    except Exception as e:
        print(f"❌ PubMed数据获取失败 PMID {pmid}: {str(e)}")
        return None


@retry(stop=stop_after_attempt(3), wait=wait_exponential(multiplier=1, min=4, max=10))
def get_abstract_from_pubmed(pmid):
    """获取结构化摘要"""
    try:
        handle = Entrez.efetch(db="pubmed", id=pmid, retmode="xml")
        data = handle.read()
        handle.close()

        root = ET.fromstring(data)
        abstract = []

        # 保留结构化标签
        for elem in root.findall('.//AbstractText'):
            label = elem.get('Label', '').upper()
            text = elem.text.strip() if elem.text else ''
            if text:
                abstract.append(f"{label + ': ' if label else ''}{text}")

        return '\n\n'.join(abstract) or "未找到摘要"
    except Exception as e:
        print(f"❌ 摘要获取失败 PMID {pmid}: {str(e)}")
        return "摘要获取失败"

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


# ==== 修改邮件发送函数 ====
@retry(stop=stop_after_attempt(3), wait=wait_exponential(multiplier=1, min=4, max=10))
def send_summary_email(content):
    try:
        msg = MIMEText(content, 'html', 'utf-8')
        msg['Subject'] = Header('📚 每日论文摘要推送', 'utf-8')
        msg['From'] = Header(f"论文助手 <{EMAIL}>", 'utf-8')
        msg['To'] = RECIPIENT_EMAIL  # 改为使用环境变量

        context = ssl.create_default_context()
        with smtplib.SMTP_SSL(SMTP_SERVER, 465, timeout=60, context=context) as server:
            server.login(EMAIL, PASSWORD)
            server.sendmail(EMAIL, [RECIPIENT_EMAIL], msg.as_string())  # 同步修改这里
        print("📧 邮件发送成功！")
    except Exception as e:
        print(f"❌ 邮件发送失败: {str(e)}")
        raise


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

                # 获取并处理摘要
                abstract = get_abstract_from_pubmed(paper['pmid'])
                zh_title = translator.translate(pubmed_data['title'])
                zh_abstract = translator.translate(abstract) if abstract else "无可用摘要"

                # 构建数据
                full_data = {
                    **pubmed_data,
                    **paper,
                    'abstract': abstract,
                    'zh_title': zh_title,
                    'zh_abstract': zh_abstract
                }

                # 在生成HTML的部分修改为：
                all_translations.append(f"""
                <div style="margin-bottom: 2rem; padding: 20px; 
                            background: #ffffff; border-radius: 12px;
                            box-shadow: 0 2px 12px rgba(0,0,0,0.1);
                            border: 1px solid #e0e0e0;">
                    <!-- 标题区 -->
                    <div style="border-left: 4px solid #1a73e8; padding-left: 16px; margin-bottom: 20px;">
                        <div style="font-size: 18px; color: #202124; 
                                 line-height: 1.4; font-weight: 600;
                                 margin-bottom: 8px;">
                            {full_data['zh_title']}
                        </div>
                        <div style="font-size: 16px; color: #5f6368;">
                            {full_data['title']}
                        </div>
                    </div>

                    <!-- 元信息 -->
                    <div style="margin-bottom: 20px; padding-left: 20px;">
                        <div style="display: flex; align-items: center; margin-bottom: 12px;">
                            <span style="font-size: 14px; color: #1a73e8; 
                                      margin-right: 8px;">👤</span>
                            <div style="font-size: 14px; color: #5f6368;">
                                {full_data['authors']}
                            </div>
                        </div>
                        <div style="display: flex; align-items: center;">
                            <span style="font-size: 14px; color: #1a73e8; 
                                      margin-right: 8px;">📚</span>
                            <div style="font-size: 16px; color: #5f6368;">
                                {full_data['journal']} ({full_data['year']} IF: {full_data['impact_factor']})
                            </div>
                        </div>
                    </div>

                    <!-- 摘要区块 -->
                    <div style="padding: 16px 20px; background: #f8f9fa; 
                              border-radius: 8px; margin-bottom: 20px;">
                        <div style="display: flex; margin-bottom: 12px;">
                            <div style="width: 3px; height: 20px; 
                                     background: #1a73e8; margin-right: 12px;"></div>
                            <div style="font-size: 15px; color: #202124;
                                     font-weight: 500;">
                                摘要
                            </div>
                        </div>
                        <div style="font-size: 16px; color: #5f6368;
                                 line-height: 1.8; text-align: justify;
                                 white-space: normal;">
                            {full_data['zh_abstract']}
                        </div>
                    </div>

                    <!-- 操作按钮 -->
                    <div style="display: flex; gap: 12px; padding-left: 20px;">
                        <a href="https://pubmed.ncbi.nlm.nih.gov/{full_data['pmid']}" 
                           target="_blank"
                           style="display: inline-flex; align-items: center;
                                  padding: 8px 16px; background: #e8f0fe;
                                  color: #1a73e8; border-radius: 20px;
                                  text-decoration: none; font-size: 14px;
                                  border: 1px solid #cbd5e1;">
                            <span style="margin-right: 6px;">📖</span>
                            PubMed链接
                        </a>
                        {f'<a href="https://doi.org/{full_data["doi"]}" ...>📚 杂志链接</a>' if full_data['doi'] else ''}
                    </div>
                </div>
                """)

                # 邮件主体样式调整
                html_content = f"""
                <html>
                    <body style="font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
                                background: #f1f3f4; margin: 0; padding: 20px;">
                        <div style="max-width: 600px; margin: 0 auto;">
                            <div style="text-align: center; margin-bottom: 30px;">
                                <h1 style="color: #202124; font-size: 22px; 
                                         margin: 24px 0 16px 0;">
                                    📰 每日论文摘要
                                </h1>
                                <div style="font-size: 14px; color: #5f6368;">
                                    {time.strftime("%Y-%m-%d")} 更新
                                </div>
                            </div>
                            {"".join(all_translations)}
                        </div>
                    </body>
                </html>
                """

        if all_translations:
            html_content = f"""
            <html>
                <body style="font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, sans-serif; 
                            max-width: 800px; margin: 0 auto; padding: 2rem 1rem; background-color: #f7fafc;">
                    <header style="text-align: center; margin-bottom: 2.5rem;">
                        <h1 style="color: #1a365d; margin: 0 0 0.5rem 0; font-size: 1.8rem;">
                            📰 今日文献推送 ({len(all_translations)}篇)
                        </h1>
                    </header>
                    {"".join(all_translations)}
                    <footer style="margin-top: 3rem; text-align: center; color: #718096; 
                             font-size: 0.9rem; padding-top: 1.5rem; border-top: 1px solid #e2e8f0;">
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


# 其他辅助函数保持原样...

if __name__ == "__main__":
    main()