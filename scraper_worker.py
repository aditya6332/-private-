import configparser
import os
import smtplib
import socket
import threading
import time
from email.mime.text import MIMEText

import psycopg2
from Bio import Entrez
from psycopg2 import pool

# --- Global State for Email Reporting ---
articles_processed_since_last_report = 0
current_keyword_being_processed = "None"

def get_env_variable(var_name):
    """Gets an environment variable or raises an error if it's not set."""
    value = os.getenv(var_name)
    if value is None:
        raise ValueError(f"Error: Environment variable '{var_name}' is not set.")
    return value

def load_config():
    """Loads non-sensitive settings from config.ini."""
    config = configparser.ConfigParser()
    config.read('config.ini')
    return config

def create_db_pool(db_config):
    """Creates a thread-safe PostgreSQL connection pool from a dict of credentials."""
    try:
        return psycopg2.pool.SimpleConnectionPool(
            1, 5,
            host=db_config['host'],
            port=db_config['port'],
            dbname=db_config['dbname'],
            user=db_config['user'],
            password=db_config['password']
        )
    except (psycopg2.OperationalError, ValueError) as e:
        print(f"Error creating DB pool for {db_config.get('host', 'N/A')}: {e}")
        return None

def check_local_db_availability(local_db_config):
    """Checks if the local DB server (via ngrok tunnel) is reachable."""
    if not all(k in local_db_config and local_db_config[k] for k in ['host', 'port']):
        return False
    try:
        with socket.create_connection((local_db_config['host'], int(local_db_config['port'])), timeout=5):
            print("✅ Local DB (via ngrok) is reachable.")
            return True
    except (socket.timeout, socket.gaierror, ConnectionRefusedError, OSError) as e:
        print(f"⚠️ Local DB (via ngrok) is not reachable: {e}. Using cloud fallback.")
        return False

def init_db_schema(conn):
    """Ensures the articles table exists with the correct schema."""
    with conn.cursor() as cursor:
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS articles (
                pmid TEXT PRIMARY KEY, doi TEXT, title TEXT, abstract TEXT, authors TEXT,
                journal TEXT, publication_date TEXT, publication_year INTEGER,
                keywords_mesh TEXT, pmc_id TEXT, scraped_at TIMESTAMPTZ DEFAULT NOW()
            );
            CREATE UNIQUE INDEX IF NOT EXISTS doi_unique_idx ON articles (doi) WHERE doi IS NOT NULL AND doi != '';
        """)
        conn.commit()

def save_articles_to_db(db_pool, articles):
    """Saves a list of articles to the specified database using INSERT ON CONFLICT."""
    global articles_processed_since_last_report
    if not articles or not db_pool: return 0

    conn = db_pool.getconn()
    try:
        init_db_schema(conn)
        with conn.cursor() as cursor:
            insert_query = """
                INSERT INTO articles (pmid, doi, title, abstract, authors, journal,
                    publication_date, publication_year, keywords_mesh, pmc_id)
                VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
                ON CONFLICT (pmid) DO NOTHING;
            """
            data_to_insert = [(a['pmid'], a['doi'], a['title'], a['abstract'], a['authors'], a['journal'],
                               a['publication_date'], a['publication_year'], a['keywords_mesh'], a['pmc_id'])
                              for a in articles if a]
            if not data_to_insert: return 0
            
            # Using execute_values is more efficient for psycopg2 v2.7+
            psycopg2.extras.execute_batch(cursor, insert_query, data_to_insert)
            saved_count = cursor.rowcount
            conn.commit()
            articles_processed_since_last_report += saved_count
            return saved_count
    except Exception as e:
        print(f"DB save error: {e}")
        conn.rollback()
        return 0
    finally:
        db_pool.putconn(conn)

def send_status_email(email_conf, total_articles_scraped):
    """Sends a status update email."""
    global articles_processed_since_last_report, current_keyword_being_processed
    subject = "PubMed Scraper Status Update"
    body = f"""Hello,

This is an automated update from your PubMed scraper running on Koyeb.

-   Current Keyword: {current_keyword_being_processed}
-   New Articles (last 30 mins): {articles_processed_since_last_report}
-   Total Scraped (this session): {total_articles_scraped}

The scraper is running. You can check the logs in your Koyeb dashboard.

Best,
Your Automated Scraper Bot
"""
    msg = MIMEText(body)
    msg['Subject'], msg['From'], msg['To'] = subject, email_conf['sender'], email_conf['recipient']
    try:
        with smtplib.SMTP(email_conf['server'], int(email_conf['port'])) as server:
            server.starttls()
            server.login(email_conf['sender'], email_conf['password'])
            server.send_message(msg)
        print(f"📬 Status email sent to {email_conf['recipient']}.")
        articles_processed_since_last_report = 0
    except Exception as e:
        print(f"❌ Failed to send status email: {e}")

def parse_article_details(article_data):
    try:
        citation = article_data['MedlineCitation']
        article_info = citation['Article']
        pmid = str(citation['PMID'])
        title = article_info.get('ArticleTitle', 'No Title Found')
        abstract_parts = article_info.get('Abstract', {}).get('AbstractText', [])
        abstract = ' '.join(abstract_parts) if isinstance(abstract_parts, list) else str(abstract_parts or '')
        authors_list = [f"{a.get('ForeName', '')} {a.get('LastName', '')}".strip() for a in article_info.get('AuthorList', [])]
        authors = ", ".join(authors_list)
        journal = article_info.get('Journal', {}).get('Title', 'N/A')
        pub_date = article_info.get('Journal', {}).get('JournalIssue', {}).get('PubDate', {})
        pub_year = pub_date.get('Year', pub_date.get('MedlineDate', 'N/A').split(' ')[0])
        publication_date = f"{pub_year}-{pub_date.get('Month', 'N/A')}-{pub_date.get('Day', 'N/A')}"
        doi = next((str(eloc) for eloc in article_info.get('ELocationID', []) if eloc.attributes.get('EIdType') == 'doi'), "")
        mesh_terms_list = [h['DescriptorName'] for h in citation.get('MeshHeadingList', [])]
        pmc_id = next((str(oid) for oid in citation.get('OtherID', []) if str(oid).startswith('PMC')), "")
        return {"pmid": pmid, "doi": doi, "title": title, "abstract": abstract, "authors": authors, "journal": journal,
                "publication_date": publication_date, "publication_year": int(pub_year) if pub_year.isdigit() else None,
                "keywords_mesh": ", ".join(mesh_terms_list), "pmc_id": pmc_id}
    except (KeyError, IndexError): return None

def search_and_expand_pubmed(seed_keyword, email, max_initial, max_final):
    Entrez.email = email
    print(f"🔬 Initial search for '{seed_keyword}'...")
    try:
        handle = Entrez.esearch(db="pubmed", term=seed_keyword, retmax=max_initial, sort="relevance")
        id_list_initial = Entrez.read(handle)["IdList"]
        handle.close()
        if not id_list_initial: return []
        time.sleep(0.4)
        handle = Entrez.efetch(db="pubmed", id=id_list_initial, rettype="medline", retmode="xml")
        records = Entrez.read(handle)['PubmedArticle']
        handle.close()
        mesh_terms = {term for article in records for term in (parse_article_details(article) or {}).get('keywords_mesh', '').split(', ') if term}
        expanded_keyword = f'("{seed_keyword}") OR ({" OR ".join([f""""{term}"[MeSH Terms]""" for term in mesh_terms])})' if mesh_terms else seed_keyword
        print(f"🔍 Expanded search with {len(mesh_terms)} extra terms...")
        time.sleep(0.4)
        handle = Entrez.esearch(db="pubmed", term=expanded_keyword, retmax=max_final, sort="relevance")
        id_list_final = Entrez.read(handle)["IdList"]
        handle.close()
        if not id_list_final: return []
        print(f"Fetching {len(id_list_final)} final papers...")
        time.sleep(0.4)
        handle = Entrez.efetch(db="pubmed", id=id_list_final, rettype="medline", retmode="xml")
        final_records = Entrez.read(handle)['PubmedArticle']
        handle.close()
        return [parse_article_details(rec) for rec in final_records]
    except Exception as e:
        print(f"An error occurred during PubMed search: {e}")
        return []

def periodic_email_scheduler(email_conf, stop_event, get_total_count_func):
    """Thread function to send emails every 30 minutes."""
    while not stop_event.is_set():
        time.sleep(1800)
        if not stop_event.is_set():
            send_status_email(email_conf, get_total_count_func())

def main():
    global current_keyword_being_processed
    try:
        ncbi_email = get_env_variable('NCBI_EMAIL')
        email_conf = {'server': get_env_variable('SMTP_SERVER'), 'port': get_env_variable('SMTP_PORT'), 'sender': get_env_variable('SENDER_EMAIL'),
                      'password': get_env_variable('SENDER_PASSWORD'), 'recipient': get_env_variable('RECIPIENT_EMAIL')}
        cloud_db_conf = {'host': get_env_variable('CLOUD_DB_HOST'), 'port': get_env_variable('CLOUD_DB_PORT'), 'dbname': get_env_variable('CLOUD_DB_NAME'),
                         'user': get_env_variable('CLOUD_DB_USER'), 'password': get_env_variable('CLOUD_DB_PASS')}
        local_db_conf = {'host': os.getenv('LOCAL_DB_HOST'), 'port': os.getenv('LOCAL_DB_PORT'), 'dbname': os.getenv('LOCAL_DB_NAME'),
                         'user': os.getenv('LOCAL_DB_USER'), 'password': os.getenv('LOCAL_DB_PASS')}
    except ValueError as e:
        print(e)
        return

    config = load_config()
    cloud_db_pool = create_db_pool(cloud_db_conf)
    if not cloud_db_pool:
        print("❌ Critical error: Could not connect to the cloud fallback database. Aborting.")
        return

    total_scraped = 0
    stop_event = threading.Event()
    email_thread = threading.Thread(target=periodic_email_scheduler, args=(email_conf, stop_event, lambda: total_scraped), daemon=True)
    email_thread.start()

    try:
        keywords = [kw.strip() for kw in config.get('Search_Parameters', 'seed_keywords').split(',') if kw.strip()]
        while True:
            for keyword in keywords:
                current_keyword_being_processed = keyword
                print(f"\n{'='*50}\nProcessing: {keyword}\n{'='*50}")
                articles = search_and_expand_pubmed(keyword, ncbi_email, config.getint('Search_Parameters', 'max_papers_initial'), config.getint('Search_Parameters', 'max_papers_final'))
                
                if articles:
                    target_pool = None
                    local_pool = None
                    if check_local_db_availability(local_db_conf):
                        local_pool = create_db_pool(local_db_conf)
                        if local_pool: target_pool = local_pool
                    
                    if not target_pool: target_pool = cloud_db_pool
                    
                    saved_count = save_articles_to_db(target_pool, articles)
                    total_scraped += saved_count
                    print(f"✅ Saved {saved_count} new articles for '{keyword}'.")
                    
                    if local_pool: local_pool.closeall()
                time.sleep(1)
            print("\nCycle complete. Restarting in 60 minutes...")
            time.sleep(3600)
    except KeyboardInterrupt:
        print("\n🛑 Shutdown signal.")
    finally:
        stop_event.set()
        email_thread.join(2)
        if cloud_db_pool: cloud_db_pool.closeall()
        print("Shutdown complete.")

if __name__ == "__main__":
    main()

