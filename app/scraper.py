# pubmed_scraper/app/scraper.py
import os
import time
import logging
import yaml
import xml.etree.ElementTree as ET
from typing import List, Dict
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from app.db import save_papers

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("pubmed_scraper.scraper")

CONFIG_PATH = os.getenv("CONFIG_PATH", "config.yaml")
with open(CONFIG_PATH, "r") as f:
    CONFIG = yaml.safe_load(f)

if os.getenv("PUBMED_EMAIL"):
    CONFIG["pubmed"]["email"] = os.getenv("PUBMED_EMAIL")
if os.getenv("ITERATIONS"):
    CONFIG["iterations"] = int(os.getenv("ITERATIONS"))

BASE_URL = CONFIG["pubmed"]["base_url"]
REQUEST_TIMEOUT = CONFIG["runtime"]["request_timeout"]
RETMAX = CONFIG["pubmed"]["retmax"]

session = requests.Session()
retries = Retry(total=5, backoff_factor=1, status_forcelist=[429, 500, 502, 503, 504])
session.mount("https://", HTTPAdapter(max_retries=retries))

def fetch_pubmed_ids(query: str) -> List[str]:
    params = {
        "db": "pubmed", "term": query, "retmax": RETMAX,
        "retmode": "xml", "email": CONFIG["pubmed"]["email"],
    }
    response = session.get(f"{BASE_URL}esearch.fcgi", params=params, timeout=REQUEST_TIMEOUT)
    response.raise_for_status()
    root = ET.fromstring(response.text)
    return [elem.text for elem in root.findall(".//Id")]

def fetch_metadata(pubmed_ids: List[str]) -> List[Dict]:
    if not pubmed_ids: return []
    params = {
        "db": "pubmed", "id": ",".join(pubmed_ids),
        "retmode": "xml", "email": CONFIG["pubmed"]["email"],
    }
    response = session.get(f"{BASE_URL}efetch.fcgi", params=params, timeout=REQUEST_TIMEOUT)
    response.raise_for_status()
    root = ET.fromstring(response.text)
    papers = []
    for article in root.findall(".//PubmedArticle"):
        pmid = article.findtext(".//PMID")
        title = article.findtext(".//ArticleTitle") or ""
        keywords = [kw.text for kw in article.findall(".//Keyword") if kw.text]
        papers.append({"pmid": pmid, "title": title, "keywords": keywords})
    return papers

def iterative_scrape() -> Dict:
    seed_keywords = set(CONFIG.get("seed_keywords", []))
    iterations = int(CONFIG.get("iterations", 3))
    seen_queries = seed_keywords
    all_papers_count = 0
    for i in range(iterations):
        new_keywords = set()
        for query in list(seen_queries):
            ids = fetch_pubmed_ids(query)
            if not ids: continue
            papers = fetch_metadata(ids)
            if papers:
                save_papers(papers)
                all_papers_count += len(papers)
                for paper in papers:
                    for kw in paper.get("keywords", []):
                        if kw.lower() not in seen_queries:
                            new_keywords.add(kw.lower())
            time.sleep(0.5)
        if not new_keywords: break
        seen_queries.update(new_keywords)
    return {"total_papers_found": all_papers_count, "total_queries": len(seen_queries)}