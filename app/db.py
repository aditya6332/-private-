# pubmed_scraper/app/db.py
import os
import logging
import yaml
from pymongo import MongoClient, errors

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("pubmed_scraper.db")

CONFIG_PATH = os.getenv("CONFIG_PATH", "config.yaml")
with open(CONFIG_PATH, "r") as f:
    CONFIG = yaml.safe_load(f)

# Allow environment variables to override config.yaml
CLOUD_URI = os.getenv("MONGODB_CLOUD_URI") or CONFIG["database"].get("cloud_uri")
LOCAL_URI = os.getenv("MONGODB_LOCAL_URI") or CONFIG["database"].get("local_uri")
DB_NAME = os.getenv("MONGODB_DB_NAME") or CONFIG["database"].get("db_name", "pubmed_scraper")
MODE = os.getenv("DB_MODE") or CONFIG["database"].get("mode", "auto")

def try_connect(uri, timeout_ms=3000):
    if not uri:
        return None
    try:
        client = MongoClient(uri, serverSelectionTimeoutMS=timeout_ms)
        client.server_info()
        logger.info(f"Successfully connected to MongoDB at {uri.split('@')[-1]}")
        return client
    except Exception as e:
        logger.warning(f"Failed to connect to MongoDB: {e}")
        return None

def get_db():
    if MODE == "local":
        client = try_connect(LOCAL_URI)
        if not client:
            raise RuntimeError("Local DB mode selected, but connection failed.")
        return client[DB_NAME]
    
    if MODE == "cloud":
        client = try_connect(CLOUD_URI)
        if not client:
            raise RuntimeError("Cloud DB mode selected, but connection failed.")
        return client[DB_NAME]

    client = try_connect(LOCAL_URI)
    if client:
        return client[DB_NAME]
    
    logger.info("Local DB not available, falling back to cloud DB.")
    client = try_connect(CLOUD_URI)
    if client:
        return client[DB_NAME]
        
    raise RuntimeError("Could not connect to any MongoDB instance (local or cloud).")

def save_papers(papers: list):
    if not papers:
        return {"inserted": 0, "skipped": 0}
    db = get_db()
    db.papers.create_index("pmid", unique=True, sparse=True)
    inserted_count = 0
    skipped_count = 0
    for paper in papers:
        filter_doc = {"pmid": paper.get("pmid")}
        if filter_doc["pmid"]:
            try:
                result = db.papers.update_one(filter_doc, {"$setOnInsert": paper}, upsert=True)
                if result.upserted_id:
                    inserted_count += 1
                else:
                    skipped_count += 1
            except errors.DuplicateKeyError:
                skipped_count += 1
    return {"inserted": inserted_count, "skipped": skipped_count}

def transfer_atlas_to_local():
    logger.info("Attempting to transfer data from Atlas to local DB...")
    local_client = try_connect(LOCAL_URI, timeout_ms=2000)
    if not local_client:
        message = "Local database is not active. Cannot start transfer."
        return {"status": "error", "message": message}
    cloud_client = try_connect(CLOUD_URI, timeout_ms=5000)
    if not cloud_client:
        message = "MongoDB Atlas is not reachable. Cannot start transfer."
        return {"status": "error", "message": message}
    try:
        local_db = local_client[DB_NAME]
        cloud_db = cloud_client[DB_NAME]
        local_papers_collection = local_db.papers
        cloud_papers_collection = cloud_db.papers
        local_papers_collection.create_index("pmid", unique=True, sparse=True)
        cloud_papers = list(cloud_papers_collection.find({}))
        if not cloud_papers:
            return {"status": "success", "message": "No documents found in Atlas to transfer.", "transferred": 0}
        transferred_count = 0
        skipped_count = 0
        for paper in cloud_papers:
            filter_doc = {"pmid": paper.get("pmid")}
            if filter_doc["pmid"]:
                result = local_papers_collection.update_one(filter_doc, {"$set": paper}, upsert=True)
                if result.upserted_id or result.modified_count > 0:
                    transferred_count += 1
                else:
                    skipped_count += 1
        message = f"Data transfer complete. Transferred/Updated: {transferred_count}, Skipped: {skipped_count}."
        return {"status": "success", "message": message, "transferred": transferred_count}
    finally:
        if local_client:
            local_client.close()
        if cloud_client:
            cloud_client.close()
            