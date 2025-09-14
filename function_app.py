# This file defines the serverless functions that power the application.
# It uses the Azure Functions v2 Python programming model.

import azure.functions as func
from azure.data.tables import TableServiceClient
import os
import logging
import json
from Bio import Entrez
import time
from pathlib import Path

# Helper function to parse article details from PubMed XML
def parse_article_details(article_data):
    """Safely parses the XML from Entrez into a flat dictionary."""
    try:
        citation = article_data['MedlineCitation']
        article_info = citation['Article']
        pmid = str(citation['PMID'])
        title = article_info.get('ArticleTitle', 'No Title Found')
        
        # Abstract
        abstract_parts = article_info.get('Abstract', {}).get('AbstractText', [])
        abstract = ' '.join(abstract_parts) if isinstance(abstract_parts, list) else str(abstract_parts or '')
        
        # Authors
        authors_list = [f"{a.get('ForeName', '')} {a.get('LastName', '')}".strip() for a in article_info.get('AuthorList', [])]
        authors = ", ".join(filter(None, authors_list))
        
        # Journal & Publication
        journal = article_info.get('Journal', {}).get('Title', 'N/A')
        pub_date = article_info.get('Journal', {}).get('JournalIssue', {}).get('PubDate', {})
        pub_year_str = pub_date.get('Year', pub_date.get('MedlineDate', 'N/A').split(' ')[0])
        publication_date = f"{pub_year_str}-{pub_date.get('Month', 'N/A')}-{pub_date.get('Day', 'N/A')}"
        
        # Identifiers
        doi = next((str(eloc) for eloc in article_info.get('ELocationID', []) if eloc.attributes.get('EIdType') == 'doi'), "")
        mesh_terms_list = [h['DescriptorName'] for h in citation.get('MeshHeadingList', [])]
        
        return {
            "PartitionKey": "pubmed", # Required for Table Storage
            "RowKey": pmid,           # Required for Table Storage, unique ID
            "doi": doi,
            "title": title,
            "abstract": abstract[:30000], # Truncate to avoid size limits
            "authors": authors,
            "journal": journal,
            "publication_date": publication_date,
            "publication_year": int(pub_year_str) if pub_year_str.isdigit() else 0,
            "keywords_mesh": ", ".join(mesh_terms_list),
        }
    except (KeyError, IndexError, TypeError) as e:
        logging.error(f"Error parsing article with PMID {citation.get('PMID', 'N/A')}: {e}")
        return None

# --- AZURE FUNCTIONS ---

app = func.FunctionApp()

# 1. Negotiate Function: Provides connection info for SignalR to the web client.
@app.function_name(name="negotiate")
@app.route(route="negotiate", auth_level=func.AuthLevel.ANONYMOUS)
@app.generic_input_binding(arg_name="connectionInfo", type="signalRConnectionInfo", hubName="serverless")
def negotiate(req: func.HttpRequest, connectionInfo: str) -> func.HttpResponse:
    """Provides SignalR connection credentials to the client."""
    return func.HttpResponse(connectionInfo, mimetype='application/json')

# 2. Serve Frontend Function: Serves the main dashboard HTML page.
@app.function_name(name="serve_frontend")
@app.route(route="dashboard", auth_level=func.AuthLevel.ANONYMOUS)
def serve_frontend(req: func.HttpRequest) -> func.HttpResponse:
    """Serves the static index.html file as the dashboard."""
    try:
        # Assumes index.html is in a 'frontend' directory relative to the function app root
        html_path = Path(__file__).parent.parent / 'frontend' / 'index.html'
        with open(html_path, 'r') as f:
            html_content = f.read()
        return func.HttpResponse(html_content, mimetype='text/html')
    except FileNotFoundError:
        return func.HttpResponse("Dashboard file not found.", status_code=404)


# 3. Scraper Function: Runs on a timer to scrape PubMed and update clients.
@app.function_name(name="ScraperTimerTrigger")
@app.timer_trigger(schedule="0 0 * * * *", arg_name="myTimer", run_on_startup=False) # Runs every hour
@app.generic_output_binding(arg_name="signalRMessages", type="signalR", hubName="serverless")
def timer_trigger(myTimer: func.TimerRequest, signalRMessages: func.Out[str]) -> None:
    """Main scraper function, runs on a schedule."""
    if myTimer.past_due:
        logging.info('The timer is past due!')
    logging.info('Python timer trigger function executed.')

    try:
        ENTREZ_EMAIL = os.environ["NCBI_EMAIL"]
        KEYWORDS = [kw.strip() for kw in os.environ["SEED_KEYWORDS"].split(',')]
        TABLE_STORAGE_CONN_STR = os.environ["TABLE_STORAGE_CONN_STR"]
        TABLE_NAME = "pubmedarticles"
        Entrez.email = ENTREZ_EMAIL
    except KeyError as e:
        logging.error(f"Missing required application setting: {e}")
        return
        
    table_service_client = TableServiceClient.from_connection_string(conn_str=TABLE_STORAGE_CONN_STR)
    table_client = table_service_client.get_table_client(table_name=TABLE_NAME)
    try:
        table_client.create_table()
        logging.info(f"Table '{TABLE_NAME}' created.")
    except Exception:
        logging.info(f"Table '{TABLE_NAME}' already exists.")

    for keyword in KEYWORDS:
        status_message = f"Starting search for keyword: {keyword}"
        logging.info(status_message)
        signalRMessages.set(json.dumps({'target': 'statusUpdate', 'arguments': [status_message]}))
        time.sleep(1)

        try:
            handle = Entrez.esearch(db="pubmed", term=keyword, retmax=20, sort="relevance")
            initial_ids = Entrez.read(handle)["IdList"]
            handle.close()
            if not initial_ids:
                logging.warning(f"No initial articles found for '{keyword}'.")
                continue

            time.sleep(0.5)
            handle = Entrez.efetch(db="pubmed", id=initial_ids, rettype="medline", retmode="xml")
            initial_articles = Entrez.read(handle)['PubmedArticle']
            handle.close()

            mesh_terms = {term for article in initial_articles for term in (parse_article_details(article) or {}).get('keywords_mesh', '').split(', ') if term}
            
            expanded_keyword = f'("{keyword}") OR ({" OR ".join([f""""{term}"[MeSH Terms]""" for term in mesh_terms])})'
            status_message = f"Found {len(mesh_terms)} MeSH terms. Performing expanded search for '{keyword}'..."
            logging.info(status_message)
            signalRMessages.set(json.dumps({'target': 'statusUpdate', 'arguments': [status_message]}))
            time.sleep(1)

            handle = Entrez.esearch(db="pubmed", term=expanded_keyword, retmax=100, sort="relevance")
            final_ids = Entrez.read(handle)["IdList"]
            handle.close()
            if not final_ids:
                logging.warning("No articles found in expanded search.")
                continue
            
            time.sleep(0.5)
            handle = Entrez.efetch(db="pubmed", id=final_ids, rettype="medline", retmode="xml")
            final_articles_xml = Entrez.read(handle)['PubmedArticle']
            handle.close()

            new_articles_count = 0
            for article_xml in final_articles_xml:
                article_entity = parse_article_details(article_xml)
                if article_entity:
                    try:
                        table_client.upsert_entity(entity=article_entity, mode='merge')
                        new_articles_count += 1
                        signalRMessages.set(json.dumps({'target': 'newArticle', 'arguments': [article_entity]}))
                    except Exception as e:
                        logging.warning(f"Could not add article {article_entity['RowKey']}. Error: {e}")

            final_status = f"Finished keyword '{keyword}'. Found and saved {new_articles_count} new/updated articles."
            logging.info(final_status)
            signalRMessages.set(json.dumps({'target': 'statusUpdate', 'arguments': [final_status]}))
            time.sleep(1)

        except Exception as e:
            error_message = f"An error occurred during PubMed search for '{keyword}': {e}"
            logging.error(error_message, exc_info=True)
            signalRMessages.set(json.dumps({'target': 'statusUpdate', 'arguments': [error_message]}))

    logging.info("Scraping cycle complete.")
    signalRMessages.set(json.dumps({'target': 'statusUpdate', 'arguments': ["Scraping cycle complete. Waiting for next run."]}))

# 4. Get Initial Data Function: Fetches recent articles for the web page on load.
@app.function_name(name="get_initial_data")
@app.route(route="get-initial-data", auth_level=func.AuthLevel.ANONYMOUS)
def get_initial_data(req: func.HttpRequest) -> func.HttpResponse:
    """Returns the 20 most recently scraped articles."""
    try:
        TABLE_STORAGE_CONN_STR = os.environ["TABLE_STORAGE_CONN_STR"]
        TABLE_NAME = "pubmedarticles"
        table_service_client = TableServiceClient.from_connection_string(conn_str=TABLE_STORAGE_CONN_STR)
        table_client = table_service_client.get_table_client(table_name=TABLE_NAME)

        entities = list(table_client.query_entities(query_filter="", results_per_page=50))
        sorted_entities = sorted(entities, key=lambda x: x.metadata['timestamp'], reverse=True)
        
        results = [{key: value for key, value in entity.items() if not key.startswith('odata')} for entity in sorted_entities[:20]]

        return func.HttpResponse(json.dumps(results), mimetype="application/json")
    except Exception as e:
        logging.error(f"Error getting initial data: {e}")
        return func.HttpResponse("Error fetching initial articles.", status_code=500)

