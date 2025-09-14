import csv
import time
import configparser
from pathlib import Path
from Bio import Entrez

def load_config(config_file='config.ini'):
    """Loads settings from the configuration file."""
    config = configparser.ConfigParser()
    config.read(config_file)
    return config

def parse_article_details(article_data):
    """
    Parses the complex dictionary from Entrez into a flat dictionary with all requested fields.
    """
    try:
        citation = article_data['MedlineCitation']
        article_info = citation['Article']
        
        pmid = str(citation['PMID'])
        title = article_info.get('ArticleTitle', 'No Title Found')

        # Abstract
        abstract = ""
        if 'Abstract' in article_info and 'AbstractText' in article_info['Abstract']:
            abstract_parts = article_info['Abstract']['AbstractText']
            abstract = ' '.join(abstract_parts) if isinstance(abstract_parts, list) else abstract_parts

        # Authors
        authors_list = []
        if 'AuthorList' in article_info:
            for author in article_info['AuthorList']:
                lastname = author.get('LastName', '')
                forename = author.get('ForeName', '')
                if lastname and forename:
                    authors_list.append(f"{forename} {lastname}")
        authors = ", ".join(authors_list)
        
        # Journal and Publication Date
        journal_info = article_info.get('Journal', {})
        journal = journal_info.get('Title', 'N/A')
        pub_date = journal_info.get('JournalIssue', {}).get('PubDate', {})
        
        pub_year = pub_date.get('Year', pub_date.get('MedlineDate', 'N/A').split(' ')[0])
        pub_month = pub_date.get('Month', 'N/A')
        pub_day = pub_date.get('Day', 'N/A')
        
        publication_date = f"{pub_year}-{pub_month}-{pub_day}"

        # DOI
        doi = ""
        if 'ELocationID' in article_info:
            for elocation in article_info['ELocationID']:
                if elocation.attributes.get('EIdType') == 'doi' and elocation.attributes.get('ValidYN') == 'Y':
                    doi = str(elocation)
                    break
        
        # Keywords (MeSH Terms)
        mesh_terms_list = []
        if 'MeshHeadingList' in citation:
            for heading in citation['MeshHeadingList']:
                mesh_terms_list.append(heading['DescriptorName'])
        
        # PMC ID (for checking free access)
        pmc_id = ""
        if 'OtherID' in citation:
            for other_id in citation['OtherID']:
               if str(other_id).startswith('PMC'):
                   pmc_id = str(other_id)
                   break

        return {
            "pmid": pmid,
            "doi": doi,
            "title": title,
            "abstract": abstract,
            "authors": authors,
            "journal": journal,
            "publication_date": publication_date,
            "publication_year": int(pub_year) if pub_year.isdigit() else None,
            "keywords_mesh": ", ".join(mesh_terms_list),
            "pmc_id": pmc_id
        }
    except (KeyError, IndexError) as e:
        print(f"Skipping an article due to parsing error: {e}")
        return None

def save_articles_to_csv(filename, articles, saved_pmids_set):
    """Saves a list of article dictionaries to a CSV file, avoiding duplicates."""
    if not articles:
        return 0

    file_exists = Path(filename).exists()
    
    new_articles_to_write = []
    for article in articles:
        if article and article['pmid'] not in saved_pmids_set:
            new_articles_to_write.append(article)
            saved_pmids_set.add(article['pmid'])

    if not new_articles_to_write:
        return 0

    fieldnames = new_articles_to_write[0].keys()
    
    try:
        with open(filename, 'a', newline='', encoding='utf-8') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            
            if not file_exists or csvfile.tell() == 0:
                writer.writeheader()
            
            writer.writerows(new_articles_to_write)
    except IOError as e:
        print(f"Error writing to CSV file '{filename}': {e}")
        return 0
        
    return len(new_articles_to_write)

def search_and_expand_pubmed(seed_keyword, email, max_initial, max_final):
    """
    Searches PubMed, expands with MeSH terms, and returns detailed article data.
    """
    Entrez.email = email
    print(f"🔬 Step 1: Performing initial search for '{seed_keyword}'...")
    
    try:
        handle = Entrez.esearch(db="pubmed", term=seed_keyword, retmax=max_initial, sort="relevance")
        record = Entrez.read(handle)
        handle.close()
    except Exception as e:
        print(f"A network error occurred during initial search: {e}. Skipping keyword.")
        return []
        
    id_list_initial = record["IdList"]

    if not id_list_initial:
        print(f"No articles found for '{seed_keyword}'.")
        return []

    print(f"Found {len(id_list_initial)} initial papers. Extracting MeSH terms...")
    time.sleep(0.5) # Pause before next request
    try:
        handle = Entrez.efetch(db="pubmed", id=id_list_initial, rettype="medline", retmode="xml")
        records = Entrez.read(handle)['PubmedArticle']
        handle.close()
    except Exception as e:
        print(f"A network error occurred fetching initial details: {e}. Skipping keyword.")
        return []
    
    mesh_terms = set()
    for article in records:
        details = parse_article_details(article)
        if details and details['keywords_mesh']:
            for term in details['keywords_mesh'].split(', '):
                mesh_terms.add(term)

    if not mesh_terms:
        expanded_keyword = seed_keyword
        print("Could not find MeSH terms for expansion. Using original keyword.")
    else:
        print(f"💡 Step 2: Found {len(mesh_terms)} unique MeSH terms for expansion.")
        mesh_query = " OR ".join([f'"{term}"[MeSH Terms]' for term in mesh_terms])
        expanded_keyword = f'("{seed_keyword}") OR ({mesh_query})'
    
    print("🔍 Step 3: Performing expanded search...")
    time.sleep(0.5)
    try:
        handle = Entrez.esearch(db="pubmed", term=expanded_keyword, retmax=max_final, sort="relevance")
        record = Entrez.read(handle)
        handle.close()
    except Exception as e:
        print(f"A network error occurred during expanded search: {e}. Skipping keyword.")
        return []

    id_list_final = record["IdList"]

    if not id_list_final:
        print("No articles found in the expanded search.")
        return []

    print(f"Found {len(id_list_final)} final papers. Fetching all metadata...")
    time.sleep(0.5)
    try:
        handle = Entrez.efetch(db="pubmed", id=id_list_final, rettype="medline", retmode="xml")
        final_records = Entrez.read(handle)['PubmedArticle']
        handle.close()
    except Exception as e:
        print(f"A network error occurred fetching final details: {e}. Skipping keyword.")
        return []

    all_articles_data = [parse_article_details(rec) for rec in final_records]
    return all_articles_data

def main():
    """Main function to run the scraper."""
    config_file = 'config.ini'
    if not Path(config_file).exists():
        print(f"Error: Configuration file '{config_file}' not found.")
        return

    config = load_config(config_file)
    
    email = config.get('NCBI_API', 'email')
    if email == 'your.email@example.com':
        print("⚠️ Please update your email address in config.ini before running.")
        return
        
    csv_filename = config.get('Output', 'csv_filename')
    keywords = [kw.strip() for kw in config.get('Search_Parameters', 'seed_keywords').split(',') if kw.strip()]
    max_initial = config.getint('Search_Parameters', 'max_papers_initial')
    max_final = config.getint('Search_Parameters', 'max_papers_final')
    
    saved_pmids = set()
    try:
        if Path(csv_filename).exists():
            with open(csv_filename, 'r', newline='', encoding='utf-8') as f:
                reader = csv.DictReader(f)
                for row in reader:
                    if 'pmid' in row:
                        saved_pmids.add(row['pmid'])
            print(f"Loaded {len(saved_pmids)} existing PMIDs from '{csv_filename}'. Will not add duplicates.")
    except Exception as e:
        print(f"Could not read existing CSV file. Starting fresh. Error: {e}")

    total_saved = 0
    for keyword in keywords:
        print(f"\n{'='*50}\nProcessing keyword: {keyword}\n{'='*50}")
        articles = search_and_expand_pubmed(keyword, email, max_initial, max_final)
        
        if articles:
            print(f"Found {len(articles)} potential articles. Saving new ones to '{csv_filename}'...")
            saved_count = save_articles_to_csv(csv_filename, articles, saved_pmids)
            total_saved += saved_count
            print(f"✅ Successfully saved {saved_count} new articles for '{keyword}'.")
        
        time.sleep(1)

    print(f"\n{'='*50}\nScraping complete. A total of {total_saved} new articles were added to your knowledge base.\n{'='*50}")

if __name__ == "__main__":
    main()
