import os
import sys

import psycopg2
from psycopg2 import extras

def get_env_variable(var_name):
    """Gets an environment variable or raises an error if it's not set."""
    value = os.getenv(var_name)
    if value is None:
        raise ValueError(f"Error: Environment variable '{var_name}' is not set. Please create a .env file or set it manually.")
    return value

def connect_to_db(db_config):
    """Establishes a connection to a PostgreSQL database from a dict."""
    try:
        conn = psycopg2.connect(**db_config)
        print(f"Successfully connected to {db_config['host']}.")
        return conn
    except psycopg2.OperationalError as e:
        print(f"❌ Failed to connect to {db_config['host']}: {e}")
        return None

def init_db_schema(conn):
    """Ensures the articles table exists in the local DB."""
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

def main():
    """Main function to run the sync process."""
    try:
        from dotenv import load_dotenv
        load_dotenv()
        print("Loaded credentials from .env file.")
    except ImportError:
        print("Warning: dotenv library not found. Assuming environment variables are set manually.")

    try:
        local_db_conf = {'host': get_env_variable('LOCAL_DB_HOST_SYNC'), 'port': get_env_variable('LOCAL_DB_PORT_SYNC'),
                         'dbname': get_env_variable('LOCAL_DB_NAME'), 'user': get_env_variable('LOCAL_DB_USER'), 'password': get_env_variable('LOCAL_DB_PASS')}
        cloud_db_conf = {'host': get_env_variable('CLOUD_DB_HOST'), 'port': get_env_variable('CLOUD_DB_PORT'), 'dbname': get_env_variable('CLOUD_DB_NAME'),
                         'user': get_env_variable('CLOUD_DB_USER'), 'password': get_env_variable('CLOUD_DB_PASS')}
    except ValueError as e:
        print(e)
        sys.exit(1)

    local_conn = connect_to_db(local_db_conf)
    cloud_conn = connect_to_db(cloud_db_conf)

    if not local_conn or not cloud_conn:
        print("Could not connect to both databases. Aborting sync.")
        return

    try:
        print("Initializing local DB schema...")
        init_db_schema(local_conn)

        with cloud_conn.cursor(cursor_factory=extras.DictCursor) as cloud_cursor, local_conn.cursor() as local_cursor:
            print("Fetching records from the cloud buffer database...")
            cloud_cursor.execute("SELECT * FROM articles;")
            buffered_articles = cloud_cursor.fetchall()

            if not buffered_articles:
                print("✅ No new articles in the cloud buffer. Local DB is up to date.")
                return

            print(f"Found {len(buffered_articles)} articles to sync.")
            insert_query = """
                INSERT INTO articles (pmid, doi, title, abstract, authors, journal,
                    publication_date, publication_year, keywords_mesh, pmc_id, scraped_at)
                VALUES (%(pmid)s, %(doi)s, %(title)s, %(abstract)s, %(authors)s, %(journal)s,
                    %(publication_date)s, %(publication_year)s, %(keywords_mesh)s, %(pmc_id)s, %(scraped_at)s)
                ON CONFLICT (pmid) DO NOTHING;
            """
            extras.execute_batch(local_cursor, insert_query, buffered_articles)
            print(f"Synced {local_cursor.rowcount} new articles to your local database.")
            local_conn.commit()

            print("Clearing synced records from the cloud buffer...")
            pmids_to_delete = tuple(article['pmid'] for article in buffered_articles)
            if pmids_to_delete:
                cloud_cursor.execute("DELETE FROM articles WHERE pmid IN %s;", (pmids_to_delete,))
                cloud_conn.commit()
                print(f"Cleared {cloud_cursor.rowcount} articles from the buffer.")

    except Exception as e:
        print(f"An error occurred during sync: {e}")
        local_conn.rollback()
        cloud_conn.rollback()
    finally:
        local_conn.close()
        cloud_conn.close()

if __name__ == "__main__":
    main()

