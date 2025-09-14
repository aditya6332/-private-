PubMed Distributed Knowledge Base Builder (Koyeb Edition)
This project provides a robust, secure, and scalable architecture for building a large-scale knowledge base by scraping PubMed. It is designed for deployment on the Koyeb serverless platform, ensuring high availability for the scraper and enhanced security for your data and credentials.

Architecture & Security
This system leverages a modern, secure architecture:

Koyeb Git-Driven App: The scraper runs as a persistent background worker on Koyeb, deployed directly from a GitHub repository.

Koyeb Managed PostgreSQL DB: A secure, managed database in the cloud acts as a temporary buffer for scraped data.

Secure Secrets: All sensitive credentials (passwords, API keys) are stored as Environment Variables in Koyeb, not in code.

Secure Tunneling with ngrok: Instead of insecure port forwarding, ngrok creates a temporary, encrypted connection to your local database when it's online.

Local Sync Script: A local script pulls data from the cloud buffer to your main knowledge base on your laptop, ensuring your home network is never exposed.

Step-by-Step Deployment Guide
Step 1: Prerequisites
Accounts:

A GitHub account.

A Koyeb account.

An ngrok account (the free tier is sufficient).

Software:

Git installed on your machine.

PostgreSQL installed on your laptop.

Python 3.8+ installed locally.

The ngrok executable downloaded to your laptop.

Step 2: Set Up the Cloud Database on Koyeb
Log in to your Koyeb dashboard.

Navigate to the Databases tab and click Create Database Service.

Choose PostgreSQL, select a region (e.g., Frankfurt, Germany), name your service (e.g., pubmed-buffer-db), and create it.

Once deployed, Koyeb will provide you with all the connection credentials (host, port, database name, user, password). Copy these to a temporary, secure note for the next steps.

Step 3: Prepare and Push Your Code to GitHub
Create a new, private repository on GitHub (e.g., pubmed-scraper).

Save the four project files (scraper_worker.py, sync_local_db.py, config.ini, requirements.txt) into a folder on your computer.

Crucially, add a .gitignore file to prevent committing sensitive local files:

# .gitignore
.env
__pycache__/

Initialize a Git repository and push the code:

git init
git add .
git commit -m "Initial commit"
git branch -M main
git remote add origin [https://github.com/YOUR_USERNAME/pubmed-scraper.git](https://github.com/YOUR_USERNAME/pubmed-scraper.git)
git push -u origin main

Step 4: Deploy the Scraper App on Koyeb
On the Koyeb dashboard, go to the Apps tab and click Create App.

Choose GitHub as the deployment method. Select your pubmed-scraper repository.

Service Type: Change the service type from Web Service to Worker. This is for background jobs that don't need a public URL.

Run Command: Ensure the run command is python scraper_worker.py.

Environment Variables: This is the most critical step for security. Click Add Variable for each of the following, using the "Secret" type for all passwords.

NCBI_EMAIL: Your actual email for the PubMed API.

Cloud DB Credentials (from Koyeb DB):

CLOUD_DB_HOST, CLOUD_DB_PORT, CLOUD_DB_NAME, CLOUD_DB_USER, CLOUD_DB_PASS

Email Notification Credentials:

SMTP_SERVER (e.g., smtp.gmail.com)

SMTP_PORT (e.g., 587)

SENDER_EMAIL (your email)

SENDER_PASSWORD (your email "App Password" if using 2FA)

RECIPIENT_EMAIL (where you want to receive updates)

Local DB Tunnel Credentials (leave these blank for now):

LOCAL_DB_HOST, LOCAL_DB_PORT, LOCAL_DB_NAME, LOCAL_DB_USER, LOCAL_DB_PASS

Click Deploy. Koyeb will pull your code, install dependencies, and start your worker. It will use the cloud DB by default since the local DB variables are empty.

Step 5: Local Laptop Usage
A. Running the Secure Tunnel with ngrok
When you want your laptop to be available for direct data transfer:

Authenticate ngrok: ngrok config add-authtoken YOUR_NGROK_TOKEN

Start a TCP tunnel to your local PostgreSQL port (usually 5432):

ngrok tcp 5432

ngrok will give you a forwarding address, like tcp://2.tcp.ngrok.io:12345.

The Host is 2.tcp.ngrok.io.

The Port is 12345.

Go back to your App's settings in Koyeb. Update the LOCAL_DB_HOST and LOCAL_DB_PORT environment variables with these new values to direct traffic to your laptop. Also, fill in your local LOCAL_DB_NAME, LOCAL_DB_USER, and LOCAL_DB_PASS. The app will automatically redeploy and begin using the tunnel.

B. Running the Sync Script
When ngrok is not running, data accumulates in the cloud buffer. To sync it:

Install local dependencies: pip install -r requirements.txt

Create a .env file in the same directory as sync_local_db.py. This file securely stores your local credentials without putting them in your code.

# .env file for sync script
# Your local DB running on your machine
LOCAL_DB_HOST_SYNC=localhost
LOCAL_DB_PORT_SYNC=5432
LOCAL_DB_NAME=knowledge_base # or whatever you named it
LOCAL_DB_USER=your_local_db_user
LOCAL_DB_PASS=your_local_db_password

# Copy the cloud DB credentials here from Koyeb
CLOUD_DB_HOST=...
CLOUD_DB_PORT=...
CLOUD_DB_NAME=...
CLOUD_DB_USER=...
CLOUD_DB_PASS=...

Run the sync script from your terminal:

python sync_local_db.py

This will securely connect to both databases, transfer the new data, and clear the buffer.