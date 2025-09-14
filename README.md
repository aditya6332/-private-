PubMed Knowledge Base Builder (Azure Serverless Edition) This project
provides a robust, secure, and scalable serverless architecture for
building a large-scale knowledge base by scraping PubMed. It is designed
for easy deployment on Microsoft Azure using modern, cloud-native
services.

This version is a complete refactor of the original project, simplifying
the architecture and removing the need for local database tunneling and
manual sync scripts.

Architecture Overview The system is composed of four key cloud
components that work together seamlessly. This design is highly
efficient, scalable, and cost-effective, running almost entirely within
Azure\'s free service tiers for typical use.

Azure Functions (Serverless Compute):

A Python-based Function App acts as the brain of the operation.

Timer Trigger: The main scraper function runs automatically on a
schedule (e.g., every hour). It scrapes PubMed, expands keywords, and
saves results.

HTTP Triggers: Other functions provide APIs for the dashboard to connect
to SignalR, get the initial list of articles, and serve the dashboard
HTML itself.

Azure Table Storage (NoSQL Database):

A highly scalable and extremely low-cost NoSQL database stores all the
scraped article data.

It\'s perfect for this use case as it handles simple, structured data
without the overhead of a relational database.

Azure SignalR Service (Real-time Messaging):

This managed service provides a real-time communication channel between
the scraper function and the web dashboard.

When the scraper finds a new article, it sends a message via SignalR,
which is immediately pushed to all connected dashboard clients, creating
a live feed.

Static HTML Dashboard (Frontend):

A single, dependency-free HTML file acts as the user interface. It\'s
served directly by one of the Azure Functions.

It connects to SignalR to receive live status updates and new articles,
and fetches the most recent articles on initial load.

Deployment Guide Follow these steps to deploy the entire application to
your Azure account.

Prerequisites An Azure account (a free account is sufficient to start).

Azure CLI installed on your machine. (Installation Guide)

Azure Functions Core Tools installed. (npm install -g
azure-functions-core-tools@4)

Python 3.9+ and Git installed locally.

Step 1: Set Up Azure Resources First, log in to Azure and create a
resource group to hold all your services. Then, create the necessary
resources using the Azure CLI.

\# 1. Log in to your Azure account az login

\# 2. Create a resource group to contain all services az group create
\--name PubmedScraperRG \--location \"EastUS\"

\# 3. Create a general-purpose storage account for Table Storage and the
Function App az storage account create \--name pubmedstoragedata\$RANDOM
\--location \"EastUS\" \--resource-group PubmedScraperRG \--sku
Standard_LRS

\# 4. Create a SignalR Service instance (Free Tier) az signalr create
\--name pubmed-scraper-signalr \--resource-group PubmedScraperRG \--sku
Free_F1

\# 5. Create the serverless Function App az functionapp create
\--resource-group PubmedScraperRG \--consumption-plan-location
\"EastUS\" \--runtime python \--runtime-version 3.9 \--functions-version
4 \--name PubmedScraperApp\$RANDOM \--storage-account
\<your-new-storage-account-name\>

Note: Replace \<your-new-storage-account-name\> with the name of the
storage account you just created.

Step 2: Configure the Function App Your Function App needs connection
strings and settings to communicate with other services. These are
stored securely as Application Settings.

Go to the Azure Portal and navigate to your Function App.

Go to Settings -\> Configuration.

Click + New application setting for each of the following:

Setting Name

Value

Description

NCBI_EMAIL

Your actual email address.

Required by the PubMed API for identification.

SEED_KEYWORDS

neuroscience, cognitive psychology

Comma-separated list of keywords to scrape.

TABLE_STORAGE_CONN_STR

Connection string from your Storage Account

Found in your Storage Account under Access keys.

AzureSignalRConnectionString

Connection string from your SignalR Service

Found in your SignalR Service under Keys.

Click Save after adding all the settings.

Step 3: Deploy the Code With the cloud infrastructure in place, you can
now deploy the application code.

\# 1. Clone your repository containing the function app code \# git
clone \...

\# 2. Navigate into the root directory of the project \# cd
your-project-folder

\# 3. Deploy the function app (replace with your Function App\'s name)
func azure functionapp publish \<your-function-app-name\>

Step 4: Access Your Dashboard The application is now live. The scraper
will trigger automatically on its schedule (once per hour).

Find your dashboard URL in the Azure Portal by navigating to your
Function App, clicking on the serve_frontend function, and clicking Get
Function Url.

The URL will look something like:
https://\<your-app-name\>.azurewebsites.net/api/dashboard

Local Development You can run the entire application on your local
machine for testing before deploying.

Fill in local.settings.json: Use the provided local.settings.json file.
Copy the connection strings from your created Azure resources into this
file.

Start the Function App: Open a terminal in the project root and run:

func start

The functions will be available at http://localhost:7071. You can access
the dashboard at http://localhost:7071/api/dashboard.
