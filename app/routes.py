# pubmed_scraper/app/routes.py
import logging
from fastapi import FastAPI, Request, BackgroundTasks
from fastapi.responses import HTMLResponse
from fastapi.templating import Jinja2Templates
from app.scraper import iterative_scrape
from app.db import transfer_atlas_to_local

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("pubmed_scraper.routes")

app = FastAPI()
templates = Jinja2Templates(directory="app/templates")

SCRAPE_STATUS = {"status": "idle", "last_result": None}
TRANSFER_STATUS = {"status": "idle", "last_result": None}

def run_scrape_task():
    global SCRAPE_STATUS
    SCRAPE_STATUS["status"] = "running"
    result = iterative_scrape()
    SCRAPE_STATUS["last_result"] = result
    SCRAPE_STATUS["status"] = "completed"

def run_transfer_task():
    global TRANSFER_STATUS
    TRANSFER_STATUS["status"] = "running"
    result = transfer_atlas_to_local()
    TRANSFER_STATUS["last_result"] = result
    TRANSFER_STATUS["status"] = "completed" if result.get("status") == "success" else "error"

@app.get("/", response_class=HTMLResponse)
async def home(request: Request):
    return templates.TemplateResponse("index.html", {"request": request, "scrape_status": SCRAPE_STATUS, "transfer_status": TRANSFER_STATUS})

@app.post("/scrape", response_class=HTMLResponse)
async def start_scrape(request: Request, background_tasks: BackgroundTasks):
    if SCRAPE_STATUS["status"] != "running":
        background_tasks.add_task(run_scrape_task)
        SCRAPE_STATUS["status"] = "running"
    return templates.TemplateResponse("index.html", {"request": request, "scrape_status": SCRAPE_STATUS, "transfer_status": TRANSFER_STATUS})

@app.post("/transfer-data", response_class=HTMLResponse)
async def start_transfer(request: Request, background_tasks: BackgroundTasks):
    if TRANSFER_STATUS["status"] != "running":
        background_tasks.add_task(run_transfer_task)
        TRANSFER_STATUS["status"] = "running"
    return templates.TemplateResponse("index.html", {"request": request, "scrape_status": SCRAPE_STATUS, "transfer_status": TRANSFER_STATUS})

@app.get("/status")
async def get_status():
    return {"scraper": SCRAPE_STATUS, "transfer": TRANSFER_STATUS}