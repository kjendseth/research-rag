import requests, os, hashlib, time, subprocess, json, openai, pathlib, re, logging
from playwright.sync_api import sync_playwright
from config import settings

def oa_pdf(doi: str) -> str | None:
    u = f"https://api.unpaywall.org/v2/{doi}?email={settings.ncbi_email}"
    data = requests.get(u, timeout=10).json()
    return data.get("best_oa_location", {}).get("url_for_pdf")

def crossref_pdf(doi: str) -> str | None:
    """
    Safely fetch a PDF URL from Crossref's /works/{doi} endpoint.
    Returns the first link with content-type application/pdf, or None.
    """
    url = f"https://api.crossref.org/works/{doi}"
    headers = {
        "User-Agent": f"research-rag/0.1 (mailto:{os.getenv('NCBI_EMAIL','user@example.com')})",
        "Accept": "application/json"
    }
    try:
        resp = requests.get(url, headers=headers, timeout=15)
        if resp.status_code != 200:
            logging.warning("Crossref returned %d for DOI %s", resp.status_code, doi)
            return None
        payload = resp.json()
        msg = payload.get("message", {})
        for link in msg.get("link", []):
            if link.get("content-type") == "application/pdf":
                return link.get("URL")
    except requests.RequestException as ex:
        logging.warning("Crossref request exception for DOI %s: %s", doi, ex)
    except ValueError as ex:
        logging.warning("Invalid JSON from Crossref for DOI %s: %s", doi, ex)
    return None
    
def fetch_pdf(doi: str, dest_dir="pdfs") -> str | None:
    pathlib.Path(dest_dir).mkdir(exist_ok=True)
    url = oa_pdf(doi) or crossref_pdf(doi)
    if not url and settings.use_headless_browser:
        url = settings.ezproxy_prefix + f"https://doi.org/{doi}"
        with sync_playwright() as p:
            b = p.chromium.launch(headless=True)
            page = b.new_page()
            page.goto(url)
            links = page.locator("a[href$='.pdf']")
            if links.count() > 0:
                url = links.first.get_attribute("href")
            b.close()
    if not url:
        return None
    fn = dest_dir + "/" + re.sub(r'[^\w.-]', '_', doi) + ".pdf"
    if os.path.exists(fn):
        return fn
    r = requests.get(url, timeout=30, stream=True)
    if r.status_code == 200:
        with open(fn, "wb") as f:
            for c in r.iter_content(65536):
                f.write(c)
    return fn if os.path.exists(fn) else None
