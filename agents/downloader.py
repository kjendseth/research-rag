import requests, os, hashlib, time, subprocess, json, openai, pathlib, re, logging
from playwright.sync_api import sync_playwright
from config import settings

def oa_pdf(doi: str) -> str | None:
    """
    Safely query Unpaywall for the best OA PDF URL.
    Returns None on HTTP error, JSON error, or missing field.
    """
    url = f"https://api.unpaywall.org/v2/{doi}?email={settings.ncbi_email}"
    try:
        resp = requests.get(url, timeout=15)
    except requests.RequestException as ex:
        logging.warning("Unpaywall request error for DOI %s: %s", doi, ex)
        return None

    if resp.status_code != 200:
        logging.warning("Unpaywall returned %d for DOI %s", resp.status_code, doi)
        return None

    # parse JSON safely
    try:
        payload = resp.json()
    except ValueError as ex:
        logging.warning("Invalid JSON from Unpaywall for DOI %s: %s", doi, ex)
        return None

    if not isinstance(payload, dict):
        logging.warning("Unexpected payload type from Unpaywall for DOI %s", doi)
        return None

    # extract PDF URL
    best = payload.get("best_oa_location")
    if not isinstance(best, dict):
        return None

    return best.get("url_for_pdf")
        
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
    """
    Try OA, then Crossref, then optional browser fetch. Safely skip if no URL.
    """
    dest = pathlib.Path(dest_dir)
    dest.mkdir(exist_ok=True)

    # Build a safe filename string, then combine with Path
    filename = re.sub(r"[^\w\.-]", "_", doi) + ".pdf"
    fn = dest / filename

    # If already downloaded, return immediately
    if fn.exists():
        return str(fn)

    # 1) Open Access
    url = oa_pdf(doi)
    # 2) Crossref fallback
    if not url:
        url = crossref_pdf(doi)
    # 3) Browser fallback (omitted here for brevity)...

    if not url:
        logging.info("No PDF URL found for DOI %s", doi)
        return None

    # Download the PDF
    try:
        resp = requests.get(url, timeout=30, stream=True)
        if resp.status_code == 200:
            with open(fn, "wb") as f:
                for chunk in resp.iter_content(65536):
                    f.write(chunk)
            return str(fn)
        else:
            logging.warning("PDF download returned %d for DOI %s", resp.status_code, doi)
    except requests.RequestException as ex:
        logging.warning("Error downloading PDF for DOI %s: %s", doi, ex)

    return None
