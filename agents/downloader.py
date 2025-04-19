import requests, os, hashlib, time, subprocess, json, openai, pathlib, re
from playwright.sync_api import sync_playwright
from config import settings

def oa_pdf(doi: str) -> str | None:
    u = f"https://api.unpaywall.org/v2/{doi}?email={settings.ncbi_email}"
    data = requests.get(u, timeout=10).json()
    return data.get("best_oa_location", {}).get("url_for_pdf")

def crossref_pdf(doi: str) -> str | None:
    r = requests.get(f"https://api.crossref.org/works/{doi}").json()["message"]
    for l in r.get("link", []):
        if l["content-type"] == "application/pdf":
            return l["URL"]
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
