#!/usr/bin/env python3
"""
Longevity Evidence Scout v1.2
Automated discovery of healthspan and blood biomarker research from PubMed

Based on ToxEcology Evidence Scout v2.2 architecture
Adapted for longevity science and blood biomarkers

v1.2 CHANGES:
- Added auto-linking to Health_Conditions and Symptom_Clusters
- Domain → Condition mapping for proper linking
- Caching system for related tables
- Improved scoring for cross-sectional/NHANES studies
"""

import os
import sys
import json
import time
import yaml
import requests
import xml.etree.ElementTree as ET
from datetime import datetime
from anthropic import Anthropic

# ============================================================================
# DOMAIN MAPPING - Must match toxin_domain single select options exactly
# ============================================================================

DOMAIN_MAP = {
    "Inflammation": "DOM_INFLAMMATION",
    "Cardiovascular": "DOM_LIPID",
    "Kidney_Liver": "DOM_KIDNEY",
    "Thyroid_Hormones": "DOM_THYROID",
    "Hormones": "DOM_HORMONE",
    "Aging_Metabolism": "DOM_METABOLIC",
    "Blood_Counts": "DOM_HEMATOLOGY",
    "Vitamins_Minerals": "DOM_NUTRIENT",
    "Oxidative_Stress": "DOM_INFLAMMATION",
    "Longevity_Biomarkers": "DOM_AGING",
    "General_Longevity": "DOM_AGING",
}

# ============================================================================
# DOMAIN → CONDITION MAPPING (for auto-linking)
# Maps toxin_domain values to Health_Conditions condition_id
# ============================================================================

DOMAIN_TO_CONDITION = {
    "DOM_INFLAMMATION": "COND_005",    # Chronic Inflammation
    "DOM_LIPID": "COND_010",           # Dyslipidemia
    "DOM_HEMATOLOGY": "COND_006",      # Immune Suppression (blood markers)
    "DOM_METABOLIC": "COND_003",       # Metabolic Syndrome
    "DOM_HORMONE": "COND_002",         # Hormonal Imbalance
    "DOM_THYROID": "COND_001",         # Thyroid Dysfunction
    "DOM_KIDNEY": "COND_020",          # Kidney Dysfunction (also covers liver via COND_013)
    "DOM_NUTRIENT": "COND_012",        # Detoxification Impairment (nutrient cofactors)
    "DOM_AGING": "COND_017",           # Mitochondrial Dysfunction (aging/longevity)
}

# ============================================================================
# CONFIGURATION
# ============================================================================

# API Keys from environment
ANTHROPIC_KEY = os.environ.get("ANTHROPIC_KEY")
AIRTABLE_API_KEY = os.environ.get("AIRTABLE_API_KEY")

# Airtable configuration
AIRTABLE_BASE_ID = os.environ.get("AIRTABLE_CLINICAL_BASE_ID", "app42HAczcSBeZOxD")
AIRTABLE_TABLE_NAME = "Clinical_Evidence"

# Related tables for auto-linking
HEALTH_CONDITIONS_TABLE = "Health_Conditions"
SYMPTOM_CLUSTERS_TABLE = "Symptom_Clusters"

# Domain detection keywords for longevity research
DOMAIN_KEYWORDS = {
    "Aging_Metabolism": [
        "glucose", "insulin", "hba1c", "fasting glucose", "ogtt",
        "metabolic", "diabetes", "prediabetes", "insulin resistance",
        "triglycerides", "cholesterol", "ldl", "hdl", "apob"
    ],
    "Inflammation": [
        "crp", "c-reactive protein", "hscrp", "interleukin", "il-6",
        "tnf-alpha", "inflammation", "inflammatory", "cytokine",
        "inflammaging", "chronic inflammation", "nf-kb"
    ],
    "Oxidative_Stress": [
        "oxidative stress", "reactive oxygen", "ros", "antioxidant",
        "glutathione", "malondialdehyde", "8-ohdg", "isoprostanes",
        "lipid peroxidation", "superoxide dismutase", "catalase"
    ],
    "Cardiovascular": [
        "homocysteine", "lp(a)", "lipoprotein", "apolipoprotein",
        "arterial stiffness", "pulse wave", "carotid", "atherosclerosis",
        "endothelial", "blood pressure", "hypertension", "cardiovascular"
    ],
    "Kidney_Liver": [
        "creatinine", "egfr", "bun", "cystatin c", "uric acid",
        "alt", "ast", "ggt", "albumin", "bilirubin", "liver function",
        "kidney function", "renal", "hepatic"
    ],
    "Thyroid_Hormones": [
        "tsh", "t3", "t4", "thyroid", "free t3", "free t4",
        "testosterone", "estrogen", "dhea", "cortisol", "hormone",
        "endocrine", "igf-1", "growth hormone"
    ],
    "Blood_Counts": [
        "hemoglobin", "hematocrit", "rbc", "wbc", "platelet",
        "neutrophil", "lymphocyte", "monocyte", "complete blood count",
        "anemia", "red blood cell", "white blood cell", "rdw"
    ],
    "Longevity_Biomarkers": [
        "telomere", "telomerase", "epigenetic clock", "dna methylation",
        "biological age", "senescence", "senolytics", "nad+", "sirtuin",
        "ampk", "mtor", "autophagy", "mitochondrial function"
    ],
    "Vitamins_Minerals": [
        "vitamin d", "25-hydroxyvitamin", "vitamin b12", "folate",
        "ferritin", "iron", "zinc", "magnesium", "selenium",
        "omega-3", "vitamin k", "deficiency"
    ]
}

# High-impact journals for longevity/aging research (expanded)
TOP_JOURNALS = [
    "nature aging", "cell metabolism", "aging cell",
    "nature medicine", "lancet healthy longevity", "jama",
    "nejm", "lancet", "bmj", "journals of gerontology",
    "geroscience", "aging", "nature", "science", "cell",
    "annals of internal medicine", "circulation", "diabetes care",
    "plos one", "plos medicine", "scientific reports",
    "frontiers in aging", "experimental gerontology"
]

# =============================================================================
# VALID EVIDENCE TYPES (must match Airtable single select options)
# =============================================================================

VALID_EVIDENCE_TYPES = [
    "Meta-analysis", "Systematic Review", "RCT", "Prospective Cohort",
    "Cohort", "Case-control", "Cross-sectional", "NHANES", "Biomonitoring",
    "Case Series", "In Vitro", "Animal", "Mechanistic", "Other"
]

def sanitize_evidence_type(evidence_type):
    """
    Sanitize evidence_type to match valid Airtable single select options.
    Maps common variations and invalid values to valid options.
    """
    if not evidence_type:
        return "Other"
    
    et = str(evidence_type).strip().strip('"').strip("'")
    et_lower = et.lower()
    
    # Direct match (case-insensitive)
    for valid in VALID_EVIDENCE_TYPES:
        if et_lower == valid.lower():
            return valid
    
    # Map common variations
    mapping = {
        "meta analysis": "Meta-analysis",
        "systematic review and meta-analysis": "Meta-analysis",
        "systematic review": "Systematic Review",
        "review": "Systematic Review",
        "randomized controlled trial": "RCT",
        "randomised controlled trial": "RCT",
        "randomized clinical trial": "RCT",
        "clinical trial": "RCT",
        "prospective cohort": "Prospective Cohort",
        "prospective study": "Prospective Cohort",
        "prospective": "Prospective Cohort",
        "longitudinal": "Prospective Cohort",
        "cohort study": "Cohort",
        "retrospective cohort": "Cohort",
        "retrospective": "Cohort",
        "case control": "Case-control",
        "case-control study": "Case-control",
        "cross sectional": "Cross-sectional",
        "cross-sectional study": "Cross-sectional",
        "population-based": "Cross-sectional",
        "observational": "Cross-sectional",
        "nhanes analysis": "NHANES",
        "nhanes study": "NHANES",
        "biomonitoring study": "Biomonitoring",
        "case series": "Case Series",
        "case report": "Case Series",
        "in vitro": "In Vitro",
        "in-vitro": "In Vitro",
        "animal study": "Animal",
        "animal model": "Animal",
        "in vivo": "Animal",
        "mechanistic study": "Mechanistic",
        "unknown": "Other",
        "not reported": "Other",
        "n/a": "Other",
        "": "Other",
    }
    
    if et_lower in mapping:
        return mapping[et_lower]
    
    # Partial matching
    if "meta" in et_lower and "analy" in et_lower:
        return "Meta-analysis"
    if "systematic" in et_lower:
        return "Systematic Review"
    if "random" in et_lower or "rct" in et_lower:
        return "RCT"
    if "prospective" in et_lower or "longitudinal" in et_lower:
        return "Prospective Cohort"
    if "cohort" in et_lower:
        return "Cohort"
    if "case-control" in et_lower or "case control" in et_lower:
        return "Case-control"
    if "cross" in et_lower and "section" in et_lower:
        return "Cross-sectional"
    if "nhanes" in et_lower:
        return "NHANES"
    if "biomonitor" in et_lower:
        return "Biomonitoring"
    if "vitro" in et_lower:
        return "In Vitro"
    if "animal" in et_lower or "mouse" in et_lower or "rat" in et_lower:
        return "Animal"
    
    return "Other"


# Statistics tracking
stats = {
    "total_searched": 0,
    "abstracts_fetched": 0,
    "duplicates_skipped": 0,
    "below_threshold": 0,
    "added_to_airtable": 0,
    "errors": 0
}

# ============================================================================
# AIRTABLE CACHING FUNCTIONS
# ============================================================================

_existing_titles_cache = None
_health_conditions_cache = None
_symptom_clusters_cache = None

def get_airtable_headers():
    """Return standard Airtable API headers"""
    return {
        "Authorization": f"Bearer {AIRTABLE_API_KEY}",
        "Content-Type": "application/json"
    }

def get_existing_titles():
    """Fetch all existing study titles for duplicate detection"""
    global _existing_titles_cache
    if _existing_titles_cache is not None:
        return _existing_titles_cache
    
    url = f"https://api.airtable.com/v0/{AIRTABLE_BASE_ID}/{AIRTABLE_TABLE_NAME}"
    headers = get_airtable_headers()
    
    titles = set()
    offset = None
    
    try:
        while True:
            params = {"fields[]": "study_title", "pageSize": 100}
            if offset:
                params["offset"] = offset
            
            response = requests.get(url, headers=headers, params=params, timeout=30)
            
            if response.status_code == 404:
                print(f"📋 Table {AIRTABLE_TABLE_NAME} not found - will create records fresh")
                _existing_titles_cache = set()
                return _existing_titles_cache
            
            response.raise_for_status()
            data = response.json()
            
            for record in data.get("records", []):
                title = record.get("fields", {}).get("study_title", "")
                if title:
                    titles.add(title.lower().strip())
            
            offset = data.get("offset")
            if not offset:
                break
        
        print(f"📋 Loaded {len(titles)} existing study titles for dedup")
        _existing_titles_cache = titles
        return _existing_titles_cache
    
    except requests.exceptions.RequestException as e:
        print(f"⚠️ Warning loading existing titles: {e}")
        _existing_titles_cache = set()
        return _existing_titles_cache


def load_health_conditions():
    """Load Health_Conditions table into cache for auto-linking"""
    global _health_conditions_cache
    if _health_conditions_cache is not None:
        return _health_conditions_cache
    
    print("📋 Loading Health_Conditions table...")
    url = f"https://api.airtable.com/v0/{AIRTABLE_BASE_ID}/{HEALTH_CONDITIONS_TABLE}"
    headers = get_airtable_headers()
    records = []
    offset = None
    
    try:
        while True:
            params = {"pageSize": 100}
            if offset:
                params["offset"] = offset
            
            response = requests.get(url, headers=headers, params=params, timeout=30)
            if response.status_code != 200:
                print(f"   ⚠️ Error loading Health_Conditions: {response.status_code}")
                _health_conditions_cache = {}
                return _health_conditions_cache
            
            data = response.json()
            records.extend(data.get("records", []))
            
            offset = data.get("offset")
            if not offset:
                break
        
        # Build cache: condition_id → Airtable record ID
        _health_conditions_cache = {}
        for rec in records:
            cond_id = rec["fields"].get("condition_id")
            if cond_id:
                _health_conditions_cache[cond_id] = rec["id"]
        
        print(f"   ✅ Loaded {len(_health_conditions_cache)} health conditions")
    except Exception as e:
        print(f"   ⚠️ Exception loading Health_Conditions: {e}")
        _health_conditions_cache = {}
    
    return _health_conditions_cache


def load_symptom_clusters():
    """Load Symptom_Clusters table into cache for auto-linking"""
    global _symptom_clusters_cache
    if _symptom_clusters_cache is not None:
        return _symptom_clusters_cache
    
    print("📋 Loading Symptom_Clusters table...")
    url = f"https://api.airtable.com/v0/{AIRTABLE_BASE_ID}/{SYMPTOM_CLUSTERS_TABLE}"
    headers = get_airtable_headers()
    records = []
    offset = None
    
    try:
        while True:
            params = {"pageSize": 100}
            if offset:
                params["offset"] = offset
            
            response = requests.get(url, headers=headers, params=params, timeout=30)
            if response.status_code != 200:
                print(f"   ⚠️ Error loading Symptom_Clusters: {response.status_code}")
                _symptom_clusters_cache = {}
                return _symptom_clusters_cache
            
            data = response.json()
            records.extend(data.get("records", []))
            
            offset = data.get("offset")
            if not offset:
                break
        
        # Build cache: primary_condition_id → [Airtable record IDs]
        _symptom_clusters_cache = {}
        for rec in records:
            cond_id = rec["fields"].get("primary_condition_id")
            if cond_id:
                if cond_id not in _symptom_clusters_cache:
                    _symptom_clusters_cache[cond_id] = []
                _symptom_clusters_cache[cond_id].append(rec["id"])
        
        print(f"   ✅ Loaded {len(records)} symptom clusters")
    except Exception as e:
        print(f"   ⚠️ Exception loading Symptom_Clusters: {e}")
        _symptom_clusters_cache = {}
    
    return _symptom_clusters_cache


# ============================================================================
# LINKING FUNCTIONS
# ============================================================================

def find_health_condition_link(condition_id):
    """Find the Airtable record ID for a given condition_id"""
    if not condition_id:
        return None
    cache = load_health_conditions()
    return cache.get(condition_id)


def find_symptom_cluster_links(condition_id):
    """Find all Symptom_Clusters that match the given condition_id"""
    if not condition_id:
        return []
    cache = load_symptom_clusters()
    return cache.get(condition_id, [])


def get_condition_from_domain(domain):
    """Map toxin_domain value to condition_id for linking"""
    return DOMAIN_TO_CONDITION.get(domain)


# ============================================================================
# PUBMED API FUNCTIONS
# ============================================================================

def search_pubmed(query, max_results=30, date_after="2024-01-01"):
    """Search PubMed and return list of PMIDs"""
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    
    full_query = f"({query}) AND human[MeSH Terms] AND (longevity OR aging OR biomarker OR healthspan)"
    
    params = {
        "db": "pubmed",
        "term": full_query,
        "retmax": max_results,
        "retmode": "json",
        "sort": "date",
        "mindate": date_after.replace("-", "/"),
        "maxdate": datetime.now().strftime("%Y/%m/%d"),
        "datetype": "pdat"
    }
    
    try:
        response = requests.get(base_url, params=params, timeout=30)
        response.raise_for_status()
        data = response.json()
        pmids = data.get("esearchresult", {}).get("idlist", [])
        stats["total_searched"] += len(pmids)
        return pmids
    except requests.exceptions.RequestException as e:
        print(f"⚠️ PubMed search error: {e}")
        stats["errors"] += 1
        return []


def fetch_pubmed_abstract(pmid):
    """Fetch article metadata and abstract from PubMed with robust encoding"""
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    params = {
        "db": "pubmed",
        "id": pmid,
        "retmode": "xml"
    }
    
    try:
        response = requests.get(base_url, params=params, timeout=30)
        response.raise_for_status()
        
        # Robust encoding handling
        try:
            content = response.content.decode('utf-8', errors='replace')
        except:
            content = response.text
        
        # Remove invalid XML characters
        content = ''.join(char for char in content if char.isprintable() or char in '\n\r\t')
        
        root = ET.fromstring(content.encode('utf-8'))
        article = root.find(".//PubmedArticle")
        
        if article is None:
            return None
        
        title = article.findtext(".//ArticleTitle", "")
        
        # Handle multiple AbstractText elements (structured abstracts)
        abstract_parts = []
        for abstract_elem in article.findall(".//AbstractText"):
            if abstract_elem.text:
                label = abstract_elem.get("Label", "")
                if label:
                    abstract_parts.append(f"{label}: {abstract_elem.text}")
                else:
                    abstract_parts.append(abstract_elem.text)
        abstract = " ".join(abstract_parts) if abstract_parts else ""
        
        journal = article.findtext(".//Journal/Title", "")
        year = article.findtext(".//PubDate/Year", "")
        
        if not year:
            medline_date = article.findtext(".//PubDate/MedlineDate", "")
            if medline_date:
                year = medline_date[:4]
        
        if not year:
            year = str(datetime.now().year)
        
        authors = []
        for author in article.findall(".//Author"):
            lastname = author.findtext("LastName", "")
            initials = author.findtext("Initials", "")
            if lastname:
                authors.append(f"{lastname} {initials}")
        
        stats["abstracts_fetched"] += 1
        
        return {
            "pmid": pmid,
            "title": title,
            "abstract": abstract,
            "journal": journal,
            "year": year,
            "authors": ", ".join(authors[:5]) + (" et al." if len(authors) > 5 else ""),
            "url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
        }
    
    except ET.ParseError as e:
        print(f"   ⚠️ XML parse error for {pmid}: {e}")
        stats["errors"] += 1
        return None
    except Exception as e:
        print(f"⚠️ Error fetching PMID {pmid}: {e}")
        stats["errors"] += 1
        return None


# ============================================================================
# DOMAIN DETECTION
# ============================================================================

def detect_domain(query, abstract):
    """Auto-detect longevity domain based on keywords"""
    text = (query + " " + abstract).lower()
    
    domain_scores = {}
    for domain, keywords in DOMAIN_KEYWORDS.items():
        score = sum(1 for kw in keywords if kw.lower() in text)
        if score > 0:
            domain_scores[domain] = score
    
    if domain_scores:
        return max(domain_scores, key=domain_scores.get)
    return "General_Longevity"


# ============================================================================
# EVIDENCE SCORING (v1.2 - Calibrated for longevity/wellness literature)
# ============================================================================

def calculate_stars(evidence_type, sample_size, journal, effect_size_reported):
    """
    Calculate evidence strength score (1-5 stars)
    
    v1.2 - Calibrated for longevity research:
    - Cross-sectional/population studies score higher (common in aging research)
    - Large cohort studies (NHANES, UK Biobank) get sample size bonus
    - Expanded journal list for aging-specific high-impact journals
    """
    score = 0
    
    evidence_type_lower = evidence_type.lower() if evidence_type else ""
    
    # Evidence type scoring - longevity calibrated
    if "meta-analysis" in evidence_type_lower or "systematic review" in evidence_type_lower:
        score += 2.5
    elif "randomized" in evidence_type_lower or "rct" in evidence_type_lower:
        score += 2.5
    elif "prospective cohort" in evidence_type_lower or "longitudinal" in evidence_type_lower:
        score += 2.0
    elif "cohort" in evidence_type_lower or "case-control" in evidence_type_lower:
        score += 1.5
    elif "cross-sectional" in evidence_type_lower or "nhanes" in evidence_type_lower or "population" in evidence_type_lower:
        score += 1.0  # Key change: cross-sectional now scores 1.0 (was 0.5)
    elif "case series" in evidence_type_lower or "case report" in evidence_type_lower:
        score += 0.5
    else:
        score += 0.5
    
    # Sample size scoring - calibrated for large population studies
    if sample_size:
        try:
            size_str = str(sample_size).replace(",", "").replace(" ", "")
            n = int(''.join(filter(str.isdigit, size_str.split()[0] if size_str.split() else size_str)))
            
            if n >= 10000:
                score += 1.5  # Large population (UK Biobank, NHANES scale)
            elif n >= 1000:
                score += 1.0
            elif n >= 500:
                score += 0.75
            elif n >= 100:
                score += 0.5
            elif n >= 50:
                score += 0.25
        except (ValueError, IndexError):
            pass
    
    # Journal quality scoring
    journal_lower = journal.lower() if journal else ""
    if any(top in journal_lower for top in TOP_JOURNALS):
        score += 1.0
    else:
        score += 0.5
    
    # Effect size reporting bonus
    if effect_size_reported and effect_size_reported.lower() not in ["not reported", "n/a", "none", ""]:
        score += 0.5
        # Bonus for dose-response or quantified relationships
        if any(x in effect_size_reported.lower() for x in ["per", "dose", "quartile", "tertile", "trend"]):
            score += 0.25
    
    return min(max(round(score), 1), 5)
    

def format_stars(num):
    """Format star rating as emoji string."""
    return f"{num} " + "⭐️" * num


# ============================================================================
# CLAUDE AI EXTRACTION
# ============================================================================

def ask_claude(article, domain):
    """Use Claude to extract structured metadata from abstract"""
    client = Anthropic(api_key=ANTHROPIC_KEY)
    
    prompt = f"""Analyze this longevity/healthspan research article and extract structured information.

TITLE: {article['title']}

ABSTRACT: {article['abstract']}

JOURNAL: {article['journal']}

DETECTED DOMAIN: {domain}

Extract the following in JSON format:
{{
    "evidence_type": "Study design (e.g., 'Meta-analysis', 'RCT', 'Prospective cohort', 'Cross-sectional', 'NHANES')",
    "sample_size": "Number of participants or 'Not reported'",
    "biomarkers_studied": ["list", "of", "biomarkers"],
    "key_findings": "2-3 sentence summary of main findings relevant to longevity/healthspan",
    "effect_size": "Quantified effect (HR, OR, correlation, β, etc.) or 'Not reported'",
    "clinical_relevance": "Practical implications for healthspan optimization",
    "limitations": "Key study limitations"
}}

Return ONLY valid JSON, no markdown formatting."""

    try:
        response = client.messages.create(
            model="claude-sonnet-4-20250514",
            max_tokens=1500,
            messages=[{"role": "user", "content": prompt}]
        )
        
        text = response.content[0].text.strip()
        
        if text.startswith("```"):
            text = text.split("\n", 1)[1]
        if text.endswith("```"):
            text = text.rsplit("\n", 1)[0]
        text = text.replace("```json", "").replace("```", "").strip()
        
        return json.loads(text)
    
    except json.JSONDecodeError as e:
        print(f"⚠️ JSON parse error: {e}")
        stats["errors"] += 1
        return None
    except Exception as e:
        print(f"⚠️ Claude API error: {e}")
        stats["errors"] += 1
        return None


# ============================================================================
# AIRTABLE UPLOAD WITH AUTO-LINKING
# ============================================================================

def get_next_evidence_id():
    """Generate unique evidence ID"""
    now = datetime.now()
    return f"LONG_{now.strftime('%m%d%H%M%S')}{stats['added_to_airtable']:03d}"


def add_to_airtable(article, extracted, stars, domain):
    """Upload study to Airtable Clinical_Evidence table with auto-linking
    
    v1.2: Now auto-links to Health_Conditions and Symptom_Clusters based on domain
    """
    url = f"https://api.airtable.com/v0/{AIRTABLE_BASE_ID}/{AIRTABLE_TABLE_NAME}"
    headers = get_airtable_headers()
    
    # Map domain to exact Airtable single select option
    mapped_domain = DOMAIN_MAP.get(domain, "DOM_AGING")
    
    # Get condition_id from domain for linking
    condition_id = get_condition_from_domain(mapped_domain)
    
    # Format biomarkers as comma-separated string
    biomarkers = extracted.get("biomarkers_studied", [])
    if isinstance(biomarkers, list):
        biomarkers_str = ", ".join(str(b) for b in biomarkers)
    else:
        biomarkers_str = str(biomarkers) if biomarkers else ""
    
    # Build record with exact field names from Clinical_Evidence table
    record = {
        "fields": {
            "evidence_id": get_next_evidence_id(),
            "study_title": str(article.get("title", ""))[:500],
            "authors_year": f"{article.get('year', '2024')}-01-01",
            "journal": str(article.get("journal", ""))[:200],
            "toxin_domain": mapped_domain,
            "condition_id": condition_id or "",  # NEW: Store condition_id
            "evidence_type": sanitize_evidence_type(extracted.get("evidence_type", "")),
            "sample_size": str(extracted.get("sample_size", ""))[:50],
            "markers_covered": biomarkers_str[:1000],
            "key_findings": str(extracted.get("key_findings", ""))[:2000],
            "effect_size": str(extracted.get("effect_size", ""))[:200],
            "clinical_relevance": str(extracted.get("clinical_relevance", ""))[:1000],
            "limitations": str(extracted.get("limitations", ""))[:1000],
            "evidence_strength_score": format_stars(stars),
            "source_url": str(article.get("url", ""))
        }
    }
    
    # =========================================================================
    # AUTO-LINKING (v1.2)
    # =========================================================================
    try:
        # 1. Link to Health_Conditions
        if condition_id:
            health_condition_rec_id = find_health_condition_link(condition_id)
            if health_condition_rec_id:
                record["fields"]["health_condition_link"] = [health_condition_rec_id]
                print(f"   🔗 Linked to Health_Conditions: {condition_id}")
        
        # 2. Link to Symptom_Clusters
        if condition_id:
            symptom_cluster_rec_ids = find_symptom_cluster_links(condition_id)
            if symptom_cluster_rec_ids:
                record["fields"]["symptom_clusters_link"] = symptom_cluster_rec_ids
                print(f"   🔗 Linked to {len(symptom_cluster_rec_ids)} Symptom_Clusters")
    
    except Exception as e:
        print(f"   ⚠️ Linking error (continuing anyway): {e}")
    
    # Remove empty string values to avoid Airtable errors
    record["fields"] = {k: v for k, v in record["fields"].items() if v and (not isinstance(v, str) or v.strip())}
    
    try:
        response = requests.post(url, headers=headers, json=record, timeout=30)
        response.raise_for_status()
        stats["added_to_airtable"] += 1
        print(f"   ✅ Added: {format_stars(stars)} | {mapped_domain}")
        return True
    except requests.exceptions.RequestException as e:
        print(f"   ❌ Airtable error: {e}")
        if hasattr(e, 'response') and e.response is not None:
            print(f"   Response: {e.response.text[:500]}")
        stats["errors"] += 1
        return False


# ============================================================================
# MAIN EXECUTION
# ============================================================================

def load_config():
    """Load configuration from config.yaml"""
    try:
        with open("config.yaml", "r") as f:
            return yaml.safe_load(f)
    except FileNotFoundError:
        print("⚠️ config.yaml not found, using defaults")
        return {
            "search": {
                "keywords": [
                    "longevity biomarkers blood",
                    "biological age blood markers",
                    "healthspan inflammation CRP",
                    "aging metabolism glucose insulin",
                    "cardiovascular biomarkers mortality"
                ],
                "date_after": "2024-06-01",
                "max_results": 25
            },
            "scoring": {
                "min_score_to_save": 3
            }
        }


def main():
    """Main execution flow"""
    print("=" * 60)
    print("🧬 LONGEVITY EVIDENCE SCOUT v1.2 (with Auto-Linking)")
    print("=" * 60)
    print(f"⏰ Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print()
    
    # Validate environment
    if not ANTHROPIC_KEY:
        print("❌ ANTHROPIC_KEY not set!")
        sys.exit(1)
    if not AIRTABLE_API_KEY:
        print("❌ AIRTABLE_API_KEY not set!")
        sys.exit(1)
    
    print(f"📊 Airtable Base: {AIRTABLE_BASE_ID}")
    print(f"📋 Table: {AIRTABLE_TABLE_NAME}")
    print()
    
    # Load config
    config = load_config()
    keywords = config.get("search", {}).get("keywords", [])
    date_after = config.get("search", {}).get("date_after", "2024-06-01")
    max_results = config.get("search", {}).get("max_results", 25)
    min_score = config.get("scoring", {}).get("min_score_to_save", 3)
    
    print(f"📅 Date filter: After {date_after}")
    print(f"🔍 Keyword groups: {len(keywords)}")
    print(f"⭐ Min score threshold: {min_score}")
    print()
    
    # =========================================================================
    # PRE-LOAD CACHES FOR AUTO-LINKING (v1.2)
    # =========================================================================
    print("📥 Loading related tables for auto-linking...")
    load_health_conditions()
    load_symptom_clusters()
    print()
    
    # Load existing titles for dedup
    existing_titles = get_existing_titles()
    
    # Process each keyword group
    for i, keyword in enumerate(keywords, 1):
        print(f"\n{'='*60}")
        print(f"🔎 [{i}/{len(keywords)}] Searching: {keyword[:50]}...")
        print("=" * 60)
        
        pmids = search_pubmed(keyword, max_results, date_after)
        print(f"   Found {len(pmids)} articles")
        
        for pmid in pmids:
            # Rate limiting for PubMed
            time.sleep(0.4)
            
            # Fetch abstract
            article = fetch_pubmed_abstract(pmid)
            if not article or not article.get("abstract"):
                continue
            
            # Check for duplicates
            title_lower = article["title"].lower().strip()
            if title_lower in existing_titles:
                print(f"   ⏭️ Duplicate: {article['title'][:50]}...")
                stats["duplicates_skipped"] += 1
                continue
            
            # Detect domain
            domain = detect_domain(keyword, article["abstract"])
            
            # Extract with Claude
            print(f"   🤖 Analyzing: {article['title'][:50]}...")
            extracted = ask_claude(article, domain)
            if not extracted:
                continue
            
            # Score evidence
            stars = calculate_stars(
                extracted.get("evidence_type", ""),
                extracted.get("sample_size", ""),
                article["journal"],
                extracted.get("effect_size", "")
            )
            
            # Check threshold
            if stars < min_score:
                print(f"   ⏭️ Below threshold ({format_stars(stars)})")
                stats["below_threshold"] += 1
                continue
            
            # Add to Airtable (now with auto-linking)
            if add_to_airtable(article, extracted, stars, domain):
                existing_titles.add(title_lower)
    
    # Print summary
    print("\n" + "=" * 60)
    print("📊 SUMMARY")
    print("=" * 60)
    print(f"   🔍 Articles searched: {stats['total_searched']}")
    print(f"   📄 Abstracts fetched: {stats['abstracts_fetched']}")
    print(f"   ⏭️ Duplicates skipped: {stats['duplicates_skipped']}")
    print(f"   📉 Below threshold: {stats['below_threshold']}")
    print(f"   ✅ Added to Airtable: {stats['added_to_airtable']}")
    print(f"   ❌ Errors: {stats['errors']}")
    print(f"   ⏰ Completed: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 60)
    
    sys.exit(0)


if __name__ == "__main__":
    main()
