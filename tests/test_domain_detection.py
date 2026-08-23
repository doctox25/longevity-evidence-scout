"""
Tests for detect_domain() and _detect_domain_fallback().

v1.4 changelog:
- FIXED: Sex hormones no longer classified as DOM_THYROID
- FIXED: RDW studies correctly go to DOM_HEMATOLOGY
- NEW: Direct DOM_ convention (no intermediate mapping)
- NEW: Title weighting (3x default) for domain detection
- NEW: Negative keyword exclusions to prevent misclassification

The domain classification determines which Airtable table and which LLM prompt
template the study is routed to. A misclassification sends the study to the
wrong domain expert, extracts wrong biomarkers, and pollutes the evidence base.

detect_domain reads domain_keywords from config.yaml (loaded once, cached).
These tests load the actual config to verify the production classification
contract, not a synthetic fixture.
"""
import unittest
import os
import yaml
import main

# Load the production config.yaml that detect_domain uses
_CONFIG_PATH = os.path.join(os.path.dirname(os.path.dirname(__file__)), "config.yaml")
with open(_CONFIG_PATH) as f:
    CONFIG = yaml.safe_load(f)

# Reset the module-level cache so our config is loaded fresh
main._domain_keywords_cache = None


class TestDetectDomain(unittest.TestCase):

    def test_testosterone_classified_as_hormone_not_thyroid(self):
        """v1.4 fix: Sex hormones no longer classified as DOM_THYROID."""
        title = "Effects of testosterone replacement on muscle mass in aging men"
        abstract = "We measured testosterone levels and SHBG in 200 community-dwelling men."
        domain = main.detect_domain(title, abstract, "", CONFIG)
        self.assertEqual(domain, "DOM_HORMONE",
            "Testosterone study must classify as DOM_HORMONE, not DOM_THYROID")

    def test_estrogen_classified_as_hormone_not_thyroid(self):
        """v1.4 fix: Sex hormones no longer classified as DOM_THYROID."""
        title = "Estrogen and cognitive decline in postmenopausal women"
        abstract = "Serum estradiol and estrogen levels were measured..."
        domain = main.detect_domain(title, abstract, "", CONFIG)
        self.assertEqual(domain, "DOM_HORMONE")

    def test_rdw_classified_as_hematology(self):
        """v1.4 fix: RDW studies correctly go to DOM_HEMATOLOGY."""
        title = "Red cell distribution width predicts mortality in older adults"
        abstract = "RDW and anisocytosis were measured in 1000 participants."
        domain = main.detect_domain(title, abstract, "", CONFIG)
        self.assertEqual(domain, "DOM_HEMATOLOGY")

    def test_tsh_classified_as_thyroid_not_hormone(self):
        """Negative keyword exclusion: thyroid terms exclude DOM_HORMONE."""
        title = "TSH and thyroid function in older adults"
        abstract = "We measured TSH, free T3, and free T4 levels..."
        domain = main.detect_domain(title, abstract, "", CONFIG)
        self.assertEqual(domain, "DOM_THYROID")

    def test_thyroid_keyword_classified_as_thyroid(self):
        title = "Subclinical hypothyroidism and cardiovascular risk"
        abstract = "Thyroid function was assessed via TSH and thyroxine..."
        domain = main.detect_domain(title, abstract, "", CONFIG)
        self.assertEqual(domain, "DOM_THYROID")

    def test_metabolic_markers_classified_correctly(self):
        title = "HbA1c and insulin resistance in prediabetic adults"
        abstract = "Glucose and HbA1c were measured in 500 adults."
        domain = main.detect_domain(title, abstract, "", CONFIG)
        self.assertEqual(domain, "DOM_METABOLIC")

    def test_lipid_markers_classified_correctly(self):
        title = "LDL cholesterol and cardiovascular events"
        abstract = "Lipid profiles including LDL and cholesterol were analyzed."
        domain = main.detect_domain(title, abstract, "", CONFIG)
        self.assertEqual(domain, "DOM_LIPID")

    def test_default_fallback_to_aging_when_no_match(self):
        # Use words that do NOT appear in any domain_keywords config entry.
        # "aging" alone matches DOM_LIVER's keywords; use unrelated terms.
        title = "Survey of weather patterns and seasonal variation"
        abstract = "We reviewed non-biomarker indicators in a general context."
        domain = main.detect_domain(title, abstract, "", CONFIG)
        self.assertEqual(domain, "DOM_AGING",
            "No matching domain keywords must fall back to DOM_AGING")

    def test_title_weighting_3x(self):
        """Title keywords count 3x (title_weight=3) vs abstract (1x)."""
        # Title has one keyword (score 3), abstract has two from a different domain (score 2).
        # Title domain should win because 3 > 2.
        title = "Testosterone supplementation study"  # DOM_HORMONE: testosterone → 3 pts
        abstract = "We measured glucose insulin HbA1c in participants."  # DOM_METABOLIC: 3 keywords → 3 pts
        domain = main.detect_domain(title, abstract, "", CONFIG)
        # Both score 3, but the order of iteration matters. Let's use a clearer example:
        # Title with 1 thyroid keyword (3 pts) vs abstract with 2 metabolic keywords (2 pts)
        title = "TSH reference ranges"  # DOM_THYROID: TSH in title → 3 pts
        abstract = "We measured glucose and insulin in participants."  # DOM_METABOLIC: 2 keywords → 2 pts
        domain = main.detect_domain(title, abstract, "", CONFIG)
        self.assertEqual(domain, "DOM_THYROID",
            "Title keyword (3pts) must outweigh 2 abstract keywords (2pts)")

    def test_negative_keyword_excludes_domain(self):
        """If a negative keyword for DOM_HORMONE appears, DOM_HORMONE is skipped."""
        # A study that mentions both testosterone AND thyroid should NOT
        # classify as DOM_HORMONE (negative keyword 'thyroid' excludes it).
        # It will instead fall through to DOM_THYROID (positive match).
        title = "Testosterone and thyroid interaction in aging men"
        abstract = "We measured testosterone, TSH, and free T4 levels."
        domain = main.detect_domain(title, abstract, "", CONFIG)
        self.assertNotEqual(domain, "DOM_HORMONE",
            "DOM_HORMONE must be excluded when thyroid negative keyword matches")


class TestDetectDomainFallback(unittest.TestCase):
    """The fallback runs when config has no domain_keywords."""

    def setUp(self):
        main._domain_keywords_cache = None

    def tearDown(self):
        main._domain_keywords_cache = None

    def test_fallback_testosterone_to_hormone(self):
        title = "Testosterone and aging"
        abstract = ""
        domain = main._detect_domain_fallback(title, abstract, "")
        self.assertEqual(domain, "DOM_HORMONE")

    def test_fallback_tsh_to_thyroid(self):
        title = "TSH reference ranges"
        abstract = ""
        domain = main._detect_domain_fallback(title, abstract, "")
        self.assertEqual(domain, "DOM_THYROID")

    def test_fallback_rdw_to_hematology(self):
        title = "RDW as a biomarker"
        abstract = ""
        domain = main._detect_domain_fallback(title, abstract, "")
        self.assertEqual(domain, "DOM_HEMATOLOGY")

    def test_fallback_default_aging(self):
        title = "General wellness study"
        abstract = "No specific biomarkers mentioned here."
        domain = main._detect_domain_fallback(title, abstract, "")
        self.assertEqual(domain, "DOM_AGING")

    def test_fallback_crp_to_inflammation(self):
        title = "CRP and longevity"
        abstract = ""
        domain = main._detect_domain_fallback(title, abstract, "")
        self.assertEqual(domain, "DOM_INFLAMMATION")

    def test_fallback_vitamin_d_to_nutrient(self):
        title = "Vitamin D deficiency in older adults"
        abstract = ""
        domain = main._detect_domain_fallback(title, abstract, "")
        self.assertEqual(domain, "DOM_NUTRIENT")


if __name__ == "__main__":
    unittest.main()
