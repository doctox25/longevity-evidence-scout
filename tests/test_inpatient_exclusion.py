"""
Tests for is_inpatient_study() — v1.5 inpatient/acute-care exclusion filter.

This function runs BEFORE the LLM API call to save API costs. If it fails to
match an inpatient study, that study gets sent to the LLM and potentially
uploaded to Airtable's Clinical_Evidence table with the wrong population type
(inpatient instead of community/outpatient). If it false-positives a community
study, that study is silently excluded from the evidence base.

main.py is importable (only stdlib + requests + pyyaml at module scope).
"""
import unittest
import main


class TestInpatientExclusion(unittest.TestCase):

    def test_icu_study_is_excluded(self):
        title = "Mortality predictors in ICU patients with sepsis"
        abstract = "We retrospectively reviewed intensive care unit admissions..."
        excluded, phrase = main.is_inpatient_study(title, abstract)
        self.assertTrue(excluded, "ICU study must be excluded")
        self.assertIn(phrase, ["icu patients", "intensive care unit"])

    def test_hemodialysis_study_is_excluded(self):
        title = "Inflammatory markers in hemodialysis patients"
        abstract = "Long-term outcomes for end-stage renal disease..."
        excluded, phrase = main.is_inpatient_study(title, abstract)
        self.assertTrue(excluded)
        self.assertIn(phrase, ["hemodialysis patients", "end-stage renal disease", "end stage renal"])

    def test_postoperative_study_is_excluded(self):
        title = "Postoperative outcomes in cardiac surgery patients"
        abstract = "Post-operative recovery was assessed..."
        excluded, phrase = main.is_inpatient_study(title, abstract)
        self.assertTrue(excluded)
        # Should match one of the surgical/post-op phrases
        self.assertIn(phrase, ["post-operative", "postoperative", "cardiac surgery patients", "surgical patients"])

    def test_community_study_is_not_excluded(self):
        title = "Longitudinal biomarkers of aging in community-dwelling adults"
        abstract = "We followed 500 healthy participants over 10 years measuring..."
        excluded, phrase = main.is_inpatient_study(title, abstract)
        self.assertFalse(excluded, "Community/outpatient study must NOT be excluded")
        self.assertEqual(phrase, "")

    def test_rct_with_outpatient_population_is_not_excluded(self):
        title = "Effect of metformin on biological aging in overweight adults"
        abstract = "A randomized controlled trial of 300 outpatient adults aged 40-65..."
        excluded, phrase = main.is_inpatient_study(title, abstract)
        self.assertFalse(excluded)

    def test_empty_title_and_abstract(self):
        excluded, phrase = main.is_inpatient_study("", "")
        self.assertFalse(excluded)
        self.assertEqual(phrase, "")

    def test_case_insensitive_matching(self):
        """Phrases are stored lowercase; combined text is lowercased."""
        title = "Outcomes in HOSPITALIZED PATIENTS with acute illness"
        abstract = ""
        excluded, phrase = main.is_inpatient_study(title, abstract)
        self.assertTrue(excluded)
        self.assertEqual(phrase, "hospitalized patients")

    def test_matched_phrase_returned_for_audit_trail(self):
        """The function returns the matched phrase for logging/debugging."""
        title = "Sepsis patients in the emergency department"
        abstract = ""
        excluded, phrase = main.is_inpatient_study(title, abstract)
        self.assertTrue(excluded)
        self.assertIn(phrase, ["sepsis patients", "emergency department"])

    def test_acute_kidney_injury_excluded(self):
        title = "Biomarker trends in acute kidney injury"
        abstract = "Hospitalized patients with acute kidney injury..."
        excluded, phrase = main.is_inpatient_study(title, abstract)
        self.assertTrue(excluded)

    def test_critically_ill_excluded(self):
        title = "Inflammation in critically ill patients"
        abstract = ""
        excluded, phrase = main.is_inpatient_study(title, abstract)
        self.assertTrue(excluded)
        self.assertEqual(phrase, "critically ill")


if __name__ == "__main__":
    unittest.main()
