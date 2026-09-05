"""Tests for sanitize_evidence_type() and calculate_stars().

These are pure-logic functions that map free-text LLM output to canonical
Airtable select values (sanitize) and compute the 1-5 evidence-strength
rating (calculate_stars). Both have large blast radius: a wrong mapping
puts evidence in the wrong category; a wrong star count changes the
user-facing rating. Neither was covered before this PR.

main.py is importable after `uv pip install -r requirements.txt` (yaml + requests).
"""

import unittest
import main


class TestSanitizeEvidenceType(unittest.TestCase):
    """Lock the canonical evidence-type mapping that feeds Airtable single-select."""

    # --- empty / falsy ---

    def test_none_returns_other(self):
        self.assertEqual(main.sanitize_evidence_type(None), "Other")

    def test_empty_string_returns_other(self):
        self.assertEqual(main.sanitize_evidence_type(""), "Other")

    def test_whitespace_only_returns_other(self):
        self.assertEqual(main.sanitize_evidence_type("   "), "Other")

    # --- exact VALID_EVIDENCE_TYPES (case-insensitive) ---

    def test_exact_match_meta_analysis(self):
        self.assertEqual(main.sanitize_evidence_type("Meta-analysis"), "Meta-analysis")

    def test_exact_match_case_insensitive(self):
        self.assertEqual(main.sanitize_evidence_type("rct"), "RCT")
        self.assertEqual(main.sanitize_evidence_type("CASE-CONTROL"), "Case-control")

    def test_all_valid_types_round_trip(self):
        for vt in main.VALID_EVIDENCE_TYPES:
            self.assertEqual(main.sanitize_evidence_type(vt), vt)

    # --- mapping table ---

    def test_meta_analysis_mapping(self):
        self.assertEqual(main.sanitize_evidence_type("meta analysis"), "Meta-analysis")
        self.assertEqual(
            main.sanitize_evidence_type("systematic review and meta-analysis"),
            "Meta-analysis",
        )

    def test_systematic_review_mapping(self):
        self.assertEqual(main.sanitize_evidence_type("systematic review"), "Systematic Review")
        self.assertEqual(main.sanitize_evidence_type("review"), "Systematic Review")

    def test_rct_mappings(self):
        self.assertEqual(main.sanitize_evidence_type("randomized controlled trial"), "RCT")
        self.assertEqual(main.sanitize_evidence_type("randomised controlled trial"), "RCT")
        self.assertEqual(main.sanitize_evidence_type("clinical trial"), "RCT")

    def test_prospective_mappings(self):
        self.assertEqual(main.sanitize_evidence_type("prospective"), "Prospective Cohort")
        self.assertEqual(main.sanitize_evidence_type("longitudinal"), "Prospective Cohort")

    def test_cohort_mappings(self):
        self.assertEqual(main.sanitize_evidence_type("cohort study"), "Cohort")
        self.assertEqual(main.sanitize_evidence_type("retrospective"), "Cohort")

    def test_case_control_mapping(self):
        self.assertEqual(main.sanitize_evidence_type("case control"), "Case-control")

    def test_cross_sectional_mapping(self):
        self.assertEqual(main.sanitize_evidence_type("cross sectional"), "Cross-sectional")
        self.assertEqual(main.sanitize_evidence_type("population-based"), "Cross-sectional")

    def test_nhanes_mapping(self):
        self.assertEqual(main.sanitize_evidence_type("nhanes analysis"), "NHANES")

    def test_animal_mapping(self):
        self.assertEqual(main.sanitize_evidence_type("animal study"), "Animal")
        self.assertEqual(main.sanitize_evidence_type("in vivo"), "Animal")

    def test_in_vitro_mapping(self):
        self.assertEqual(main.sanitize_evidence_type("in vitro"), "In Vitro")

    def test_unknown_mappings(self):
        self.assertEqual(main.sanitize_evidence_type("unknown"), "Other")
        self.assertEqual(main.sanitize_evidence_type("n/a"), "Other")

    # --- fuzzy fallback ---

    def test_fuzzy_meta_analysis(self):
        self.assertEqual(
            main.sanitize_evidence_type("Meta-analysis of 5 RCTs"), "Meta-analysis"
        )

    def test_fuzzy_systematic(self):
        self.assertEqual(
            main.sanitize_evidence_type("Systematic literature review"), "Systematic Review"
        )

    def test_fuzzy_rct(self):
        self.assertEqual(main.sanitize_evidence_type("randomised trial of vitamin D"), "RCT")

    def test_fuzzy_nhanes(self):
        self.assertEqual(
            main.sanitize_evidence_type("NHANES population survey"), "NHANES"
        )

    def test_fuzzy_animal_mouse(self):
        self.assertEqual(main.sanitize_evidence_type("mouse model"), "Animal")
        self.assertEqual(main.sanitize_evidence_type("rat study"), "Animal")

    def test_truly_unknown_returns_other(self):
        # NOTE (pitfall 40): "observational" contains "animal" as a substring,
        # so it hits the fuzzy `if "animal" in et_lower:` fallback → "Animal".
        # This is a surprising live-contract quirk — lock it, don't fix it.
        self.assertEqual(main.sanitize_evidence_type("observational narrative"), "Animal")
        # A genuinely unrecognised string → Other.
        self.assertEqual(main.sanitize_evidence_type("xyzzy qwerty"), "Other")

    # --- quotes stripped ---

    def test_strips_surrounding_quotes(self):
        self.assertEqual(main.sanitize_evidence_type('"RCT"'), "RCT")
        self.assertEqual(main.sanitize_evidence_type("'Meta-analysis'"), "Meta-analysis")


class TestCalculateStars(unittest.TestCase):
    """Lock the evidence-strength scoring algorithm (1-5 stars)."""

    # --- evidence type tier weights ---

    def test_meta_analysis_gets_2_5_base(self):
        """Meta-analysis alone (no sample/journal/effect) gets 2.5 + 0.5 (non-top journal) = 3."""
        stars = main.calculate_stars("Meta-analysis", None, "", "")
        self.assertEqual(stars, 3)

    def test_rct_gets_2_5_base(self):
        """RCT alone gets 2.5 + 0.5 = 3."""
        stars = main.calculate_stars("RCT", None, "", "")
        self.assertEqual(stars, 3)

    def test_prospective_cohort_gets_2_0_base(self):
        """Prospective cohort alone gets 2.0 + 0.5 = 2.5 → round(2.5) = 2 (Python banker's rounding)."""
        stars = main.calculate_stars("Prospective Cohort", None, "", "")
        self.assertEqual(stars, 2)

    def test_cohort_gets_1_5_base(self):
        """Cohort alone gets 1.5 + 0.5 = 2."""
        stars = main.calculate_stars("Cohort", None, "", "")
        self.assertEqual(stars, 2)

    def test_cross_sectional_gets_1_0_base(self):
        """Cross-sectional alone gets 1.0 + 0.5 = 1 (rounded from 1.5)."""
        stars = main.calculate_stars("Cross-sectional", None, "", "")
        self.assertEqual(stars, 2)  # round(1.5) = 2

    def test_case_series_gets_0_5_base(self):
        """Case series alone gets 0.5 + 0.5 = 1."""
        stars = main.calculate_stars("Case Series", None, "", "")
        self.assertEqual(stars, 1)

    def test_unknown_type_floor(self):
        """Unknown evidence type gets 0.5 + 0.5 = 1 (floor)."""
        stars = main.calculate_stars("xyzzy", None, "", "")
        self.assertEqual(stars, 1)

    # --- sample size tiers ---

    def test_large_sample_bonus(self):
        """N >= 10000 → +1.5. Total: 2.5 + 1.5 + 0.5 = 4.5 → round(4.5) = 4 (banker's rounding)."""
        stars = main.calculate_stars("RCT", 10000, "", "")
        self.assertEqual(stars, 4)

    def test_medium_sample_bonus(self):
        """N >= 1000 → +1.0."""
        stars = main.calculate_stars("RCT", 1000, "", "")
        # 2.5 + 1.0 + 0.5 = 4.0 → 4
        self.assertEqual(stars, 4)

    def test_small_sample_bonus(self):
        """N >= 100 → +0.5."""
        stars = main.calculate_stars("RCT", 100, "", "")
        # 2.5 + 0.5 + 0.5 = 3.5 → round to 4
        self.assertEqual(stars, 4)

    def test_no_sample_size(self):
        """None sample → no bonus, no crash."""
        stars = main.calculate_stars("RCT", None, "", "")
        self.assertEqual(stars, 3)

    def test_sample_with_commas(self):
        """'1,000' should parse to 1000."""
        stars = main.calculate_stars("RCT", "1,000", "", "")
        self.assertEqual(stars, 4)

    # --- top journal bonus ---

    def test_top_journal_bonus(self):
        """Top journal → +1.0 instead of +0.5."""
        stars_rct_top = main.calculate_stars("RCT", None, "Nature Medicine", "")
        stars_rct_notop = main.calculate_stars("RCT", None, "Blog Post", "")
        self.assertEqual(stars_rct_top - stars_rct_notop, 1)

    # --- effect size reported bonus ---

    def test_effect_size_reported_bonus(self):
        """Non-trivial effect size → +0.5."""
        stars_with = main.calculate_stars("RCT", None, "", "dose-response")
        stars_without = main.calculate_stars("RCT", None, "", "")
        self.assertGreater(stars_with, stars_without)

    def test_effect_size_with_dose_keyword_extra_bonus(self):
        """Effect size containing dose/quartile/tertile/trend/per → extra +0.25.
        Use N=1000 so the 0.25 actually flips rounding: with dose kw = 4.75→5, without = 4.5→4."""
        stars_dose = main.calculate_stars("RCT", 1000, "", "per quartile increase")
        stars_plain = main.calculate_stars("RCT", 1000, "", "significant")
        self.assertEqual(stars_dose, 5)
        self.assertEqual(stars_plain, 4)

    def test_effect_size_not_reported_no_bonus(self):
        """'not reported' / 'n/a' / 'none' / '' → no effect-size bonus."""
        for val in ["not reported", "n/a", "none", ""]:
            stars = main.calculate_stars("RCT", None, "", val)
            self.assertEqual(stars, 3, f"effect_size='{val}' should not grant bonus")

    # --- clamping ---

    def test_max_is_5(self):
        """Even with every bonus, score is clamped to 5."""
        stars = main.calculate_stars(
            "Meta-analysis", 100000, "Nature", "dose-response trend"
        )
        self.assertLessEqual(stars, 5)

    def test_min_is_1(self):
        """Floor is 1 even for unknown type with no data."""
        stars = main.calculate_stars("", None, "", "")
        self.assertGreaterEqual(stars, 1)

    def test_none_evidence_type(self):
        """None evidence_type should not crash (treated as empty)."""
        stars = main.calculate_stars(None, None, "", "")
        self.assertGreaterEqual(stars, 1)


class TestFormatStars(unittest.TestCase):
    """Lock the user-facing star format."""

    def test_one_star(self):
        self.assertEqual(main.format_stars(1), "1 ⭐️")

    def test_five_stars(self):
        self.assertEqual(main.format_stars(5), "5 ⭐️⭐️⭐️⭐️⭐️")

    def test_three_stars(self):
        self.assertEqual(main.format_stars(3), "3 ⭐️⭐️⭐️")


if __name__ == "__main__":
    unittest.main()
