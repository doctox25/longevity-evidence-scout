"""
Tests for ask_llm() — LLM extraction with fence-stripping and JSON parsing.

PR #3 (feat: switch extraction to xAI Grok 4.6):
- Renamed ask_claude → ask_llm
- Changed API endpoint from Anthropic to xAI (https://api.x.ai/v1/chat/completions)
- Bearer auth with XAI_API_KEY
- Model from XAI_MODEL env var with empty-string-safe default
- Fence-stripping and JSON parse logic unchanged from the Claude era

The fence-stripping logic is the critical path: if it fails to strip markdown
code fences, json.loads() throws JSONDecodeError, the study is silently
dropped (returns None), and the evidence base loses a real study. If it
strips too aggressively, valid JSON content is lost.

ask_llm is monkeypatchable: requests.post is the only external call. We
replace it with a fake that returns canned responses to test the
fence-stripping and JSON parsing without network access.

Also includes a source-parsing test for the XAI_MODEL empty-string fix:
  PR #3 review fix: "empty XAI_MODEL env var must not override the grok-4.6 default"
  The fix changed from os.environ.get("XAI_MODEL", "grok-4.6") (which returns
  "" when the env var is set but empty) to os.environ.get("XAI_MODEL") or "grok-4.6"
  (which uses Python's truthiness: "" is falsy, so "or" falls through to the default).
"""
import unittest
import json
import os
import re
from unittest.mock import patch, MagicMock
import main


class FakeResponse:
    """Mimics requests.Response with .json() and .raise_for_status()."""
    def __init__(self, content_text, status_code=200):
        self._content_text = content_text
        self.status_code = status_code

    def raise_for_status(self):
        if self.status_code >= 400:
            raise Exception(f"HTTP {self.status_code}")

    def json(self):
        return {
            "choices": [{
                "message": {"content": self._content_text}
            }]
        }


class TestAskLlmFenceStripping(unittest.TestCase):

    def setUp(self):
        # ask_llm reads XAI_API_KEY at module scope; set it so the function
        # doesn't blow up trying to format the Authorization header with None
        self._orig_key = os.environ.get("XAI_API_KEY")
        os.environ["XAI_API_KEY"] = "fake-test-key"
        main.XAI_API_KEY = "fake-test-key"
        # Save and reset stats
        self._orig_stats = dict(main.stats)
        main.stats["errors"] = 0

    def tearDown(self):
        if self._orig_key is None:
            os.environ.pop("XAI_API_KEY", None)
        else:
            os.environ["XAI_API_KEY"] = self._orig_key
        main.XAI_API_KEY = os.environ.get("XAI_API_KEY")
        main.stats.clear()
        main.stats.update(self._orig_stats)

    def _make_article(self):
        return {
            "title": "Test Study",
            "abstract": "This is a test abstract about biomarkers.",
            "journal": "Aging Cell",
        }

    @patch("main.requests.post")
    def test_plain_json_no_fences(self, mock_post):
        """LLM returns plain JSON without markdown fences → parsed directly."""
        raw_json = '{"evidence_type": "RCT", "sample_size": "100"}'
        mock_post.return_value = FakeResponse(raw_json)
        result = main.ask_llm(self._make_article(), "DOM_AGING")
        self.assertIsNotNone(result)
        self.assertEqual(result["evidence_type"], "RCT")
        self.assertEqual(result["sample_size"], "100")
        self.assertEqual(main.stats["errors"], 0)

    @patch("main.requests.post")
    def test_json_with_markdown_code_fences(self, mock_post):
        """LLM wraps JSON in ```json ... ``` fences → fences stripped, JSON parsed."""
        raw = '```json\n{"evidence_type": "Meta-analysis", "sample_size": "500"}\n```'
        mock_post.return_value = FakeResponse(raw)
        result = main.ask_llm(self._make_article(), "DOM_AGING")
        self.assertIsNotNone(result, "Fenced JSON must be parsed successfully")
        self.assertEqual(result["evidence_type"], "Meta-analysis")
        self.assertEqual(result["sample_size"], "500")
        self.assertEqual(main.stats["errors"], 0)

    @patch("main.requests.post")
    def test_json_with_bare_code_fences_no_language(self, mock_post):
        """LLM wraps JSON in bare ``` fences (no 'json' prefix)."""
        raw = '```\n{"evidence_type": "Cross-sectional", "sample_size": "200"}\n```'
        mock_post.return_value = FakeResponse(raw)
        result = main.ask_llm(self._make_article(), "DOM_AGING")
        self.assertIsNotNone(result)
        self.assertEqual(result["evidence_type"], "Cross-sectional")

    @patch("main.requests.post")
    def test_inline_fence_no_newline_returns_none(self, mock_post):
        """Inline fences without a newline after the opening fence are NOT
        handled by the production code — text.startswith('```') triggers
        split('\\n', 1)[1], which raises IndexError when there's no newline.
        This is caught by the except Exception handler and returns None.
        Documenting this as the live contract (pitfall 40: test the actual
        behavior, not the hoped-for one)."""
        raw = '```json{"evidence_type": "Prospective cohort"}```'
        mock_post.return_value = FakeResponse(raw)
        result = main.ask_llm(self._make_article(), "DOM_AGING")
        self.assertIsNone(result, "Inline fence without newline is not parseable")
        self.assertEqual(main.stats["errors"], 1)

    @patch("main.requests.post")
    def test_malformed_json_returns_none_and_increments_errors(self, mock_post):
        """Malformed JSON → JSONDecodeError → returns None, stats['errors'] incremented."""
        mock_post.return_value = FakeResponse("{not valid json")
        result = main.ask_llm(self._make_article(), "DOM_AGING")
        self.assertIsNone(result, "Malformed JSON must return None")
        self.assertEqual(main.stats["errors"], 1)

    @patch("main.requests.post")
    def test_api_error_returns_none_and_increments_errors(self, mock_post):
        """API failure → Exception → returns None, stats['errors'] incremented."""
        mock_post.side_effect = Exception("Connection refused")
        result = main.ask_llm(self._make_article(), "DOM_AGING")
        self.assertIsNone(result)
        self.assertEqual(main.stats["errors"], 1)

    @patch("main.requests.post")
    def test_http_error_returns_none_and_increments_errors(self, mock_post):
        """HTTP error (raise_for_status) → Exception → returns None, errors++."""
        mock_post.return_value = FakeResponse("error", status_code=429)
        result = main.ask_llm(self._make_article(), "DOM_AGING")
        self.assertIsNone(result)
        self.assertEqual(main.stats["errors"], 1)

    @patch("main.requests.post")
    def test_sends_bearer_auth_with_xai_api_key(self, mock_post):
        """Authorization header must be 'Bearer <XAI_API_KEY>'."""
        mock_post.return_value = FakeResponse('{"evidence_type": "RCT"}')
        main.ask_llm(self._make_article(), "DOM_AGING")
        call_kwargs = mock_post.call_args.kwargs
        self.assertEqual(
            call_kwargs["headers"]["Authorization"],
            "Bearer fake-test-key",
            "xAI adapter must send Bearer token auth"
        )

    @patch("main.requests.post")
    def test_sends_model_from_xai_model_constant(self, mock_post):
        """The model in the request body must come from the XAI_MODEL constant."""
        mock_post.return_value = FakeResponse('{"evidence_type": "RCT"}')
        main.ask_llm(self._make_article(), "DOM_AGING")
        call_kwargs = mock_post.call_args.kwargs
        self.assertEqual(
            call_kwargs["json"]["model"],
            main.XAI_MODEL,
            "Request body model must match XAI_MODEL constant"
        )

    @patch("main.requests.post")
    def test_sends_to_xai_base_url(self, mock_post):
        """URL must be {XAI_BASE_URL}/chat/completions."""
        mock_post.return_value = FakeResponse('{"evidence_type": "RCT"}')
        main.ask_llm(self._make_article(), "DOM_AGING")
        call_args = mock_post.call_args.args
        self.assertIn("/chat/completions", call_args[0],
            "URL must include /chat/completions endpoint")
        self.assertIn(main.XAI_BASE_URL, call_args[0],
            "URL must start with XAI_BASE_URL")


class TestXaiModelDefault(unittest.TestCase):
    """
    Source-parsing test for the XAI_MODEL empty-string fix (PR #3 review).

    The fix changed from:
        os.environ.get("XAI_MODEL", "grok-4.6")   # returns "" if env is set but empty
    to:
        os.environ.get("XAI_MODEL") or "grok-4.6"  # "" is falsy → "grok-4.6"

    We source-parse (not import) because the constant is set at module scope
    from the environment, and we need to verify the EXPRESSION, not the
    resolved value (which depends on the test environment's env vars).

    Per pitfall 29: strip comments/docstrings before asserting, and pin the
    executable contract, not a vocabulary word.
    """

    def test_xai_model_uses_or_pattern_not_dict_default(self):
        """The `or` operator is what makes empty-string fall through to default.
        `os.environ.get("XAI_MODEL", "grok-4.6")` would return "" for an empty
        env var — the production bug that PR #3 review fixed."""
        main_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), "main.py")
        with open(main_path) as f:
            source = f.read()

        # Strip comments and docstrings (pitfall 29)
        # Remove # comments
        lines = []
        in_docstring = False
        for line in source.split("\n"):
            stripped = line.lstrip()
            if '"""' in stripped or "'''" in stripped:
                # Toggle docstring state (simplified: handles single-line and multi-line)
                count = stripped.count('"""') + stripped.count("'''")
                if count >= 2:
                    # Single-line docstring on one line — skip this line
                    continue
                in_docstring = not in_docstring
                continue
            if in_docstring:
                continue
            if stripped.startswith("#"):
                continue
            lines.append(line)
        executable = "\n".join(lines)

        # Pin the executable contract: `or "grok-4.6"` (not dict default)
        # This is the expression that makes empty string fall through.
        self.assertRegex(
            executable,
            r'XAI_MODEL\s*=\s*os\.environ\.get\(\s*["\']XAI_MODEL["\']\s*\)\s*or\s*["\']grok-4\.6["\']',
            "XAI_MODEL must use `os.environ.get('XAI_MODEL') or 'grok-4.6'` "
            "(the `or` pattern), NOT `os.environ.get('XAI_MODEL', 'grok-4.6')` "
            "(the dict-default pattern that returns empty string for an empty env var)"
        )

        # Negative check: the old buggy pattern must NOT be present
        self.assertNotRegex(
            executable,
            r'XAI_MODEL\s*=\s*os\.environ\.get\(\s*["\']XAI_MODEL["\']\s*,\s*["\']grok-4\.6["\']\s*\)',
            "The dict-default pattern os.environ.get('XAI_MODEL', 'grok-4.6') "
            "must NOT be used — it returns empty string for an empty env var"
        )

    def test_xai_base_url_uses_or_pattern(self):
        """XAI_BASE_URL also uses the `or` pattern for empty-string safety."""
        main_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), "main.py")
        with open(main_path) as f:
            source = f.read()

        # Strip comments (single-line only for this check — simpler)
        executable = "\n".join(
            line for line in source.split("\n")
            if not line.lstrip().startswith("#")
        )

        self.assertRegex(
            executable,
            r'XAI_BASE_URL\s*=\s*os\.environ\.get\(\s*["\']XAI_BASE_URL["\']\s*\)\s*or\s*["\']https://api\.x\.ai/v1["\']',
            "XAI_BASE_URL must use the `or` pattern for empty-string safety"
        )


if __name__ == "__main__":
    unittest.main()
