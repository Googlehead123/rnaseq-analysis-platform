"""Tests for comparison validation functions from rnaseq_analysis_platform.py"""

import pytest
from unittest.mock import MagicMock, patch
from rnaseq_analysis_platform import (
    validate_comparison,
    add_comparison,
    remove_comparison,
)


@pytest.fixture
def mock_streamlit():
    """Mock Streamlit session state for testing."""
    with patch("rnaseq_analysis_platform.st") as mock_st:
        mock_st.session_state = {"comparisons": []}
        yield mock_st


class TestValidateComparison:
    """Tests for validate_comparison function."""

    def test_validate_comparison_same_condition(self):
        """Test that same test and ref conditions are rejected."""
        is_valid, err = validate_comparison(
            ("Control", "Control"), ["Control", "Treatment"]
        )
        assert not is_valid
        assert "must be different" in err.lower()

    def test_validate_comparison_missing_test_condition(self):
        """Test that missing test condition is rejected."""
        is_valid, err = validate_comparison(
            ("Missing", "Control"), ["Control", "Treatment"]
        )
        assert not is_valid
        assert "not valid" in err.lower()

    def test_validate_comparison_missing_ref_condition(self):
        """Test that missing ref condition is rejected."""
        is_valid, err = validate_comparison(
            ("Treatment", "Missing"), ["Control", "Treatment"]
        )
        assert not is_valid
        assert "not valid" in err.lower()

    def test_validate_comparison_success(self):
        """Test that valid comparison passes validation."""
        is_valid, err = validate_comparison(
            ("Treatment", "Control"), ["Control", "Treatment"]
        )
        assert is_valid
        assert err == ""

    def test_validate_comparison_reverse_order(self):
        """Test that comparison order doesn't matter for validation."""
        is_valid, err = validate_comparison(
            ("Control", "Treatment"), ["Control", "Treatment"]
        )
        assert is_valid
        assert err == ""

    def test_validate_comparison_with_multiple_conditions(self):
        """Test validation with more than 2 conditions."""
        is_valid, err = validate_comparison(
            ("Treated_A", "Control"), ["Control", "Treated_A", "Treated_B", "Treated_C"]
        )
        assert is_valid
        assert err == ""

    def test_validate_comparison_empty_conditions_list(self):
        """Test validation with empty available conditions."""
        is_valid, err = validate_comparison(("Treatment", "Control"), [])
        assert not is_valid
        assert "not valid" in err.lower()


class TestAddComparison:
    """Tests for add_comparison function."""

    def test_add_comparison_new(self, mock_streamlit):
        """Test adding a new comparison to session state."""
        with patch("rnaseq_analysis_platform.st", mock_streamlit):
            add_comparison("Treatment", "Control")
            assert ("Treatment", "Control") in mock_streamlit.session_state[
                "comparisons"
            ]

    def test_add_comparison_no_duplicates(self, mock_streamlit):
        """Test that duplicate comparisons are not added."""
        with patch("rnaseq_analysis_platform.st", mock_streamlit):
            add_comparison("Treatment", "Control")
            add_comparison("Treatment", "Control")
            # Should only have one entry
            assert (
                mock_streamlit.session_state["comparisons"].count(
                    ("Treatment", "Control")
                )
                == 1
            )

    def test_add_comparison_multiple_different(self, mock_streamlit):
        """Test adding multiple different comparisons."""
        with patch("rnaseq_analysis_platform.st", mock_streamlit):
            add_comparison("Treatment", "Control")
            add_comparison("Treated_A", "Control")
            add_comparison("Treated_B", "Control")

            assert len(mock_streamlit.session_state["comparisons"]) == 3
            assert ("Treatment", "Control") in mock_streamlit.session_state[
                "comparisons"
            ]
            assert ("Treated_A", "Control") in mock_streamlit.session_state[
                "comparisons"
            ]
            assert ("Treated_B", "Control") in mock_streamlit.session_state[
                "comparisons"
            ]

    def test_add_comparison_preserves_existing(self, mock_streamlit):
        """Test that adding comparison preserves existing ones."""
        with patch("rnaseq_analysis_platform.st", mock_streamlit):
            mock_streamlit.session_state["comparisons"] = [("Old", "Ref")]
            add_comparison("Treatment", "Control")

            assert len(mock_streamlit.session_state["comparisons"]) == 2
            assert ("Old", "Ref") in mock_streamlit.session_state["comparisons"]
            assert ("Treatment", "Control") in mock_streamlit.session_state[
                "comparisons"
            ]


class TestRemoveComparison:
    """Tests for remove_comparison function."""

    def test_remove_comparison_existing(self, mock_streamlit):
        """Test removing an existing comparison."""
        with patch("rnaseq_analysis_platform.st", mock_streamlit):
            mock_streamlit.session_state["comparisons"] = [("Treatment", "Control")]
            remove_comparison(("Treatment", "Control"))

            assert ("Treatment", "Control") not in mock_streamlit.session_state[
                "comparisons"
            ]
            assert len(mock_streamlit.session_state["comparisons"]) == 0

    def test_remove_comparison_nonexistent(self, mock_streamlit):
        """Test removing a non-existent comparison (should not error)."""
        with patch("rnaseq_analysis_platform.st", mock_streamlit):
            mock_streamlit.session_state["comparisons"] = [("Treatment", "Control")]
            # Should not raise error
            remove_comparison(("Missing", "Ref"))

            # Original comparison should still be there
            assert ("Treatment", "Control") in mock_streamlit.session_state[
                "comparisons"
            ]

    def test_remove_comparison_from_multiple(self, mock_streamlit):
        """Test removing one comparison from multiple."""
        with patch("rnaseq_analysis_platform.st", mock_streamlit):
            mock_streamlit.session_state["comparisons"] = [
                ("Treatment", "Control"),
                ("Treated_A", "Control"),
                ("Treated_B", "Control"),
            ]
            remove_comparison(("Treated_A", "Control"))

            assert len(mock_streamlit.session_state["comparisons"]) == 2
            assert ("Treatment", "Control") in mock_streamlit.session_state[
                "comparisons"
            ]
            assert ("Treated_A", "Control") not in mock_streamlit.session_state[
                "comparisons"
            ]
            assert ("Treated_B", "Control") in mock_streamlit.session_state[
                "comparisons"
            ]

    def test_remove_comparison_empty_list(self, mock_streamlit):
        """Test removing from empty comparisons list."""
        with patch("rnaseq_analysis_platform.st", mock_streamlit):
            mock_streamlit.session_state["comparisons"] = []
            # Should not raise error
            remove_comparison(("Treatment", "Control"))

            assert len(mock_streamlit.session_state["comparisons"]) == 0


class TestComparisonWorkflow:
    """Integration tests for comparison workflow."""

    def test_add_and_remove_workflow(self, mock_streamlit):
        """Test typical add/remove workflow."""
        with patch("rnaseq_analysis_platform.st", mock_streamlit):
            # Add comparisons
            add_comparison("Treatment", "Control")
            add_comparison("Treated_A", "Control")
            assert len(mock_streamlit.session_state["comparisons"]) == 2

            # Remove one
            remove_comparison(("Treatment", "Control"))
            assert len(mock_streamlit.session_state["comparisons"]) == 1
            assert ("Treated_A", "Control") in mock_streamlit.session_state[
                "comparisons"
            ]

    def test_validate_before_add(self, mock_streamlit):
        """Test validating comparison before adding."""
        with patch("rnaseq_analysis_platform.st", mock_streamlit):
            available = ["Control", "Treatment"]

            # Validate before adding
            is_valid, err = validate_comparison(("Treatment", "Control"), available)
            assert is_valid

            # Add if valid
            if is_valid:
                add_comparison("Treatment", "Control")

            assert ("Treatment", "Control") in mock_streamlit.session_state[
                "comparisons"
            ]

    def test_validate_reject_before_add(self, mock_streamlit):
        """Test that invalid comparisons are not added."""
        with patch("rnaseq_analysis_platform.st", mock_streamlit):
            available = ["Control", "Treatment"]

            # Validate before adding
            is_valid, err = validate_comparison(("Control", "Control"), available)
            assert not is_valid

            # Don't add if invalid
            if is_valid:
                add_comparison("Control", "Control")

            assert ("Control", "Control") not in mock_streamlit.session_state[
                "comparisons"
            ]
