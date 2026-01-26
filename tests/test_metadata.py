"""
Tests for metadata validation functions from rnaseq_analysis_platform.py

Tests cover:
- get_unique_conditions()
- count_samples_per_condition()
- validate_metadata_assignment()
"""

import pytest
import pandas as pd
from rnaseq_analysis_platform import (
    get_unique_conditions,
    count_samples_per_condition,
    validate_metadata_assignment,
)


# ============================================================================
# Tests for get_unique_conditions()
# ============================================================================


def test_get_unique_conditions():
    """Extract unique conditions from metadata."""
    metadata = {
        "S1": {"condition": "Control", "group": "control"},
        "S2": {"condition": "Control", "group": "control"},
        "S3": {"condition": "Treatment", "group": "treatment"},
    }
    conditions = get_unique_conditions(metadata)
    assert conditions == ["Control", "Treatment"]  # Sorted


def test_get_unique_conditions_empty():
    """Handle empty metadata."""
    metadata = {}
    conditions = get_unique_conditions(metadata)
    assert conditions == []


def test_get_unique_conditions_single():
    """Handle single condition."""
    metadata = {
        "S1": {"condition": "Control", "group": "control"},
        "S2": {"condition": "Control", "group": "control"},
    }
    conditions = get_unique_conditions(metadata)
    assert conditions == ["Control"]


# ============================================================================
# Tests for count_samples_per_condition()
# ============================================================================


def test_count_samples_per_condition():
    """Count samples per condition."""
    metadata = {
        "S1": {"condition": "Control", "group": "control"},
        "S2": {"condition": "Control", "group": "control"},
        "S3": {"condition": "Treatment", "group": "treatment"},
    }
    counts = count_samples_per_condition(metadata)
    assert counts == {"Control": 2, "Treatment": 1}


def test_count_samples_per_condition_empty():
    """Handle empty metadata."""
    metadata = {}
    counts = count_samples_per_condition(metadata)
    assert counts == {}


def test_count_samples_per_condition_multiple_conditions():
    """Count samples across multiple conditions."""
    metadata = {
        f"S{i}": {"condition": f"Cond_{i % 3}", "group": "test"} for i in range(1, 10)
    }
    counts = count_samples_per_condition(metadata)
    assert counts["Cond_1"] == 3
    assert counts["Cond_2"] == 3
    assert counts["Cond_0"] == 3


# ============================================================================
# Tests for validate_metadata_assignment()
# ============================================================================


def test_validate_metadata_min_samples(sample_counts_df):
    """Fail when condition has < 2 samples."""
    metadata = {"sample_1": {"condition": "Control", "group": "control"}}
    is_valid, errors = validate_metadata_assignment(metadata, sample_counts_df)
    assert not is_valid
    assert any("need ≥2" in error for error in errors)


def test_validate_metadata_max_conditions(sample_counts_df):
    """Fail when > 10 unique conditions."""
    # Create metadata with 11 unique conditions
    metadata = {
        f"sample_{i}": {"condition": f"Cond_{i}", "group": "treatment"}
        for i in range(1, 12)
    }
    is_valid, errors = validate_metadata_assignment(metadata, sample_counts_df)
    assert not is_valid
    assert any("Maximum 10 conditions" in error for error in errors)


def test_validate_metadata_success(sample_counts_df):
    """Valid assignment passes all checks."""
    metadata = {
        f"sample_{i}": {
            "condition": "Control" if i <= 5 else "Treatment",
            "group": "control" if i <= 5 else "treatment",
        }
        for i in range(1, 11)
    }
    is_valid, errors = validate_metadata_assignment(metadata, sample_counts_df)
    assert is_valid
    assert errors == []


def test_validate_metadata_unassigned_samples(sample_counts_df):
    """Fail when samples are unassigned."""
    # Only assign 5 samples, but sample_counts_df has 10
    metadata = {
        f"sample_{i}": {"condition": "Control", "group": "control"} for i in range(1, 6)
    }
    is_valid, errors = validate_metadata_assignment(metadata, sample_counts_df)
    assert not is_valid
    assert any("Unassigned samples" in error for error in errors)


def test_validate_metadata_multiple_errors(sample_counts_df):
    """Detect multiple validation errors."""
    # Only 1 sample assigned (unassigned + min samples error)
    metadata = {"sample_1": {"condition": "Control", "group": "control"}}
    is_valid, errors = validate_metadata_assignment(metadata, sample_counts_df)
    assert not is_valid
    assert len(errors) >= 2  # Unassigned + min samples


def test_validate_metadata_exactly_two_per_condition(sample_counts_df):
    """Pass with exactly 2 samples per condition (boundary case)."""
    metadata = {
        "sample_1": {"condition": "Control", "group": "control"},
        "sample_2": {"condition": "Control", "group": "control"},
        "sample_3": {"condition": "Treatment", "group": "treatment"},
        "sample_4": {"condition": "Treatment", "group": "treatment"},
        "sample_5": {"condition": "Other", "group": "other"},
        "sample_6": {"condition": "Other", "group": "other"},
        "sample_7": {"condition": "A", "group": "a"},
        "sample_8": {"condition": "A", "group": "a"},
        "sample_9": {"condition": "B", "group": "b"},
        "sample_10": {"condition": "B", "group": "b"},
    }
    is_valid, errors = validate_metadata_assignment(metadata, sample_counts_df)
    assert is_valid
    assert errors == []


def test_validate_metadata_exactly_ten_conditions(sample_counts_df):
    """Pass with exactly 10 conditions (boundary case)."""
    metadata = {
        f"sample_{i}": {
            "condition": f"Cond_{i}",
            "group": f"group_{i}",
        }
        for i in range(1, 11)
    }
    is_valid, errors = validate_metadata_assignment(metadata, sample_counts_df)
    assert is_valid
    assert errors == []
