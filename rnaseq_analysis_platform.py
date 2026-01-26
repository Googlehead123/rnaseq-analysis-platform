"""
RNA-seq Analysis Platform
A Streamlit-based application for RNA-seq data analysis and visualization.
"""

import streamlit as st
import pandas as pd
from typing import Dict, List, Tuple

# Allowed group types (fixed enum)
ALLOWED_GROUPS = ["control", "positive_control", "treatment"]


def get_unique_conditions(sample_metadata: Dict[str, Dict]) -> List[str]:
    """Get list of unique conditions from metadata."""
    conditions = set(m["condition"] for m in sample_metadata.values())
    return sorted(list(conditions))


def count_samples_per_condition(sample_metadata: Dict[str, Dict]) -> Dict[str, int]:
    """Count how many samples are assigned to each condition."""
    counts = {}
    for meta in sample_metadata.values():
        condition = meta["condition"]
        counts[condition] = counts.get(condition, 0) + 1
    return counts


def validate_metadata_assignment(
    sample_metadata: Dict[str, Dict], counts_df: pd.DataFrame
) -> Tuple[bool, List[str]]:
    """
    Validate metadata assignment.

    Checks:
    - All samples assigned
    - Max 10 conditions (guardrail)
    - Min 2 samples per condition

    Returns:
        (is_valid, list_of_errors)
    """
    errors = []

    # Check all samples assigned
    unassigned = [s for s in counts_df.index if s not in sample_metadata]
    if unassigned:
        errors.append(f"Unassigned samples: {', '.join(unassigned)}")

    # Get unique conditions
    conditions = get_unique_conditions(sample_metadata)

    # Guardrail: Max 10 conditions
    if len(conditions) > 10:
        errors.append(
            "Maximum 10 conditions supported. Please reduce the number of unique conditions."
        )

    # Check minimum 2 samples per condition
    sample_counts = count_samples_per_condition(sample_metadata)
    for condition, count in sample_counts.items():
        if count < 2:
            errors.append(
                f"Condition '{condition}' has only {count} sample(s), need ≥2 for DE analysis"
            )

    return len(errors) == 0, errors


# Page configuration
st.set_page_config(
    page_title="RNA-seq Analysis Platform",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

# Main title
st.title("🧬 RNA-seq Analysis Platform")
st.markdown("---")

# Placeholder content
st.info("RNA-seq Analysis Platform - Coming Soon")
st.write(
    "This application will provide comprehensive RNA-seq data analysis capabilities."
)
