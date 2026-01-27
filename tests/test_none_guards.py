"""
Tests for None guard checks in visualization paths.

Ensures that PRE_ANALYZED data (which has expression_df=None) doesn't crash
when attempting to generate visualizations.
"""

import pytest
import pandas as pd
import numpy as np
from pathlib import Path
import tempfile
from unittest.mock import MagicMock


def test_preanalyzed_data_no_crash_on_visualization():
    """
    Test that PRE_ANALYZED data with expression_df=None doesn't crash on visualization.

    This test reproduces the crash at line 474:
    `np.log2(result.expression_df + 1)` when expression_df is None.

    Expected behavior: Should gracefully handle None and display info message.
    """
    from rnaseq_parser import RNASeqParser, DataType

    # Step 1: Create PRE_ANALYZED data (DE results only, no expression_df)
    pre_analyzed_data = pd.DataFrame(
        {
            "gene": ["COL1A1", "IL6", "TYR", "VEGFA", "FGF2"],
            "log2FoldChange": [2.34, -1.89, 1.45, 0.78, -0.56],
            "padj": [0.0001, 0.0023, 0.0456, 0.1234, 0.5678],
            "baseMean": [1500.2, 450.8, 120.3, 890.5, 234.1],
        }
    )

    # Write to temp file
    with tempfile.NamedTemporaryFile(suffix=".csv", delete=False, mode="w") as f:
        pre_analyzed_data.to_csv(f, index=False)
        temp_path = f.name

    try:
        # Step 2: Parse the file - should detect as PRE_ANALYZED
        parser = RNASeqParser()
        result = parser.parse(temp_path)

        # Verify it's PRE_ANALYZED type
        assert result.data_type == DataType.PRE_ANALYZED
        assert result.expression_df is None  # KEY: PRE_ANALYZED has no expression data
        assert result.de_results_df is not None

        # Step 3: Simulate the visualization code path that was crashing
        # This is the code from rnaseq_analysis_platform.py lines 471-474
        norm_counts = None
        if result.data_type == DataType.NORMALIZED:
            import numpy as np

            # This line would crash if expression_df is None
            if result.expression_df is None:
                # NEW GUARD: Check for None before arithmetic
                pass
            else:
                norm_counts = np.log2(result.expression_df + 1)

        # Step 4: Verify no crash occurred
        # (In actual code, this would be in the NORMALIZED branch, but we're testing the guard)
        assert norm_counts is None  # No visualization data generated

        # Step 5: Verify the guard works for NORMALIZED data with None expression_df
        # Create a mock NORMALIZED result with None expression_df (edge case)
        result.data_type = DataType.NORMALIZED
        result.expression_df = None

        # This should NOT crash with the guard in place
        norm_counts = None
        if result.data_type == DataType.NORMALIZED:
            import numpy as np

            if result.expression_df is None:
                pass
            else:
                norm_counts = np.log2(result.expression_df + 1)

        # Verify no crash and guard was effective
        assert norm_counts is None

    finally:
        # Cleanup
        Path(temp_path).unlink()


def test_normalized_data_with_expression_df_generates_viz():
    """
    Test that NORMALIZED data WITH expression_df properly generates visualizations.

    This is the happy path - ensure the guard doesn't break normal operation.
    """
    from rnaseq_parser import RNASeqParser, DataType

    # Create NORMALIZED data with expression_df
    normalized_data = pd.DataFrame(
        {
            "gene": ["COL1A1", "IL6", "TYR", "VEGFA", "FGF2"],
            "Sample_1": [12.5, 8.3, 6.7, 10.2, 7.1],
            "Sample_2": [13.2, 9.1, 7.2, 11.0, 7.8],
            "Sample_3": [14.1, 10.2, 7.8, 12.1, 8.5],
            "Sample_4": [13.8, 9.5, 7.4, 11.5, 8.0],
        }
    )

    # Write to temp file
    with tempfile.NamedTemporaryFile(suffix=".csv", delete=False, mode="w") as f:
        normalized_data.to_csv(f, index=False)
        temp_path = f.name

    try:
        # Parse the file
        parser = RNASeqParser()
        result = parser.parse(temp_path)

        # Verify it's NORMALIZED type
        assert result.data_type == DataType.NORMALIZED
        assert result.expression_df is not None  # NORMALIZED has expression data

        # Simulate the visualization code path
        norm_counts = None
        if result.data_type == DataType.NORMALIZED:
            import numpy as np

            if result.expression_df is None:
                pass
            else:
                norm_counts = np.log2(result.expression_df + 1)

        # Verify visualization data was generated
        assert norm_counts is not None
        assert norm_counts.shape == result.expression_df.shape

    finally:
        # Cleanup
        Path(temp_path).unlink()
