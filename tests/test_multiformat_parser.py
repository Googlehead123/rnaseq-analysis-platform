"""Tests for multi-format ParseResult support."""

import pandas as pd
import pytest
from rnaseq_parser import ParseResult, DataType


class TestParseResultMultiFormat:
    """Test ParseResult with multiple data types."""

    def test_parseresult_supports_multiple_datatypes(self):
        """Test that ParseResult can hold all three data types simultaneously."""
        # Create sample DataFrames
        expression_df = pd.DataFrame(
            {
                "Sample_1": [100, 200, 300],
                "Sample_2": [150, 250, 350],
            },
            index=["Gene_A", "Gene_B", "Gene_C"],
        )

        normalized_df = pd.DataFrame(
            {
                "Sample_1": [10.5, 20.3, 30.1],
                "Sample_2": [15.2, 25.4, 35.6],
            },
            index=["Gene_A", "Gene_B", "Gene_C"],
        )

        de_results_df = pd.DataFrame(
            {
                "gene": ["Gene_A", "Gene_B", "Gene_C"],
                "log2FoldChange": [2.5, -1.3, 0.8],
                "padj": [0.001, 0.05, 0.1],
                "baseMean": [150.0, 225.0, 325.0],
            }
        )

        # Create ParseResult with all three fields populated
        result = ParseResult(
            expression_df=expression_df,
            de_results_df=de_results_df,
            data_type=DataType.RAW_COUNTS,
            can_run_de=True,
            warnings=[],
            dropped_columns=[],
            gene_column_source="test",
            needs_user_input=False,
            gene_column_candidates=[],
            normalized_df=normalized_df,
            data_types_detected=[
                DataType.RAW_COUNTS,
                DataType.NORMALIZED,
                DataType.PRE_ANALYZED,
            ],
        )

        # Verify all fields are accessible
        assert result.expression_df is not None
        assert result.expression_df.shape == (3, 2)
        assert result.normalized_df is not None
        assert result.normalized_df.shape == (3, 2)
        assert result.de_results_df is not None
        assert result.de_results_df.shape == (3, 4)
        assert len(result.data_types_detected) == 3
        assert DataType.RAW_COUNTS in result.data_types_detected
        assert DataType.NORMALIZED in result.data_types_detected
        assert DataType.PRE_ANALYZED in result.data_types_detected

    def test_parseresult_backward_compatible(self):
        """Test that ParseResult maintains backward compatibility with old behavior."""
        # Create sample DataFrame
        expression_df = pd.DataFrame(
            {
                "Sample_1": [100, 200, 300],
                "Sample_2": [150, 250, 350],
            },
            index=["Gene_A", "Gene_B", "Gene_C"],
        )

        # Create ParseResult with only expression_df (old behavior)
        result = ParseResult(
            expression_df=expression_df,
            de_results_df=None,
            data_type=DataType.RAW_COUNTS,
            can_run_de=True,
            warnings=[],
            dropped_columns=[],
            gene_column_source="test",
            needs_user_input=False,
            gene_column_candidates=[],
            normalized_df=None,
            data_types_detected=[],
        )

        # Verify backward compatibility
        assert result.expression_df is not None
        assert result.expression_df.shape == (3, 2)
        assert result.normalized_df is None
        assert result.de_results_df is None
        assert result.data_types_detected == []
        assert result.data_type == DataType.RAW_COUNTS
        assert result.can_run_de is True
