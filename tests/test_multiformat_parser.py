# pyright: reportMissingImports=false, reportArgumentType=false
"""Tests for multi-format ParseResult support."""

from pathlib import Path

import pandas as pd
from rnaseq_parser import (
    ParseResult,
    DataType,
    detect_de_columns,
    detect_sample_columns,
    SampleColumnInfo,
    RNASeqParser,
)


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
            index=pd.Index(["Gene_A", "Gene_B", "Gene_C"]),
        )

        normalized_df = pd.DataFrame(
            {
                "Sample_1": [10.5, 20.3, 30.1],
                "Sample_2": [15.2, 25.4, 35.6],
            },
            index=pd.Index(["Gene_A", "Gene_B", "Gene_C"]),
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
            index=pd.Index(["Gene_A", "Gene_B", "Gene_C"]),
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


class TestDetectDEColumns:
    """Test pattern-based DE column detection."""

    def test_detect_de_columns_collision_raises_error(self):
        """Test that multiple comparisons with different prefixes raise error."""
        import pytest
        from rnaseq_parser import ParserValidationError

        # DataFrame with columns from two different comparisons
        df = pd.DataFrame(
            {
                "gene": ["Gene_A", "Gene_B", "Gene_C"],
                "Bt10U/none.fc": [2.5, -1.3, 0.8],
                "Bt10U/other.fc": [1.5, -0.8, 0.5],  # Different comparison
                "Bt10U/none.bh.pval": [0.001, 0.05, 0.1],
                "Bt10U/other.bh.pval": [0.002, 0.06, 0.15],  # Different comparison
            }
        )

        # Should raise ParserValidationError due to prefix mismatch
        with pytest.raises(ParserValidationError) as exc_info:
            detect_de_columns(df)

        assert "prefix mismatch" in str(exc_info.value).lower()
        assert "multiple comparisons" in str(exc_info.value).lower()

    def test_detect_de_columns_pattern_fc(self):
        """Test detection of fold change column with pattern matching."""
        df = pd.DataFrame(
            {
                "gene": ["Gene_A", "Gene_B", "Gene_C"],
                "Bt10U/none.fc": [2.5, -1.3, 0.8],
                "Bt10U/none.bh.pval": [0.001, 0.05, 0.1],
            }
        )

        column_mapping, comparison_name = detect_de_columns(df)

        assert "log2FoldChange" in column_mapping
        assert column_mapping["log2FoldChange"] == "Bt10U/none.fc"
        assert comparison_name == "Bt10U/none"

    def test_detect_de_columns_pattern_pval(self):
        """Test detection of adjusted p-value column with pattern matching."""
        df = pd.DataFrame(
            {
                "gene": ["Gene_A", "Gene_B", "Gene_C"],
                "Bt10U/none.fc": [2.5, -1.3, 0.8],
                "Bt10U/none.bh.pval": [0.001, 0.05, 0.1],
            }
        )

        column_mapping, comparison_name = detect_de_columns(df)

        assert "padj" in column_mapping
        assert column_mapping["padj"] == "Bt10U/none.bh.pval"

    def test_detect_de_columns_extracts_comparison_name(self):
        """Test extraction of comparison name from pattern."""
        df = pd.DataFrame(
            {
                "gene": ["Gene_A", "Gene_B", "Gene_C"],
                "Bt10U/none.fc": [2.5, -1.3, 0.8],
                "Bt10U/none.raw.pval": [0.001, 0.05, 0.1],
            }
        )

        column_mapping, comparison_name = detect_de_columns(df)

        assert comparison_name == "Bt10U/none"

    def test_detect_de_columns_fallback_to_standard(self):
        """Test fallback to standard column names when patterns don't match."""
        df = pd.DataFrame(
            {
                "gene": ["Gene_A", "Gene_B", "Gene_C"],
                "log2FoldChange": [2.5, -1.3, 0.8],
                "padj": [0.001, 0.05, 0.1],
            }
        )

        column_mapping, comparison_name = detect_de_columns(df)

        assert column_mapping["log2FoldChange"] == "log2FoldChange"
        assert column_mapping["padj"] == "padj"
        assert comparison_name is None


class TestDetectSampleColumns:
    """Test sample/condition auto-detection from column names."""

    def test_detect_sample_columns_read_count(self):
        """Test detection of Read_Count columns."""
        df = pd.DataFrame(
            {
                "231222_none_3d_Read_Count": [100, 200, 300],
                "231223_none_3d_Read_Count": [150, 250, 350],
            },
            index=pd.Index(["Gene_A", "Gene_B", "Gene_C"]),
        )

        result = detect_sample_columns(df)

        assert isinstance(result, SampleColumnInfo)
        assert len(result.count_columns) == 2
        assert "231222_none_3d_Read_Count" in result.count_columns
        assert "231223_none_3d_Read_Count" in result.count_columns

    def test_detect_sample_columns_fpkm(self):
        """Test detection of FPKM columns."""
        df = pd.DataFrame(
            {
                "240203_Bt10U_3d_FPKM": [10.5, 20.3, 30.1],
                "240204_Bt10U_3d_FPKM": [15.2, 25.4, 35.6],
            },
            index=pd.Index(["Gene_A", "Gene_B", "Gene_C"]),
        )

        result = detect_sample_columns(df)

        assert isinstance(result, SampleColumnInfo)
        assert len(result.fpkm_columns) == 2
        assert "240203_Bt10U_3d_FPKM" in result.fpkm_columns
        assert "240204_Bt10U_3d_FPKM" in result.fpkm_columns

    def test_detect_sample_columns_tpm(self):
        """Test detection of TPM columns."""
        df = pd.DataFrame(
            {
                "240203_Bt10U_3d_TPM": [10.5, 20.3, 30.1],
                "240204_Bt10U_3d_TPM": [15.2, 25.4, 35.6],
            },
            index=pd.Index(["Gene_A", "Gene_B", "Gene_C"]),
        )

        result = detect_sample_columns(df)

        assert isinstance(result, SampleColumnInfo)
        assert len(result.tpm_columns) == 2
        assert "240203_Bt10U_3d_TPM" in result.tpm_columns
        assert "240204_Bt10U_3d_TPM" in result.tpm_columns

    def test_detect_conditions_from_columns(self):
        """Test extraction of unique conditions from sample columns."""
        df = pd.DataFrame(
            {
                "231222_none_3d_Read_Count": [100, 200, 300],
                "231223_none_3d_Read_Count": [150, 250, 350],
                "240203_Bt10U_3d_FPKM": [10.5, 20.3, 30.1],
                "240204_Bt10U_3d_FPKM": [15.2, 25.4, 35.6],
            },
            index=pd.Index(["Gene_A", "Gene_B", "Gene_C"]),
        )

        result = detect_sample_columns(df)

        assert isinstance(result, SampleColumnInfo)
        assert "none" in result.conditions_detected
        assert "Bt10U" in result.conditions_detected
        assert len(result.conditions_detected) == 2

    def test_sample_to_condition_mapping(self):
        """Test mapping of sample columns to conditions."""
        df = pd.DataFrame(
            {
                "231222_none_3d_Read_Count": [100, 200, 300],
                "240203_Bt10U_3d_FPKM": [10.5, 20.3, 30.1],
            },
            index=pd.Index(["Gene_A", "Gene_B", "Gene_C"]),
        )

        result = detect_sample_columns(df)

        assert isinstance(result, SampleColumnInfo)
        assert result.sample_to_condition["231222_none_3d_Read_Count"] == "none"
        assert result.sample_to_condition["240203_Bt10U_3d_FPKM"] == "Bt10U"

    def test_detect_sample_columns_with_various_conditions(self):
        """Test detection of sample columns with various condition names."""
        df = pd.DataFrame(
            {
                "231222_none_3d_Read_Count": [100, 200, 300],
                "231222_Bt10U_3d_FPKM": [10.5, 20.3, 30.1],
                "240203_SwX15_5d_TPM": [5.2, 8.9, 12.4],
                "240726_my45_7d_Read_Count": [150, 250, 350],
            },
            index=pd.Index(["Gene1", "Gene2", "Gene3"]),
        )

        result = detect_sample_columns(df)

        assert isinstance(result, SampleColumnInfo)
        assert "none" in result.conditions_detected
        assert "Bt10U" in result.conditions_detected
        assert "SwX15" in result.conditions_detected
        assert "my45" in result.conditions_detected
        assert "3d" not in result.conditions_detected
        assert "5d" not in result.conditions_detected
        assert "7d" not in result.conditions_detected

        assert result.sample_to_condition["231222_none_3d_Read_Count"] == "none"
        assert result.sample_to_condition["231222_Bt10U_3d_FPKM"] == "Bt10U"
        assert result.sample_to_condition["240203_SwX15_5d_TPM"] == "SwX15"
        assert result.sample_to_condition["240726_my45_7d_Read_Count"] == "my45"


class TestParseMultiFormatExcel:
    """Integration tests for multi-format Excel parsing."""

    @staticmethod
    def _multiformat_path() -> str:
        return str(
            Path(__file__).resolve().parents[1]
            / "Reference sequencing data"
            / "data3_Bt10U_vs_none_fc2_&_raw.p.xlsx"
        )

    def test_parse_multiformat_extracts_de_results(self):
        """DE results are extracted with canonical columns."""
        parser = RNASeqParser()
        result = parser.parse_multiformat(self._multiformat_path())  # type: ignore[reportArgumentType]

        assert result.de_results_df is not None
        assert 300 <= len(result.de_results_df) <= 340
        assert set(["gene", "log2FoldChange", "padj"]).issubset(
            result.de_results_df.columns
        )

    def test_parse_multiformat_extracts_counts(self):
        """Read_Count columns are extracted into expression matrix."""
        parser = RNASeqParser()
        result = parser.parse_multiformat(self._multiformat_path())  # type: ignore[reportArgumentType]

        assert result.expression_df is not None
        assert result.expression_df.shape[0] == 6
        assert 300 <= result.expression_df.shape[1] <= 340

    def test_parse_multiformat_extracts_normalized(self):
        """TPM columns are extracted into normalized matrix."""
        parser = RNASeqParser()
        result = parser.parse_multiformat(self._multiformat_path())  # type: ignore[reportArgumentType]

        assert result.normalized_df is not None
        assert result.normalized_df.shape[0] == 6
        assert 300 <= result.normalized_df.shape[1] <= 340

    def test_parse_multiformat_all_datatypes_populated(self):
        """All three datatypes are populated with correct metadata."""
        parser = RNASeqParser()
        result = parser.parse_multiformat(self._multiformat_path())  # type: ignore[reportArgumentType]

        assert result.expression_df is not None
        assert result.de_results_df is not None
        assert result.normalized_df is not None
        assert result.data_types_detected == [
            DataType.PRE_ANALYZED,
            DataType.RAW_COUNTS,
            DataType.NORMALIZED,
        ]

    def test_parse_multiformat_backward_compatible(self):
        """Parser remains backward compatible for simple CSV files."""
        parser = RNASeqParser()
        csv_path = (
            Path(__file__).resolve().parents[1] / "tests" / "data" / "sample_counts.csv"
        )
        result = parser.parse(str(csv_path))

        assert result.expression_df is not None
        assert result.normalized_df is None
        assert result.de_results_df is None
        assert result.data_types_detected == []


def test_parse_multiformat_extracts_de_results():
    """Test that DE results are extracted with canonical column names."""
    parser = RNASeqParser()
    result = parser.parse(
        "Reference sequencing data/data3_Bt10U_vs_none_fc2_&_raw.p.xlsx"
    )

    assert result.de_results_df is not None, "DE results should be extracted"
    assert "gene" in result.de_results_df.columns, "Should have 'gene' column"
    assert "log2FoldChange" in result.de_results_df.columns, (
        "Should have canonical log2FoldChange"
    )
    assert "padj" in result.de_results_df.columns, "Should have canonical padj"
    assert len(result.de_results_df) > 0, "Should have genes"


def test_parse_multiformat_extracts_counts():
    """Test that count matrix is extracted from Read_Count columns."""
    parser = RNASeqParser()
    result = parser.parse(
        "Reference sequencing data/data3_Bt10U_vs_none_fc2_&_raw.p.xlsx"
    )

    assert result.expression_df is not None, "Expression matrix should be extracted"
    assert result.expression_df.shape[0] == 6, "Should have 6 samples (rows)"
    assert result.expression_df.shape[1] > 0, "Should have genes (columns)"


def test_parse_multiformat_extracts_normalized():
    """Test that normalized matrix is extracted from TPM columns."""
    parser = RNASeqParser()
    result = parser.parse(
        "Reference sequencing data/data3_Bt10U_vs_none_fc2_&_raw.p.xlsx"
    )

    assert result.normalized_df is not None, "Normalized matrix should be extracted"
    assert result.normalized_df.shape[0] == 6, "Should have 6 samples (rows)"
    assert result.normalized_df.shape[1] > 0, "Should have genes (columns)"


def test_parse_multiformat_all_datatypes_populated():
    """Test that all three data types are detected and populated."""
    parser = RNASeqParser()
    result = parser.parse(
        "Reference sequencing data/data3_Bt10U_vs_none_fc2_&_raw.p.xlsx"
    )

    assert DataType.PRE_ANALYZED in result.data_types_detected
    assert DataType.RAW_COUNTS in result.data_types_detected
    assert DataType.NORMALIZED in result.data_types_detected
    assert len(result.data_types_detected) == 3


def test_parse_multiformat_backward_compatible():
    """Test that simple CSV files still work as before."""
    parser = RNASeqParser()
    csv_path = (
        Path(__file__).resolve().parents[1] / "tests" / "data" / "sample_counts.csv"
    )
    result = parser.parse(str(csv_path))

    assert result.expression_df is not None
    assert result.normalized_df is None
    assert result.de_results_df is None
    assert result.data_types_detected == []
