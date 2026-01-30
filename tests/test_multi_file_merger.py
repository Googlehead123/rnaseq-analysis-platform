"""Tests for multi-file merger utility."""

import pytest
import pandas as pd
from rnaseq_parser import ParseResult, DataType
from multi_file_merger import merge_parse_results, validate_merge_compatibility


def _make_result(df, data_type=DataType.RAW_COUNTS):
    """Helper to create ParseResult from DataFrame."""
    return ParseResult(
        expression_df=df,
        normalized_df=None,
        de_results_df=None,
        data_type=data_type,
        can_run_de=(data_type == DataType.RAW_COUNTS),
        warnings=[],
        dropped_columns=[],
        gene_column_source="test",
        needs_user_input=False,
        gene_column_candidates=[],
        data_types_detected=[data_type],
    )


class TestValidateMergeCompatibility:
    def test_compatible_files(self):
        df1 = pd.DataFrame({"A": [1, 2], "B": [3, 4]}, index=["S1", "S2"])
        df2 = pd.DataFrame({"A": [5, 6], "B": [7, 8]}, index=["S3", "S4"])
        result = validate_merge_compatibility([_make_result(df1), _make_result(df2)])
        assert result["compatible"] is True
        assert result["report"]["gene_overlap_pct"] == 100.0

    def test_mixed_data_types(self):
        df1 = pd.DataFrame({"A": [1]}, index=["S1"])
        df2 = pd.DataFrame({"A": [1]}, index=["S2"])
        r1 = _make_result(df1, DataType.RAW_COUNTS)
        r2 = _make_result(df2, DataType.NORMALIZED)
        result = validate_merge_compatibility([r1, r2])
        assert any("Mixed" in i or "mixed" in i.lower() for i in result["issues"])

    def test_no_gene_overlap(self):
        df1 = pd.DataFrame({"A": [1]}, index=["S1"])
        df2 = pd.DataFrame({"B": [1]}, index=["S2"])
        result = validate_merge_compatibility([_make_result(df1), _make_result(df2)])
        assert any("No shared" in i or "no shared" in i.lower() or "overlap" in i.lower() for i in result["issues"])


class TestMergeParseResults:
    def test_single_file_passthrough(self):
        df = pd.DataFrame({"A": [1, 2]}, index=["S1", "S2"])
        result = merge_parse_results([_make_result(df)], ["file1.csv"])
        assert result.expression_df.equals(df)

    def test_two_file_merge_inner_join(self):
        df1 = pd.DataFrame({"GeneA": [100, 200], "GeneB": [300, 400]}, index=["S1", "S2"])
        df2 = pd.DataFrame({"GeneA": [500, 600], "GeneC": [700, 800]}, index=["S3", "S4"])
        result = merge_parse_results([_make_result(df1), _make_result(df2)], ["f1.csv", "f2.csv"])
        assert result.expression_df.shape == (4, 1)  # inner join: only GeneA
        assert list(result.expression_df.columns) == ["GeneA"]
        assert result.expression_df.shape[0] == 4

    def test_full_overlap_merge(self):
        df1 = pd.DataFrame({"A": [1, 2], "B": [3, 4]}, index=["S1", "S2"])
        df2 = pd.DataFrame({"A": [5, 6], "B": [7, 8]}, index=["S3", "S4"])
        result = merge_parse_results([_make_result(df1), _make_result(df2)], ["f1.csv", "f2.csv"])
        assert result.expression_df.shape == (4, 2)
        assert set(result.expression_df.index) == {"S1", "S2", "S3", "S4"}

    def test_duplicate_sample_names(self):
        df1 = pd.DataFrame({"A": [1]}, index=["Sample1"])
        df2 = pd.DataFrame({"A": [2]}, index=["Sample1"])
        result = merge_parse_results([_make_result(df1), _make_result(df2)], ["batch1.csv", "batch2.csv"])
        assert result.expression_df.shape[0] == 2
        # Sample names should be deduplicated
        assert len(set(result.expression_df.index)) == 2

    def test_mixed_type_raises(self):
        df1 = pd.DataFrame({"A": [1]}, index=["S1"])
        df2 = pd.DataFrame({"A": [2]}, index=["S2"])
        r1 = _make_result(df1, DataType.RAW_COUNTS)
        r2 = _make_result(df2, DataType.NORMALIZED)
        with pytest.raises(ValueError, match="[Dd]ifferent|[Mm]ixed"):
            merge_parse_results([r1, r2], ["f1.csv", "f2.csv"])

    def test_no_overlap_raises(self):
        df1 = pd.DataFrame({"A": [1]}, index=["S1"])
        df2 = pd.DataFrame({"B": [2]}, index=["S2"])
        with pytest.raises(ValueError, match="[Nn]o shared|[Cc]annot"):
            merge_parse_results([_make_result(df1), _make_result(df2)], ["f1.csv", "f2.csv"])

    def test_preserves_data_type(self):
        df1 = pd.DataFrame({"A": [1]}, index=["S1"])
        df2 = pd.DataFrame({"A": [2]}, index=["S2"])
        result = merge_parse_results([_make_result(df1), _make_result(df2)], ["f1.csv", "f2.csv"])
        assert result.data_type == DataType.RAW_COUNTS
        assert result.can_run_de is True
