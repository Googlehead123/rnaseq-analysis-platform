"""Tests for gene search utility."""
import pytest
import pandas as pd
from gene_search import search_genes, get_gene_summary


@pytest.fixture
def sample_df():
    return pd.DataFrame({
        "gene": ["TP53", "BRCA1", "MYC", "BRCA2", "COL1A1"],
        "log2FoldChange": [2.5, -1.8, 3.1, -0.5, 1.2],
        "padj": [0.001, 0.02, 0.0001, 0.5, 0.03],
    })


def test_search_genes_partial_match(sample_df):
    result = search_genes(sample_df, "BR")
    assert len(result) == 2
    assert "BRCA1" in result["gene"].values
    assert "BRCA2" in result["gene"].values


def test_search_genes_case_insensitive(sample_df):
    result = search_genes(sample_df, "brca")
    assert len(result) == 2


def test_search_genes_empty_query(sample_df):
    result = search_genes(sample_df, "")
    assert len(result) == len(sample_df)


def test_search_genes_no_match(sample_df):
    result = search_genes(sample_df, "ZZZZZ")
    assert len(result) == 0


def test_get_gene_summary_found(sample_df):
    result = get_gene_summary(sample_df, "TP53")
    assert result is not None
    assert result["gene"] == "TP53"


def test_get_gene_summary_not_found(sample_df):
    result = get_gene_summary(sample_df, "NONEXISTENT")
    assert result is None


def test_get_gene_summary_case_insensitive(sample_df):
    result = get_gene_summary(sample_df, "tp53")
    assert result is not None
