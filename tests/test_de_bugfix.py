"""Tests for the KeyError: 'gene' bug fix in de_analysis."""
import pytest
import pandas as pd
import numpy as np
from de_analysis import ensure_gene_column


def test_ensure_gene_column_with_named_index():
    df = pd.DataFrame(
        {"log2FoldChange": [1.5, -2.0], "padj": [0.01, 0.05]},
        index=pd.Index(["TP53", "BRCA1"], name="Gene"),
    )
    result = ensure_gene_column(df)
    assert "gene" in result.columns
    assert list(result["gene"]) == ["TP53", "BRCA1"]


def test_ensure_gene_column_with_uppercase_column():
    df = pd.DataFrame({"Gene": ["TP53", "BRCA1"], "padj": [0.01, 0.05]})
    result = ensure_gene_column(df)
    assert "gene" in result.columns


def test_ensure_gene_column_already_correct():
    df = pd.DataFrame({"gene": ["TP53"], "padj": [0.01]})
    result = ensure_gene_column(df)
    assert "gene" in result.columns
    assert len(result) == 1


def test_ensure_gene_column_symbol_alias():
    df = pd.DataFrame({"SYMBOL": ["TP53", "MYC"], "padj": [0.01, 0.05]})
    result = ensure_gene_column(df)
    assert "gene" in result.columns


def test_ensure_gene_column_idempotent():
    df = pd.DataFrame(
        {"log2FoldChange": [1.5], "padj": [0.01]},
        index=pd.Index(["TP53"], name="Gene"),
    )
    result1 = ensure_gene_column(df)
    result2 = ensure_gene_column(result1)
    assert "gene" in result2.columns
    assert list(result2["gene"]) == ["TP53"]
