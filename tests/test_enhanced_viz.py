"""Tests for enhanced visualization functions."""
import pytest
import pandas as pd
import numpy as np
from visualizations import create_volcano_plot, create_ma_plot, compute_de_summary


@pytest.fixture
def de_results_df():
    np.random.seed(42)
    n = 100
    return pd.DataFrame({
        "gene": [f"Gene_{i}" for i in range(n)],
        "log2FoldChange": np.random.randn(n) * 2,
        "padj": np.random.uniform(0, 1, n),
        "baseMean": np.random.uniform(10, 10000, n),
        "pvalue": np.random.uniform(0, 1, n),
    })


def test_volcano_plot_creates_figure(de_results_df):
    fig = create_volcano_plot(de_results_df)
    assert fig is not None
    assert hasattr(fig, "data")


def test_ma_plot_creates_figure(de_results_df):
    fig = create_ma_plot(de_results_df)
    assert fig is not None
    assert hasattr(fig, "data")


def test_compute_de_summary_keys(de_results_df):
    summary = compute_de_summary(de_results_df)
    assert summary["total_genes"] == 100
    assert "significant_genes" in summary
    assert "upregulated" in summary
    assert "downregulated" in summary
    assert "top_up_genes" in summary
    assert "top_down_genes" in summary
    assert "top_significant" in summary


def test_compute_de_summary_counts_consistent(de_results_df):
    summary = compute_de_summary(de_results_df)
    assert summary["upregulated"] + summary["downregulated"] <= summary["significant_genes"]
    assert summary["significant_genes"] <= summary["total_genes"]
