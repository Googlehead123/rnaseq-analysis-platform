"""Tests for demo dataset module."""
import pytest
import pandas as pd
from demo_data import load_demo_dataset, get_demo_description


def test_load_demo_dataset_shapes():
    counts_df, metadata_df = load_demo_dataset()
    assert counts_df.shape[0] == 6
    assert counts_df.shape[1] > 50
    assert metadata_df.shape[0] == 6


def test_load_demo_dataset_samples():
    counts_df, _ = load_demo_dataset()
    assert all("Control" in s or "Treatment" in s for s in counts_df.index)


def test_load_demo_dataset_conditions():
    _, metadata_df = load_demo_dataset()
    assert "condition" in metadata_df.columns
    assert set(metadata_df["condition"].unique()) == {"Control", "Treatment"}


def test_load_demo_dataset_integer_counts():
    counts_df, _ = load_demo_dataset()
    for col in counts_df.columns:
        assert pd.api.types.is_integer_dtype(counts_df[col])


def test_load_demo_dataset_reproducible():
    df1, _ = load_demo_dataset()
    df2, _ = load_demo_dataset()
    assert df1.equals(df2)


def test_load_demo_dataset_index_match():
    counts_df, metadata_df = load_demo_dataset()
    assert set(counts_df.index) == set(metadata_df.index)


def test_get_demo_description():
    desc = get_demo_description()
    assert isinstance(desc, str)
    assert len(desc) > 50
