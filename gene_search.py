"""
Gene search utility for RNA-seq analysis platform.

Provides functions for searching and retrieving gene information from DataFrames.
"""

from typing import Optional
import pandas as pd


def search_genes(
    df: pd.DataFrame, query: str, column: str = "gene", case_sensitive: bool = False
) -> pd.DataFrame:
    """
    Search for genes in DataFrame using partial substring matching.

    Supports flexible gene lookup with case-insensitive matching and fallback
    to index if column doesn't exist.

    Args:
        df: DataFrame with gene information (genes as rows or index)
        query: Search string (partial substring matching)
        column: Column name containing gene names (default: "gene")
        case_sensitive: Whether to perform case-sensitive matching (default: False)

    Returns:
        Filtered DataFrame containing rows matching the query

    Behavior:
        - Empty query: Returns full DataFrame
        - Missing column: Tries to search in index instead
        - Case-insensitive by default: "brca" matches "BRCA1", "Brca2"
        - Partial matching: "br" matches "BRCA1", "BRCA2", "BRAF"

    Examples:
        >>> df = pd.DataFrame({'gene': ['TP53', 'BRCA1', 'MYC', 'BRCA2']})
        >>> search_genes(df, 'br')
        gene
        1  BRCA1
        3  BRCA2

        >>> search_genes(df, 'TP')
        gene
        0  TP53

        >>> search_genes(df, '')  # Empty query returns all
        gene
        0  TP53
        1  BRCA1
        2  MYC
        3  BRCA2
    """
    # Handle empty query
    if not query:
        return df.copy()

    # Try to use specified column, fallback to index
    if column in df.columns:
        search_target = df[column]
    elif column == "gene" and df.index.name == "gene":
        # Fallback: search in index if column is "gene" and index is named "gene"
        search_target = df.index
    else:
        # Try index as last resort
        search_target = df.index

    # Perform case-sensitive or case-insensitive matching
    if case_sensitive:
        mask = search_target.astype(str).str.contains(query, regex=False)
    else:
        mask = search_target.astype(str).str.contains(
            query, case=False, regex=False
        )

    return df[mask].copy()


def get_gene_summary(
    df: pd.DataFrame, gene_name: str, column: str = "gene"
) -> Optional[pd.Series]:
    """
    Retrieve summary information for a specific gene.

    Performs exact match first, then falls back to case-insensitive matching.
    Returns None if gene not found.

    Args:
        df: DataFrame with gene information (genes as rows or index)
        gene_name: Gene name to search for
        column: Column name containing gene names (default: "gene")

    Returns:
        pd.Series with gene information if found, None otherwise

    Behavior:
        - Exact match: Returns row if gene_name matches exactly
        - Case-insensitive fallback: If exact match fails, tries case-insensitive
        - Index fallback: If column doesn't exist, searches in index
        - Returns None: If gene not found after all attempts

    Examples:
        >>> df = pd.DataFrame({'gene': ['TP53', 'BRCA1', 'MYC'], 'log2fc': [1.5, -2.0, 0.8]})
        >>> get_gene_summary(df, 'TP53')
        gene       TP53
        log2fc      1.5
        Name: 0, dtype: object

        >>> get_gene_summary(df, 'tp53')  # Case-insensitive fallback
        gene       TP53
        log2fc      1.5
        Name: 0, dtype: object

        >>> get_gene_summary(df, 'NOTFOUND')
        None
    """
    # Determine search target
    if column in df.columns:
        search_target = df[column]
        use_index = False
    else:
        # Fallback to index
        search_target = df.index
        use_index = True

    # Try exact match first
    exact_matches = search_target == gene_name
    if exact_matches.any():
        if use_index:
            return df[exact_matches].iloc[0]
        else:
            return df[exact_matches].iloc[0]

    # Fallback to case-insensitive match
    case_insensitive_matches = search_target.astype(str).str.lower() == gene_name.lower()
    if case_insensitive_matches.any():
        if use_index:
            return df[case_insensitive_matches].iloc[0]
        else:
            return df[case_insensitive_matches].iloc[0]

    # Not found
    return None
