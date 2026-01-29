"""
Pathway Enrichment Analysis Module

Provides pathway enrichment analysis using GSEApy Enrichr API.
Queries GO Biological Process and KEGG pathway databases.

Classes:
    PathwayEnrichment: Main class for pathway enrichment analysis

Author: RNA-seq Analysis Platform
"""

import gseapy as gp
import pandas as pd
from typing import List, Tuple, Optional
from de_analysis import ensure_gene_column
from de_analysis import ensure_gene_column


class PathwayEnrichment:
    """
    Pathway enrichment analysis using GSEApy Enrichr API.

    Supports:
    - Gene selection from DE results with adaptive thresholds
    - GO Biological Process enrichment
    - KEGG pathway enrichment
    - Graceful offline handling
    """

    def select_genes_for_enrichment(
        self,
        de_results: pd.DataFrame,
        padj_threshold: float = 0.05,
        lfc_threshold: float = 1.0,
    ) -> Tuple[List[str], Optional[str]]:
        """
        Select genes for pathway enrichment analysis.

        Rules:
        1. Filter by padj < padj_threshold (default 0.05)
        2. Filter by |log2FoldChange| > lfc_threshold (default 1.0)
        3. Combine UP and DOWN regulated genes into single list
        4. If < 5 genes pass filters: return all genes with padj < 0.1 (relaxed)
        5. If still < 5 genes: return empty list + warning message
        6. Cap at 500 genes maximum (Enrichr recommendation)

        Args:
            de_results: DataFrame with columns: gene, log2FoldChange, padj
            padj_threshold: Adjusted p-value threshold (default 0.05)
            lfc_threshold: Absolute log2 fold change threshold (default 1.0)

        Returns:
            Tuple of (gene_list, error_message)
            - gene_list: List of gene symbols sorted by padj ascending
            - error_message: None if successful, error string if too few genes
        """
        # Ensure gene column is properly named (handles various naming conventions)
        de_results = ensure_gene_column(de_results)

        # Drop NaN values before filtering (normal for low-count genes)
        de_results = de_results.dropna(subset=["padj", "log2FoldChange"])

        # Filter by padj and |log2FoldChange|
        sig = de_results[
            (de_results["padj"] < padj_threshold)
            & (de_results["log2FoldChange"].abs() > lfc_threshold)
        ]

        # If < 5 genes, relax threshold to padj < 0.1
        if len(sig) < 5:
            sig = de_results[de_results["padj"] < 0.1]

        # If still < 5 genes, return empty list with warning
        if len(sig) < 5:
            return [], "Too few significant genes for enrichment analysis"

        # Sort by padj ascending, take top 500
        sig = sig.sort_values("padj").head(500)

        return sig["gene"].tolist(), None

    def run_enrichment(
        self, gene_list: List[str], gene_sets: List[str], organism: str = "Human"
    ) -> Tuple[pd.DataFrame, Optional[str]]:
        """
        Run Enrichr API query for pathway enrichment.

        Args:
            gene_list: List of gene symbols
            gene_sets: List of gene set libraries (e.g., ['GO_Biological_Process_2023'])
            organism: Organism name (default 'Human')

        Returns:
            Tuple of (results_df, error_message)
            - results_df: DataFrame with columns: Term, Overlap, P-value, Adjusted P-value, Genes
            - error_message: None if successful, error string if API fails
        """
        # Handle empty gene list
        if not gene_list:
            return pd.DataFrame(), "No genes provided for enrichment"

        try:
            # Query Enrichr API
            enr = gp.enrichr(
                gene_list=gene_list,
                gene_sets=gene_sets,
                organism=organism,
                outdir=None,  # Don't save to disk
                cutoff=0.05,
            )

            # Format results
            results_df = self.format_results(enr)

            return results_df, None

        except Exception as e:
            # Graceful offline handling
            error_msg = f"Enrichment analysis failed (possibly offline): {str(e)}"
            return pd.DataFrame(), error_msg

    def get_go_enrichment(self, genes: List[str]) -> Tuple[pd.DataFrame, Optional[str]]:
        """
        Query GO Biological Process 2023 database.

        Args:
            genes: List of gene symbols

        Returns:
            Tuple of (results_df, error_message)
            - results_df: Top 20 GO terms with columns: Term, Overlap, P-value, Adjusted P-value, Genes
            - error_message: None if successful, error string if API fails
        """
        return self.run_enrichment(
            gene_list=genes, gene_sets=["GO_Biological_Process_2023"], organism="Human"
        )

    def get_kegg_enrichment(
        self, genes: List[str]
    ) -> Tuple[pd.DataFrame, Optional[str]]:
        """
        Query KEGG 2021 Human pathway database.

        Args:
            genes: List of gene symbols

        Returns:
            Tuple of (results_df, error_message)
            - results_df: Top 20 KEGG pathways with columns: Term, Overlap, P-value, Adjusted P-value, Genes
            - error_message: None if successful, error string if API fails
        """
        return self.run_enrichment(
            gene_list=genes, gene_sets=["KEGG_2021_Human"], organism="Human"
        )

    def format_results(self, enr_results) -> pd.DataFrame:
        """
        Standardize enrichr results to consistent format.

        Args:
            enr_results: GSEApy enrichr result object

        Returns:
            DataFrame with columns: Term, Overlap, P-value, Adjusted P-value, Genes
            Limited to top 20 pathways sorted by Adjusted P-value
        """
        # Extract results DataFrame from GSEApy object
        results_df = enr_results.results

        # Select and rename columns to standardized format
        standardized = results_df[
            ["Term", "Overlap", "P-value", "Adjusted P-value", "Genes"]
        ].copy()

        # Sort by Adjusted P-value ascending
        standardized = standardized.sort_values("Adjusted P-value")

        # Return top 20 pathways
        return standardized.head(20)
