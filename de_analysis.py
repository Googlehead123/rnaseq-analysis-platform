"""
Differential expression analysis using PyDESeq2.

Implements "fit once, contrast many" model for efficient multi-comparison analysis.
"""

from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional
import pandas as pd
import numpy as np
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats


@dataclass
class DEResult:
    """Result from differential expression analysis."""

    results_df: pd.DataFrame  # DE results with columns: gene, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj
    normalized_counts: Optional[
        pd.DataFrame
    ]  # samples × genes (raw normalized) - None for PRE_ANALYZED
    log_normalized_counts: Optional[
        pd.DataFrame
    ]  # samples × genes (log2(norm+1)) - None for PRE_ANALYZED
    dds: Optional[DeseqDataSet]  # Fitted model - None for PRE_ANALYZED
    comparison: Tuple[str, str]  # (test_condition, reference_condition)
    n_significant: int  # Count of genes with padj < 0.05
    warnings: List[str]  # Any warnings during analysis


class DEAnalysisEngine:
    """Differential expression analysis using PyDESeq2."""

    def fit_model(
        self,
        counts_df: pd.DataFrame,
        metadata_df: pd.DataFrame,
        design_factor: str = "condition",
    ) -> Tuple[DeseqDataSet, pd.DataFrame, pd.DataFrame]:
        """
        Fit DESeq2 model ONCE. Returns fitted model + normalized counts.

        Args:
            counts_df: samples × genes DataFrame with integer counts
            metadata_df: samples × conditions DataFrame (index must match counts_df.index)
            design_factor: Column name in metadata_df (default: "condition")

        Returns:
            dds: Fitted DeseqDataSet (reuse for multiple contrasts)
            normalized_df: samples × genes (raw normalized)
            log_normalized_df: samples × genes (log2(norm+1))

        Raises:
            Exception if fit fails (caller should abort all comparisons)
        """
        # Initialize DeseqDataSet
        dds = DeseqDataSet(
            counts=counts_df,  # samples × genes, integers
            metadata=metadata_df,  # samples × conditions
            design_factors=design_factor,  # Column name in metadata
            refit_cooks=True,
        )

        # Fit model (size factors, dispersion, GLM)
        dds.deseq2()

        # Extract normalized counts
        normalized_counts_array = dds.layers[
            "normed_counts"
        ]  # numpy array, samples × genes
        normalized_df = pd.DataFrame(
            normalized_counts_array,
            index=dds.obs_names,  # sample names
            columns=dds.var_names,  # gene names
        )

        # Log-transform for visualization
        log_normalized_df = np.log2(normalized_df + 1)

        return dds, normalized_df, log_normalized_df

    def get_comparison(
        self,
        dds: DeseqDataSet,
        test_condition: str,
        reference_condition: str,
        normalized_df: pd.DataFrame,
        log_normalized_df: pd.DataFrame,
    ) -> DEResult:
        """
        Compute single contrast from fitted model.

        Args:
            dds: Fitted DeseqDataSet
            test_condition: Test condition name
            reference_condition: Reference condition name
            normalized_df: Raw normalized counts (samples × genes)
            log_normalized_df: Log-normalized counts (samples × genes)

        Returns:
            DEResult for this comparison

        Raises:
            Exception on failure (caller should catch and mark comparison as failed)
        """
        # Get design factor from dds
        design_factor = (
            dds.design_factors[0]
            if isinstance(dds.design_factors, list)
            else dds.design_factors
        )

        # Compute statistics for this contrast
        stat_res = DeseqStats(
            dds, contrast=[design_factor, test_condition, reference_condition]
        )
        stat_res.summary()

        # Extract results DataFrame
        results_df = stat_res.results_df.reset_index().rename(columns={"index": "gene"})
        # Columns: gene, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj

        # Count significant genes
        n_sig = int((results_df["padj"] < 0.05).sum())

        return DEResult(
            results_df=results_df,
            normalized_counts=normalized_df,
            log_normalized_counts=log_normalized_df,
            dds=dds,
            comparison=(test_condition, reference_condition),
            n_significant=n_sig,
            warnings=[],
        )

    def run_all_comparisons(
        self,
        counts_df: pd.DataFrame,
        metadata_df: pd.DataFrame,
        comparisons: List[Tuple[str, str]],
        design_factor: str = "condition",
    ) -> Dict[Tuple[str, str], DEResult]:
        """
        Main entry point: fit model once, compute all contrasts.

        This is the method called by the UI. Implements "fit once, contrast many" model.

        Args:
            counts_df: samples × genes DataFrame with integer counts
            metadata_df: samples × conditions DataFrame
            comparisons: List of (test, reference) tuples
            design_factor: Column name in metadata_df (default: "condition")

        Returns:
            Dict mapping (test, ref) → DEResult
            Failed comparisons have empty results_df and warnings populated
        """
        results = {}

        try:
            # Fit model once
            dds, normalized_df, log_normalized_df = self.fit_model(
                counts_df, metadata_df, design_factor
            )
        except Exception as e:
            # Model fit failed - all comparisons fail
            for comparison in comparisons:
                results[comparison] = DEResult(
                    results_df=pd.DataFrame(),  # Empty
                    normalized_counts=None,
                    log_normalized_counts=None,
                    dds=None,
                    comparison=comparison,
                    n_significant=0,
                    warnings=[f"Model fit failed: {str(e)}"],
                )
            return results

        # Compute each comparison
        for test_cond, ref_cond in comparisons:
            try:
                result = self.get_comparison(
                    dds, test_cond, ref_cond, normalized_df, log_normalized_df
                )
                results[(test_cond, ref_cond)] = result
            except Exception as e:
                # This comparison failed, but others may succeed
                results[(test_cond, ref_cond)] = DEResult(
                    results_df=pd.DataFrame(),  # Empty
                    normalized_counts=normalized_df,  # Still provide normalized counts
                    log_normalized_counts=log_normalized_df,
                    dds=dds,
                    comparison=(test_cond, ref_cond),
                    n_significant=0,
                    warnings=[f"Comparison failed: {str(e)}"],
                )

        return results

    @staticmethod
    def filter_results(
        results_df: pd.DataFrame,
        padj_threshold: float = 0.05,
        lfc_threshold: float = 1.0,
    ) -> pd.DataFrame:
        """
        Filter DE results to significant genes.

        Args:
            results_df: DE results DataFrame
            padj_threshold: Adjusted p-value threshold (default: 0.05)
            lfc_threshold: Absolute log2 fold change threshold (default: 1.0)

        Returns:
            Filtered DataFrame with significant genes only
        """
        return results_df[
            (results_df["padj"] < padj_threshold)
            & (abs(results_df["log2FoldChange"]) > lfc_threshold)
        ].copy()
