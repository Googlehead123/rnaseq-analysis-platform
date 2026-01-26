"""
Excel export module for RNA-seq analysis results.

Exports multi-sheet Excel workbooks with DE results, enrichment analysis,
and analysis metadata. Sheet layout adapts to data type (RAW_COUNTS, NORMALIZED, PRE_ANALYZED).
"""

from dataclasses import dataclass
from typing import Dict, Tuple, Optional, List, Any
from datetime import datetime
import re
import sys
import pandas as pd
import plotly.graph_objects as go
from rnaseq_parser import DataType
from de_analysis import DEResult


@dataclass
class EnrichmentResult:
    """Pathway enrichment results for a single comparison."""

    go_results: pd.DataFrame  # Columns: Term, Overlap, P-value, Adjusted P-value, Genes
    kegg_results: pd.DataFrame  # Same columns
    genes_used: List[str]  # Genes submitted to Enrichr
    selection_note: str  # E.g., "172 genes (padj<0.05, |log2FC|>1)"
    error: Optional[str] = None  # None if successful, error message if failed


@dataclass
class ExportData:
    """Complete export data bundle - constructed by UI before calling export."""

    # Data type determines which fields are populated
    data_type: DataType

    # DE Results (RAW_COUNTS + PRE_ANALYZED)
    # For RAW_COUNTS: contains all comparisons from st.session_state['de_results']
    # For PRE_ANALYZED: single-entry dict with key ("preanalyzed", "preanalyzed")
    de_results: Dict[Tuple[str, str], DEResult]  # Empty dict for NORMALIZED

    # Expression matrix (RAW_COUNTS + NORMALIZED)
    # Used for heatmap/PCA in PDF
    expression_matrix: Optional[pd.DataFrame]  # samples × genes, log2(norm+1)

    # Enrichment results (per comparison for RAW_COUNTS, single for PRE_ANALYZED)
    enrichment_results: Dict[Tuple[str, str], EnrichmentResult]  # Empty dict if none

    # Figures (Plotly Figure objects - export converts to PNG as needed)
    figures: Dict[str, go.Figure]  # Keys: "volcano", "heatmap", "pca", "panel_{name}"

    # Settings/metadata
    settings: Dict[str, Any]  # padj_threshold, lfc_threshold, comparisons list, etc.
    sample_conditions: Dict[str, str]  # sample → condition mapping (for reports)

    # Active comparison (for title/context in single-comparison views)
    active_comparison: Optional[Tuple[str, str]] = None


class ExportEngine:
    """Excel export engine for RNA-seq analysis results."""

    def sanitize_sheet_name(self, name: str, max_length: int = 31) -> str:
        """
        Sanitize sheet name for Excel compatibility.
        
        Excel sheet name rules:
        - Max 31 characters
        - Cannot contain: [ ] : * ? / \
        - Cannot start or end with '
        
        Args:
            name: Raw sheet name
            max_length: Maximum length (default 31 for Excel)
            
        Returns:
            Sanitized sheet name
        """
        # Remove forbidden characters
        name = re.sub(r"[\[\]:*?/\\]", "_", name)
        # Remove leading/trailing quotes
        name = name.strip("'")
        # Truncate
        return name[:max_length]

    def export_excel(self, filepath: str, export_data: ExportData) -> None:
        """
        Export analysis results to multi-sheet Excel workbook.

        Sheet layout depends on data_type:
        - RAW_COUNTS: DE_{comparison}, Sig_{comparison}, GO_{comparison}, KEGG_{comparison}, Settings
        - NORMALIZED: Expression Matrix, Settings
        - PRE_ANALYZED: DE Results, Significant Genes, GO Enrichment, KEGG Enrichment, Settings

        Args:
            filepath: Output Excel file path (.xlsx)
            export_data: Complete export data bundle
        """
        with pd.ExcelWriter(filepath, engine="openpyxl") as writer:
            if export_data.data_type == DataType.RAW_COUNTS:
                # Write DE results per comparison
                for (test, ref), de_result in export_data.de_results.items():
                    # DE sheet
                    sheet_name = self.sanitize_sheet_name(f"DE_{test}_vs_{ref}")
                    de_result.results_df.to_excel(
                        writer, sheet_name=sheet_name, index=False
                    )

                    # Significant genes sheet
                    sig_sheet = self.sanitize_sheet_name(f"Sig_{test}_vs_{ref}")
                    sig_genes = de_result.results_df[
                        de_result.results_df["padj"] < 0.05
                    ]
                    sig_genes.to_excel(writer, sheet_name=sig_sheet, index=False)

                # Write enrichment results per comparison
                for (
                    test,
                    ref,
                ), enrich_result in export_data.enrichment_results.items():
                    if enrich_result.error is None:
                        go_sheet = self.sanitize_sheet_name(f"GO_{test}_vs_{ref}")
                        enrich_result.go_results.to_excel(
                            writer, sheet_name=go_sheet, index=False
                        )

                        kegg_sheet = self.sanitize_sheet_name(f"KEGG_{test}_vs_{ref}")
                        enrich_result.kegg_results.to_excel(
                            writer, sheet_name=kegg_sheet, index=False
                        )

            elif export_data.data_type == DataType.NORMALIZED:
                # Expression matrix only
                if export_data.expression_matrix is not None:
                    export_data.expression_matrix.to_excel(
                        writer, sheet_name="Expression Matrix"
                    )

            elif export_data.data_type == DataType.PRE_ANALYZED:
                # Single DE result
                if export_data.de_results:
                    de_result = list(export_data.de_results.values())[0]
                    de_result.results_df.to_excel(
                        writer, sheet_name="DE Results", index=False
                    )

                    sig_genes = de_result.results_df[
                        de_result.results_df["padj"] < 0.05
                    ]
                    sig_genes.to_excel(
                        writer, sheet_name="Significant Genes", index=False
                    )

                # Enrichment if available
                if export_data.enrichment_results:
                    enrich_result = list(export_data.enrichment_results.values())[0]
                    if enrich_result.error is None:
                        enrich_result.go_results.to_excel(
                            writer, sheet_name="GO Enrichment", index=False
                        )
                        enrich_result.kegg_results.to_excel(
                            writer, sheet_name="KEGG Enrichment", index=False
                        )

            # Settings sheet (all modes)
            self._write_settings_sheet(writer, export_data)

    def _write_settings_sheet(
        self, writer: pd.ExcelWriter, export_data: ExportData
    ) -> None:
        """
        Write Settings sheet with analysis metadata.

        Settings sheet contains key-value rows with sections:
        - Analysis Date, Data Type, Python Version, PyDESeq2 Version
        - Thresholds (padj, log2FC)
        - Comparisons (with success/failure status)
        - Enrichment Status (per comparison)
        - Sample Conditions (sample → condition mapping)

        Args:
            writer: pd.ExcelWriter instance
            export_data: Complete export data bundle
        """
        settings_data = [
            ["Parameter", "Value"],  # Header row
            ["Analysis Date", datetime.now().strftime("%Y-%m-%d %H:%M:%S")],
            ["Data Type", export_data.data_type.value],
            [
                "Python Version",
                f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}",
            ],
        ]

        # Add PyDESeq2 version if available
        try:
            import pydeseq2

            settings_data.append(["PyDESeq2 Version", pydeseq2.__version__])
        except (ImportError, AttributeError):
            settings_data.append(["PyDESeq2 Version", "N/A"])

        # Thresholds section (if available in settings)
        if export_data.settings:
            settings_data.append(["---", "---"])  # Separator
            settings_data.append(["Thresholds", ""])

            if "padj_threshold" in export_data.settings:
                settings_data.append(
                    ["padj Threshold", str(export_data.settings["padj_threshold"])]
                )
            if "lfc_threshold" in export_data.settings:
                settings_data.append(
                    ["log2FC Threshold", str(export_data.settings["lfc_threshold"])]
                )

        # Comparisons section (RAW_COUNTS and PRE_ANALYZED)
        if export_data.de_results:
            settings_data.append(["---", "---"])
            settings_data.append(["Comparisons", ""])

            for (test, ref), de_result in export_data.de_results.items():
                comparison_name = f"{test}_vs_{ref}"
                if de_result.warnings:
                    # Failed comparison
                    status = f"FAILED ({de_result.warnings[0]})"
                else:
                    # Successful comparison
                    status = f"SUCCESS ({de_result.n_significant} significant genes)"
                settings_data.append([comparison_name, status])

        # Enrichment status section
        if export_data.enrichment_results:
            settings_data.append(["---", "---"])
            settings_data.append(["Enrichment Status", ""])

            for (test, ref), enrich_result in export_data.enrichment_results.items():
                comparison_name = f"{test}_vs_{ref}"

                # GO enrichment status
                if enrich_result.error:
                    go_status = f"FAILED ({enrich_result.error})"
                else:
                    go_count = len(enrich_result.go_results)
                    go_status = f"SUCCESS ({go_count} pathways)"
                settings_data.append([f"GO_{comparison_name}", go_status])

                # KEGG enrichment status
                if enrich_result.error:
                    kegg_status = f"FAILED ({enrich_result.error})"
                else:
                    kegg_count = len(enrich_result.kegg_results)
                    kegg_status = f"SUCCESS ({kegg_count} pathways)"
                settings_data.append([f"KEGG_{comparison_name}", kegg_status])

        # Sample conditions section (RAW_COUNTS only)
        if (
            export_data.data_type == DataType.RAW_COUNTS
            and export_data.sample_conditions
        ):
            settings_data.append(["---", "---"])
            settings_data.append(["Sample Conditions", ""])

            for sample, condition in sorted(export_data.sample_conditions.items()):
                settings_data.append([sample, condition])

        # Write settings sheet
        settings_df = pd.DataFrame(settings_data)
        settings_df.to_excel(writer, sheet_name="Settings", index=False, header=False)
