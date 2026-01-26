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
import io
import pandas as pd
import plotly.graph_objects as go
from reportlab.lib.pagesizes import letter
from reportlab.platypus import SimpleDocTemplate, Paragraph, Table, Image, Spacer
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.lib.utils import ImageReader
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

    def export_figure(
        self, fig: go.Figure, filepath: str, format: str = "png", scale: int = 3
    ) -> None:
        """
        Export Plotly figure to static image file.

        Args:
            fig: Plotly Figure object
            filepath: Output file path
            format: Image format ('png', 'svg', 'pdf')
            scale: Scale factor for raster formats (default 3 for ~300 DPI)
                   - scale=3 with default 700x500 = 2100x1500 pixels ≈ 300 DPI at 7x5 inches
                   - For SVG/PDF (vector), scale doesn't matter

        Example:
            engine = ExportEngine()
            engine.export_figure(volcano_fig, "volcano.png", format="png", scale=3)
            engine.export_figure(heatmap_fig, "heatmap.svg", format="svg")
        """
        fig.write_image(filepath, format=format, scale=scale)

    def export_pdf_report(self, filepath: str, export_data: ExportData) -> None:
        """
        Generate comprehensive PDF report from ExportData bundle.

        PDF layout adapts to data_type:
        - RAW_COUNTS: Full report with DE, volcano, enrichment, heatmap, PCA, panels
        - NORMALIZED: Heatmap, PCA, panels only (no DE/enrichment)
        - PRE_ANALYZED: DE, volcano, enrichment only (no heatmap/PCA/panels)

        Figures embedded as in-memory PNG bytes (no temp files).

        Args:
            filepath: Output PDF file path
            export_data: Complete export data bundle
        """
        doc = SimpleDocTemplate(filepath, pagesize=letter)
        styles = getSampleStyleSheet()
        story = []

        # 1. Title page
        story.append(Paragraph("RNA-seq Analysis Report", styles["Title"]))
        story.append(Spacer(1, 12))
        story.append(
            Paragraph(
                f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
                styles["Normal"],
            )
        )
        story.append(Spacer(1, 24))

        # 2. Methods section
        story.append(Paragraph("Methods", styles["Heading1"]))
        story.append(
            Paragraph(f"Data Type: {export_data.data_type.value}", styles["Normal"])
        )
        story.append(Spacer(1, 6))

        if export_data.data_type == DataType.RAW_COUNTS:
            # Add thresholds for RAW_COUNTS mode
            story.append(
                Paragraph(
                    f"Thresholds: padj &lt; {export_data.settings.get('padj_threshold', 0.05)}, "
                    f"|log2FC| &gt; {export_data.settings.get('lfc_threshold', 1.0)}",
                    styles["Normal"],
                )
            )
            story.append(Spacer(1, 6))

            # List comparisons
            comparisons = export_data.settings.get("comparisons", [])
            if comparisons:
                story.append(
                    Paragraph(f"Comparisons: {len(comparisons)}", styles["Normal"])
                )
                for test, ref in comparisons:
                    story.append(Paragraph(f"  • {test} vs {ref}", styles["Normal"]))
        elif export_data.data_type == DataType.NORMALIZED:
            story.append(
                Paragraph(
                    "Normalized expression data provided for visualization only.",
                    styles["Normal"],
                )
            )
        elif export_data.data_type == DataType.PRE_ANALYZED:
            story.append(
                Paragraph(
                    "Pre-analyzed differential expression results uploaded.",
                    styles["Normal"],
                )
            )

        story.append(Spacer(1, 24))

        # 3. DE Results table (RAW_COUNTS and PRE_ANALYZED only)
        if export_data.data_type in [DataType.RAW_COUNTS, DataType.PRE_ANALYZED]:
            if export_data.de_results:
                story.append(
                    Paragraph("Top Differentially Expressed Genes", styles["Heading1"])
                )
                story.append(Spacer(1, 6))

                # Get active comparison or first available
                active = (
                    export_data.active_comparison
                    or list(export_data.de_results.keys())[0]
                )
                de_result = export_data.de_results[active]

                # Add comparison name
                if export_data.data_type == DataType.RAW_COUNTS:
                    test, ref = active
                    story.append(
                        Paragraph(f"Comparison: {test} vs {ref}", styles["Normal"])
                    )
                    story.append(Spacer(1, 6))

                # Top 20 genes by padj
                top_genes = de_result.results_df.nsmallest(20, "padj")[
                    ["gene", "log2FoldChange", "padj"]
                ].values.tolist()

                table_data = [["Gene", "log2FC", "padj"]] + [
                    [str(g), f"{fc:.2f}", f"{p:.2e}"] for g, fc, p in top_genes
                ]
                story.append(Table(table_data))
                story.append(Spacer(1, 24))
        else:
            # NORMALIZED mode - omit DE section
            story.append(
                Paragraph("Differential Expression Analysis", styles["Heading1"])
            )
            story.append(Spacer(1, 6))
            story.append(
                Paragraph(
                    "Section not available for normalized data.", styles["Normal"]
                )
            )
            story.append(Spacer(1, 24))

        # 4. Volcano plot (RAW_COUNTS and PRE_ANALYZED only)
        if export_data.data_type in [DataType.RAW_COUNTS, DataType.PRE_ANALYZED]:
            if "volcano" in export_data.figures:
                story.append(Paragraph("Volcano Plot", styles["Heading1"]))
                story.append(Spacer(1, 6))
                volcano_fig = export_data.figures["volcano"]
                png_bytes = volcano_fig.to_image(
                    format="png", scale=2, width=800, height=600
                )
                img_buffer = io.BytesIO(png_bytes)
                story.append(Image(ImageReader(img_buffer), width=400, height=300))
                story.append(Spacer(1, 24))

        # 5. Pathway enrichment (RAW_COUNTS and PRE_ANALYZED only)
        if export_data.data_type in [DataType.RAW_COUNTS, DataType.PRE_ANALYZED]:
            if export_data.enrichment_results:
                story.append(Paragraph("Pathway Enrichment", styles["Heading1"]))
                story.append(Spacer(1, 6))

                active = (
                    export_data.active_comparison
                    or list(export_data.enrichment_results.keys())[0]
                )
                enrich_result = export_data.enrichment_results[active]

                if enrich_result.error is None:
                    # GO results (top 10)
                    story.append(
                        Paragraph("GO Biological Process (Top 10)", styles["Heading2"])
                    )
                    story.append(Spacer(1, 6))
                    go_top = enrich_result.go_results.head(10)[
                        ["Term", "Adjusted P-value"]
                    ].values.tolist()
                    go_table = [["Term", "Adj. P-value"]] + [
                        [term, f"{p:.2e}"] for term, p in go_top
                    ]
                    story.append(Table(go_table))
                    story.append(Spacer(1, 12))

                    # KEGG results (top 10)
                    story.append(
                        Paragraph("KEGG Pathways (Top 10)", styles["Heading2"])
                    )
                    story.append(Spacer(1, 6))
                    kegg_top = enrich_result.kegg_results.head(10)[
                        ["Term", "Adjusted P-value"]
                    ].values.tolist()
                    kegg_table = [["Term", "Adj. P-value"]] + [
                        [term, f"{p:.2e}"] for term, p in kegg_top
                    ]
                    story.append(Table(kegg_table))
                    story.append(Spacer(1, 24))
                else:
                    story.append(
                        Paragraph(
                            f"Enrichment analysis failed: {enrich_result.error}",
                            styles["Normal"],
                        )
                    )
                    story.append(Spacer(1, 24))

        # 6. Gene panels (RAW_COUNTS and NORMALIZED only - need expression matrix)
        if export_data.data_type in [DataType.RAW_COUNTS, DataType.NORMALIZED]:
            panel_figs = {
                k: v for k, v in export_data.figures.items() if k.startswith("panel_")
            }
            if panel_figs:
                story.append(Paragraph("Gene Panel Analysis", styles["Heading1"]))
                story.append(Spacer(1, 6))
                for panel_name, panel_fig in panel_figs.items():
                    story.append(
                        Paragraph(
                            panel_name.replace("panel_", "").replace("_", " ").title(),
                            styles["Heading2"],
                        )
                    )
                    story.append(Spacer(1, 6))
                    png_bytes = panel_fig.to_image(
                        format="png", scale=2, width=800, height=600
                    )
                    img_buffer = io.BytesIO(png_bytes)
                    story.append(Image(ImageReader(img_buffer), width=400, height=300))
                    story.append(Spacer(1, 12))
                story.append(Spacer(1, 12))

        # 7. Heatmap (RAW_COUNTS and NORMALIZED only)
        if export_data.data_type in [DataType.RAW_COUNTS, DataType.NORMALIZED]:
            if "heatmap" in export_data.figures:
                story.append(Paragraph("Clustered Heatmap", styles["Heading1"]))
                story.append(Spacer(1, 6))
                heatmap_fig = export_data.figures["heatmap"]
                png_bytes = heatmap_fig.to_image(
                    format="png", scale=2, width=800, height=600
                )
                img_buffer = io.BytesIO(png_bytes)
                story.append(Image(ImageReader(img_buffer), width=400, height=300))
                story.append(Spacer(1, 24))

        # 8. PCA (RAW_COUNTS and NORMALIZED only)
        if export_data.data_type in [DataType.RAW_COUNTS, DataType.NORMALIZED]:
            if "pca" in export_data.figures:
                story.append(Paragraph("PCA Plot", styles["Heading1"]))
                story.append(Spacer(1, 6))
                pca_fig = export_data.figures["pca"]
                png_bytes = pca_fig.to_image(
                    format="png", scale=2, width=800, height=600
                )
                img_buffer = io.BytesIO(png_bytes)
                story.append(Image(ImageReader(img_buffer), width=400, height=300))
                story.append(Spacer(1, 24))

        # Build PDF
        doc.build(story)
