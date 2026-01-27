"""
Gene panel analysis for dermatology RNA-seq research.

Provides curated gene panels and scoring methods for pathway-level analysis.
"""

from typing import Dict, List, Optional
from pathlib import Path
import warnings
import yaml
import pandas as pd
import numpy as np
import plotly.graph_objects as go


class GenePanelAnalyzer:
    """
    Analyze gene expression using curated dermatology gene panels.

    Panels include:
    - Anti-aging: Collagen synthesis, ECM remodeling
    - Skin Barrier: Stratum corneum, tight junctions
    - Anti-inflammation: NF-kB pathway cytokines
    - Whitening/Melanogenesis: Melanin synthesis pathway
    - Sebum Regulation: Lipogenesis transcription factors
    - Wound Healing: Growth factors, keratins, ECM remodeling
    """

    def __init__(self, config_path: str = "config/gene_panels.yaml"):
        """
        Initialize GenePanelAnalyzer with gene panel configuration.

        Args:
            config_path: Path to YAML config file with gene panels
        """
        self.config_path = config_path
        self.panels = self.load_panels(config_path)

    def load_panels(self, config_path: str) -> Dict[str, List[str]]:
        """
        Load gene panels from YAML configuration file.

        Args:
            config_path: Path to YAML config file

        Returns:
            Dict mapping panel_name → List[gene_symbols]

        Example:
            {
                "Anti-aging": ["COL1A1", "COL3A1", "ELN", ...],
                "Skin Barrier": ["FLG", "LOR", "IVL", ...]
            }
        """
        config_file = Path(config_path)
        if not config_file.exists():
            raise FileNotFoundError(f"Gene panel config not found: {config_path}")

        with open(config_file, "r") as f:
            config = yaml.safe_load(f)

        # Extract gene lists from nested structure
        panels = {}
        for panel_name, panel_info in config["panels"].items():
            panels[panel_name] = panel_info["genes"]

        return panels

    def score_panel(
        self,
        expression_df: pd.DataFrame,
        panel_name: str,
        sample_conditions: Dict[str, str],
    ) -> Dict[str, float]:
        """
        Calculate panel score per condition using z-score normalization.

        Algorithm (from plan lines 2566-2612):
        1. Get genes in panel that exist in expression_df
        2. For each gene: compute z-score ACROSS ALL SAMPLES
           z_gene[sample] = (expr[sample] - mean(expr)) / std(expr)
        3. For each sample: compute mean z-score across panel genes
           panel_score[sample] = mean(z_gene for gene in panel)
        4. For each condition: compute mean panel score across samples
           condition_score[cond] = mean(panel_score for sample in condition)

        Args:
            expression_df: samples × genes (log2 normalized counts)
            panel_name: Name of panel from config (e.g., "Anti-aging")
            sample_conditions: Dict[sample_name → condition_name]

        Returns:
            Dict[condition_name → float score]

        Example:
            {"Control": -0.5, "Treatment": 1.2}  # Treatment upregulates panel

        Raises:
            ValueError: If expression_df is None/empty or panel not found
        """
        # Validate expression_df
        if expression_df is None or expression_df.empty:
            raise ValueError("expression_df cannot be None or empty")

        if panel_name not in self.panels:
            available = ", ".join(self.panels.keys())
            raise ValueError(f"Panel '{panel_name}' not found. Available: {available}")

        # Step 1: Get genes in panel that exist in expression_df
        panel_genes = self.panels[panel_name]
        available_genes = [g for g in panel_genes if g in expression_df.columns]

        # Warn about missing genes
        missing_genes = [g for g in panel_genes if g not in expression_df.columns]
        if missing_genes:
            warnings.warn(
                f"Panel '{panel_name}': {len(missing_genes)}/{len(panel_genes)} genes missing from expression data: {missing_genes}"
            )

        # Require minimum 2 genes for meaningful analysis
        if len(available_genes) < 2:
            raise ValueError(
                f"Panel '{panel_name}' has < 2 available genes ({len(available_genes)}/{len(panel_genes)}). "
                f"Cannot compute panel score."
            )

        # Step 2: Z-score per gene (across all samples)
        panel_data = expression_df[available_genes]
        z_scores = (panel_data - panel_data.mean()) / panel_data.std()

        # Step 3: Mean z-score per sample (across genes)
        sample_scores = z_scores.mean(axis=1)  # Series: sample → score

        # Step 4: Mean per condition
        condition_scores = {}
        for condition in set(sample_conditions.values()):
            samples_in_cond = [
                s for s, c in sample_conditions.items() if c == condition
            ]
            condition_scores[condition] = sample_scores.loc[samples_in_cond].mean()

        return condition_scores

    def plot_panel(
        self,
        expression_df: pd.DataFrame,
        panel_name: str,
        sample_conditions: Dict[str, str],
    ) -> go.Figure:
        """
        Create bar plot showing gene expression across conditions.

        Args:
            expression_df: samples × genes (log2 normalized counts)
            panel_name: Name of panel from config
            sample_conditions: Dict[sample_name → condition_name]

        Returns:
            Plotly Figure with grouped bars (genes × conditions)

        Raises:
            ValueError: If expression_df is None/empty or panel not found
        """
        # Validate expression_df
        if expression_df is None or expression_df.empty:
            raise ValueError("expression_df cannot be None or empty")

        if panel_name not in self.panels:
            available = ", ".join(self.panels.keys())
            raise ValueError(f"Panel '{panel_name}' not found. Available: {available}")

        # Get genes in panel that exist in expression_df
        panel_genes = self.panels[panel_name]
        available_genes = [g for g in panel_genes if g in expression_df.columns]

        if len(available_genes) == 0:
            raise ValueError(
                f"Panel '{panel_name}' has no genes available in expression data"
            )

        # Compute mean expression per gene per condition
        panel_data = expression_df[available_genes]

        # Create DataFrame for plotting: genes × conditions
        plot_data = []
        for condition in sorted(set(sample_conditions.values())):
            samples_in_cond = [
                s for s, c in sample_conditions.items() if c == condition
            ]
            mean_expr = panel_data.loc[samples_in_cond].mean(axis=0)  # Mean per gene

            for gene in available_genes:
                plot_data.append(
                    {
                        "Gene": gene,
                        "Condition": condition,
                        "Expression": mean_expr[gene],
                    }
                )

        plot_df = pd.DataFrame(plot_data)

        # Create grouped bar chart
        fig = go.Figure()

        for condition in sorted(set(sample_conditions.values())):
            cond_data = plot_df[plot_df["Condition"] == condition]
            fig.add_trace(
                go.Bar(
                    name=condition,
                    x=cond_data["Gene"],
                    y=cond_data["Expression"],
                    text=cond_data["Expression"].round(2),
                    textposition="auto",
                )
            )

        # Update layout
        fig.update_layout(
            title=f"Gene Panel: {panel_name}",
            xaxis_title="Gene",
            yaxis_title="Mean log2(normalized counts + 1)",
            barmode="group",
            template="plotly_white",
            hovermode="x unified",
            height=500,
        )

        return fig
