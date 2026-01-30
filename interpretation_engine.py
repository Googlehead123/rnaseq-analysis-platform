"""
Interpretation Engine for RNA-seq Analysis Platform.

Generates automated, publication-ready interpretations from DE and enrichment results.
"""

from dataclasses import dataclass, field
from typing import List, Dict, Any, Optional
import pandas as pd
import numpy as np


@dataclass
class Interpretation:
    """Single interpretation or insight."""
    category: str       # 'de', 'enrichment', 'qc', 'deconvolution'
    level: str          # 'critical', 'important', 'informative'
    title: str
    description: str
    evidence: Dict[str, Any]
    recommendations: List[str] = field(default_factory=list)


class InterpretationEngine:
    """
    Automated interpretation of RNA-seq results.
    
    Generates human-readable insights from DE results, enrichment results,
    and other analyses. Designed for PhD-level biological interpretation.
    """
    
    def __init__(self):
        self.interpretations: List[Interpretation] = []
    
    def interpret_de_results(
        self,
        de_results_df: pd.DataFrame,
        contrast_name: str,
        padj_threshold: float = 0.05,
        lfc_threshold: float = 1.0
    ) -> List[Interpretation]:
        """
        Generate interpretations from differential expression results.
        
        Args:
            de_results_df: DataFrame with columns: gene, log2FoldChange, padj, pvalue, baseMean
            contrast_name: Human-readable comparison name (e.g., "Treatment vs Control")
            padj_threshold: Adjusted p-value cutoff
            lfc_threshold: Absolute log2 fold change cutoff
        
        Returns:
            List of Interpretation objects
        """
        interpretations = []
        
        # Filter valid rows
        df = de_results_df.dropna(subset=['padj', 'log2FoldChange']).copy()
        
        n_genes = len(df)
        n_up = int(((df['padj'] < padj_threshold) & (df['log2FoldChange'] > lfc_threshold)).sum())
        n_down = int(((df['padj'] < padj_threshold) & (df['log2FoldChange'] < -lfc_threshold)).sum())
        n_sig = n_up + n_down
        
        # 1. Summary interpretation
        interpretations.append(Interpretation(
            category='de',
            level='important',
            title='Differential Expression Summary',
            description=(
                f"Analysis of {n_genes:,} genes in {contrast_name} revealed "
                f"{n_sig:,} significantly differentially expressed genes. "
                f"Of these, {n_up:,} were upregulated and {n_down:,} were downregulated "
                f"(|log₂FC| > {lfc_threshold}, padj < {padj_threshold})."
            ),
            evidence={
                'total_genes_tested': n_genes,
                'significant': n_sig,
                'upregulated': n_up,
                'downregulated': n_down,
                'thresholds': {'padj': padj_threshold, 'lfc': lfc_threshold}
            },
            recommendations=[
                "Examine top DE genes for biological relevance",
                "Perform pathway enrichment analysis on significant genes",
                "Validate key findings with qPCR or Western blot"
            ]
        ))
        
        # 2. Effect size assessment
        sig_mask = df['padj'] < padj_threshold
        if sig_mask.any():
            median_lfc = df.loc[sig_mask, 'log2FoldChange'].abs().median()
            
            if median_lfc < 0.5:
                interpretations.append(Interpretation(
                    category='de',
                    level='informative',
                    title='Modest Effect Sizes Detected',
                    description=(
                        f"The median absolute log₂ fold change among significant genes is {median_lfc:.2f}, "
                        f"indicating relatively modest effect sizes. This may suggest subtle but coordinated "
                        f"biological changes. Consider Gene Set Enrichment Analysis (GSEA) to detect "
                        f"pathway-level shifts that individual gene testing may miss."
                    ),
                    evidence={'median_abs_lfc': float(median_lfc)},
                    recommendations=[
                        "Run GSEA with full ranked gene list",
                        "Consider lowering LFC threshold for exploratory analysis"
                    ]
                ))
            elif median_lfc > 3.0:
                interpretations.append(Interpretation(
                    category='de',
                    level='important',
                    title='Large Effect Sizes Detected',
                    description=(
                        f"The median absolute log₂ fold change among significant genes is {median_lfc:.2f}, "
                        f"indicating strong transcriptional changes. This suggests a robust biological response."
                    ),
                    evidence={'median_abs_lfc': float(median_lfc)},
                    recommendations=[
                        "Prioritize top DE genes for functional validation",
                        "Check for potential technical artifacts in samples with extreme changes"
                    ]
                ))
        
        # 3. Balance assessment
        if n_sig > 0:
            ratio = n_up / n_sig if n_sig > 0 else 0.5
            if ratio > 0.75:
                interpretations.append(Interpretation(
                    category='de',
                    level='informative',
                    title='Predominantly Upregulated Response',
                    description=(
                        f"{n_up}/{n_sig} ({ratio*100:.0f}%) of significant genes are upregulated, "
                        f"suggesting a predominantly activating transcriptional response."
                    ),
                    evidence={'up_ratio': float(ratio), 'up': n_up, 'down': n_down},
                    recommendations=["Investigate activated pathways and transcription factors"]
                ))
            elif ratio < 0.25:
                interpretations.append(Interpretation(
                    category='de',
                    level='informative',
                    title='Predominantly Downregulated Response',
                    description=(
                        f"{n_down}/{n_sig} ({(1-ratio)*100:.0f}%) of significant genes are downregulated, "
                        f"suggesting a predominantly repressive transcriptional response."
                    ),
                    evidence={'down_ratio': float(1 - ratio), 'up': n_up, 'down': n_down},
                    recommendations=["Investigate repressed pathways and regulatory mechanisms"]
                ))
        
        # 4. No significant genes warning
        if n_sig == 0:
            interpretations.append(Interpretation(
                category='de',
                level='critical',
                title='No Significant Genes Detected',
                description=(
                    f"No genes passed the significance thresholds "
                    f"(padj < {padj_threshold}, |log₂FC| > {lfc_threshold}). "
                    f"This may indicate: (1) no biological effect, (2) insufficient statistical power, "
                    f"(3) high variability between replicates, or (4) overly stringent thresholds."
                ),
                evidence={'total_genes': n_genes, 'thresholds': {'padj': padj_threshold, 'lfc': lfc_threshold}},
                recommendations=[
                    "Check sample QC and replicate consistency",
                    "Consider relaxing thresholds (padj < 0.1, |LFC| > 0.5)",
                    "Try GSEA which does not require a significance cutoff",
                    "Verify experimental design and sample labeling"
                ]
            ))
        
        self.interpretations.extend(interpretations)
        return interpretations
    
    def interpret_enrichment(
        self,
        enrichment_df: pd.DataFrame,
        top_n: int = 10
    ) -> List[Interpretation]:
        """
        Generate interpretations from pathway enrichment results.
        
        Args:
            enrichment_df: DataFrame with columns: Term, Adjusted P-value (from Enrichr format)
            top_n: Number of top terms to analyze
        
        Returns:
            List of Interpretation objects
        """
        interpretations = []
        
        if enrichment_df is None or enrichment_df.empty:
            interpretations.append(Interpretation(
                category='enrichment',
                level='important',
                title='No Significant Enrichment',
                description="No significantly enriched pathways were identified.",
                evidence={},
                recommendations=[
                    "Try GSEA with full ranked gene list",
                    "Use different enrichment databases",
                    "Consider relaxing significance thresholds"
                ]
            ))
            self.interpretations.extend(interpretations)
            return interpretations
        
        n_terms = len(enrichment_df)
        interpretations.append(Interpretation(
            category='enrichment',
            level='informative',
            title='Enrichment Overview',
            description=f"Found {n_terms} significantly enriched terms/pathways.",
            evidence={'n_enriched_terms': n_terms},
            recommendations=["Review top pathways for biological relevance"]
        ))
        
        # Keyword-based biological category detection
        term_col = 'Term' if 'Term' in enrichment_df.columns else enrichment_df.columns[0]
        terms_text = ' '.join(enrichment_df[term_col].astype(str).tolist()).lower()
        
        # Immune detection
        immune_kw = ['immune', 'inflammatory', 'cytokine', 't cell', 'b cell', 
                     'leukocyte', 'lymphocyte', 'interferon', 'interleukin',
                     'nf-kb', 'toll-like', 'chemokine', 'complement']
        immune_hits = [kw for kw in immune_kw if kw in terms_text]
        if immune_hits:
            interpretations.append(Interpretation(
                category='enrichment',
                level='important',
                title='Immune/Inflammatory Response Detected',
                description=(
                    f"Enrichment analysis reveals significant immune/inflammatory pathway involvement. "
                    f"Detected keywords: {', '.join(immune_hits[:5])}."
                ),
                evidence={'immune_keywords': immune_hits},
                recommendations=[
                    "Consider cell type deconvolution analysis",
                    "Investigate specific cytokine/chemokine profiles",
                    "Check for immune cell infiltration markers"
                ]
            ))
        
        # Cell cycle / proliferation detection
        prolif_kw = ['cell cycle', 'mitotic', 'dna replication', 'proliferation', 
                     'g1/s', 'g2/m', 'checkpoint', 'mitosis']
        prolif_hits = [kw for kw in prolif_kw if kw in terms_text]
        if prolif_hits:
            interpretations.append(Interpretation(
                category='enrichment',
                level='important',
                title='Cell Cycle/Proliferation Signatures',
                description=(
                    f"Cell cycle and proliferation pathways are enriched, suggesting "
                    f"active cell division or growth regulation changes."
                ),
                evidence={'proliferation_keywords': prolif_hits},
                recommendations=[
                    "Examine key cell cycle regulators (CDK, Cyclin genes)",
                    "Consider Ki-67 or PCNA as proliferation markers"
                ]
            ))
        
        # Metabolic detection
        metab_kw = ['metabolic', 'metabolism', 'glycolysis', 'oxidative phosphorylation',
                    'fatty acid', 'lipid', 'amino acid', 'tca cycle']
        metab_hits = [kw for kw in metab_kw if kw in terms_text]
        if metab_hits:
            interpretations.append(Interpretation(
                category='enrichment',
                level='informative',
                title='Metabolic Pathway Changes',
                description=(
                    f"Metabolic pathways are enriched, indicating potential metabolic reprogramming."
                ),
                evidence={'metabolic_keywords': metab_hits},
                recommendations=[
                    "Investigate specific metabolic shifts",
                    "Consider metabolomics validation"
                ]
            ))
        
        # ECM / structural detection (relevant for dermatology)
        ecm_kw = ['extracellular matrix', 'collagen', 'ecm', 'cell adhesion',
                  'integrin', 'focal adhesion', 'basement membrane']
        ecm_hits = [kw for kw in ecm_kw if kw in terms_text]
        if ecm_hits:
            interpretations.append(Interpretation(
                category='enrichment',
                level='important',
                title='ECM/Structural Pathway Changes',
                description=(
                    f"Extracellular matrix and structural pathways are enriched, "
                    f"suggesting tissue remodeling or structural changes."
                ),
                evidence={'ecm_keywords': ecm_hits},
                recommendations=[
                    "Examine collagen and MMP gene expression",
                    "Consider tissue staining validation (e.g., Masson's trichrome)"
                ]
            ))
        
        self.interpretations.extend(interpretations)
        return interpretations
    
    def generate_methods_text(self, params: Dict) -> str:
        """
        Generate publication-ready methods section text.
        
        Args:
            params: Dictionary with analysis parameters:
                - de_tool: DE analysis tool name (default: 'PyDESeq2')
                - design: Design formula (default: '~ condition')
                - padj: P-value threshold (default: 0.05)
                - lfc: LFC threshold (default: 1.0)
                - quantification: Quantification method (default: 'gene-level counts')
                - enrichment_method: Enrichment method (default: 'over-representation analysis')
                - databases: List of databases used (default: ['GO Biological Process', 'KEGG'])
                - n_samples: Number of samples (optional)
                - n_conditions: Number of conditions (optional)
        
        Returns:
            Formatted methods text as markdown string
        """
        de_tool = params.get('de_tool', 'PyDESeq2')
        design = params.get('design', '~ condition')
        padj = params.get('padj', 0.05)
        lfc = params.get('lfc', 1.0)
        quant = params.get('quantification', 'gene-level counts')
        enr_method = params.get('enrichment_method', 'over-representation analysis (ORA)')
        databases = params.get('databases', ['GO Biological Process', 'KEGG'])
        
        methods = f"""## Methods

### Differential Expression Analysis

Differential gene expression analysis was performed using {de_tool} with the design formula: `{design}`. Raw {quant} were used as input. Genes with fewer than 10 total counts across all samples were filtered prior to analysis. Size factor normalization and dispersion estimation were performed using the {de_tool} pipeline. Genes with Benjamini-Hochberg adjusted p-value < {padj} and absolute log₂ fold change > {lfc} were considered significantly differentially expressed.

### Pathway Enrichment Analysis

Functional enrichment analysis was conducted using {enr_method} via the Enrichr API. Significantly differentially expressed genes were tested against the following databases: {', '.join(databases)}. P-values were corrected for multiple testing using the Benjamini-Hochberg method. Enriched terms with adjusted p-value < 0.05 were reported.

### Visualization and Statistical Analysis

All analyses were performed in Python. Interactive visualizations were generated using Plotly. Principal component analysis (PCA) was performed on log₂-transformed normalized counts using scikit-learn. Hierarchical clustering was performed using Ward's method with Euclidean distance. Statistical significance was determined at α = {padj} unless otherwise noted.
"""
        return methods

    def interpret_top_genes(
        self,
        de_results_df: pd.DataFrame,
        contrast_name: str,
        top_n: int = 10
    ) -> List[Interpretation]:
        """
        Identify top differentially expressed genes with dermatology panel highlights.
        """
        interpretations = []

        if de_results_df is None or de_results_df.empty:
            interpretations.append(Interpretation(
                category='de',
                level='informative',
                title='Top DE Genes Unavailable',
                description="No DE results were provided to identify top genes.",
                evidence={},
                recommendations=["Provide DE results with gene, log2FoldChange, padj columns"]
            ))
            self.interpretations.extend(interpretations)
            return interpretations

        df = de_results_df.copy()
        required_cols = {'padj', 'log2FoldChange'}
        missing = required_cols - set(df.columns)
        if missing:
            interpretations.append(Interpretation(
                category='de',
                level='informative',
                title='Top DE Genes Unavailable',
                description=(
                    "Required DE columns are missing: "
                    f"{', '.join(sorted(missing))}."
                ),
                evidence={'missing_columns': sorted(missing)},
                recommendations=["Ensure DE results include padj and log2FoldChange"]
            ))
            self.interpretations.extend(interpretations)
            return interpretations

        if 'gene' not in df.columns:
            df = df.reset_index().rename(columns={df.index.name or 'index': 'gene'})

        df = df.dropna(subset=['gene', 'padj', 'log2FoldChange'])
        if df.empty:
            interpretations.append(Interpretation(
                category='de',
                level='informative',
                title='Top DE Genes Unavailable',
                description="No valid DE rows were available after filtering missing values.",
                evidence={},
                recommendations=["Check DE results for missing values"]
            ))
            self.interpretations.extend(interpretations)
            return interpretations

        sig_df = df[df['padj'] < 0.05].copy()
        if sig_df.empty:
            interpretations.append(Interpretation(
                category='de',
                level='informative',
                title='No Significant Top Genes',
                description=(
                    f"No genes passed padj < 0.05 in {contrast_name}; "
                    "top gene listing was not generated."
                ),
                evidence={'padj_threshold': 0.05},
                recommendations=["Consider relaxed thresholds or GSEA"]
            ))
            self.interpretations.extend(interpretations)
            return interpretations

        sig_df = sig_df.assign(abs_lfc=sig_df['log2FoldChange'].abs())
        top_df = sig_df.sort_values('abs_lfc', ascending=False).head(top_n)

        up_df = top_df[top_df['log2FoldChange'] > 0]
        down_df = top_df[top_df['log2FoldChange'] < 0]

        def format_gene(row: pd.Series) -> str:
            return f"{row['gene']} (log2FC={row['log2FoldChange']:.2f}, padj={row['padj']:.3g})"

        top_up = [format_gene(row) for _, row in up_df.iterrows()]
        top_down = [format_gene(row) for _, row in down_df.iterrows()]

        collagen_genes = {
            'COL1A1', 'COL1A2', 'COL3A1', 'ELN', 'FBN1', 'MMP1', 'MMP3', 'MMP9', 'TIMP1'
        }
        barrier_genes = {
            'FLG', 'LOR', 'IVL', 'KRT1', 'KRT10', 'KRT14', 'CLDN1', 'TJP1'
        }
        inflammation_genes = {
            'IL1A', 'IL1B', 'IL6', 'IL8', 'CXCL8', 'TNF', 'PTGS2', 'NFKB1'
        }
        melanogenesis_genes = {
            'MITF', 'TYR', 'TYRP1', 'DCT', 'PMEL', 'MC1R'
        }

        top_genes_upper = {str(g).upper() for g in top_df['gene']}
        panel_hits = {
            'collagen': sorted(top_genes_upper & collagen_genes),
            'barrier': sorted(top_genes_upper & barrier_genes),
            'inflammation': sorted(top_genes_upper & inflammation_genes),
            'melanogenesis': sorted(top_genes_upper & melanogenesis_genes)
        }
        panel_hits = {k: v for k, v in panel_hits.items() if v}

        description_parts = [
            f"Top upregulated: {', '.join(top_up) if top_up else 'None'}",
            f"Top downregulated: {', '.join(top_down) if top_down else 'None'}"
        ]
        if panel_hits:
            panel_text = '; '.join(
                f"{panel}: {', '.join(genes)}" for panel, genes in panel_hits.items()
            )
            description_parts.append(f"Dermatology panel hits: {panel_text}.")

        interpretations.append(Interpretation(
            category='de',
            level='important',
            title='Top Differentially Expressed Genes',
            description=(
                f"{contrast_name} top DE genes show the strongest transcriptional shifts. "
                + ' '.join(description_parts)
            ),
            evidence={
                'top_n': top_n,
                'n_significant': int(sig_df.shape[0]),
                'top_upregulated': up_df[['gene', 'log2FoldChange', 'padj']].to_dict(orient='records'),
                'top_downregulated': down_df[['gene', 'log2FoldChange', 'padj']].to_dict(orient='records'),
                'dermatology_panel_hits': panel_hits
            },
            recommendations=[
                "Review top genes for known pathway associations",
                "Cross-check panel hits with phenotype-specific readouts"
            ]
        ))

        self.interpretations.extend(interpretations)
        return interpretations

    def interpret_volcano_pattern(
        self,
        de_results_df: pd.DataFrame,
        contrast_name: str
    ) -> List[Interpretation]:
        """
        Interpret volcano plot patterns from DE results.
        """
        interpretations = []

        if de_results_df is None or de_results_df.empty:
            interpretations.append(Interpretation(
                category='de',
                level='informative',
                title='Volcano Pattern Unavailable',
                description="No DE results were provided to interpret volcano patterns.",
                evidence={},
                recommendations=["Provide DE results with padj and log2FoldChange"]
            ))
            self.interpretations.extend(interpretations)
            return interpretations

        df = de_results_df.copy()
        if 'gene' not in df.columns:
            df = df.reset_index().rename(columns={df.index.name or 'index': 'gene'})

        required_cols = {'padj', 'log2FoldChange'}
        missing = required_cols - set(df.columns)
        if missing:
            interpretations.append(Interpretation(
                category='de',
                level='informative',
                title='Volcano Pattern Unavailable',
                description=(
                    "Required DE columns are missing: "
                    f"{', '.join(sorted(missing))}."
                ),
                evidence={'missing_columns': sorted(missing)},
                recommendations=["Ensure DE results include padj and log2FoldChange"]
            ))
            self.interpretations.extend(interpretations)
            return interpretations

        df = df.dropna(subset=['padj', 'log2FoldChange'])
        sig_df = df[df['padj'] < 0.05]

        n_sig = int(sig_df.shape[0])
        n_up = int((sig_df['log2FoldChange'] > 0).sum())
        n_down = int((sig_df['log2FoldChange'] < 0).sum())

        if n_sig < 10:
            pattern = 'minimal'
            pattern_text = 'Minimal differential signal (fewer than 10 significant genes).'
        else:
            up_ratio = n_up / n_sig if n_sig > 0 else 0
            if 0.4 <= up_ratio <= 0.6:
                pattern = 'symmetric'
                pattern_text = 'Symmetric volcano pattern with balanced up/down regulation.'
            elif up_ratio > 0.7:
                pattern = 'skewed_up'
                pattern_text = 'Skewed toward upregulation with a predominance of induced genes.'
            elif up_ratio < 0.3:
                pattern = 'skewed_down'
                pattern_text = 'Skewed toward downregulation with a predominance of repressed genes.'
            else:
                pattern = 'moderately_skewed'
                pattern_text = 'Moderately skewed volcano pattern with an imbalance of up/down genes.'

        outlier_lfc = df[df['log2FoldChange'].abs() > 5]
        extreme_sig = df[df['padj'].clip(lower=1e-300).apply(lambda x: -np.log10(x)) > 50]

        outlier_genes = outlier_lfc['gene'].astype(str).tolist()[:10]
        extreme_sig_genes = extreme_sig['gene'].astype(str).tolist()[:10]

        description_parts = [
            f"{contrast_name} shows {n_sig} significant genes ({n_up} up, {n_down} down).",
            pattern_text
        ]
        if outlier_genes:
            description_parts.append(
                f"Extreme log2FC genes (|LFC| > 5): {', '.join(outlier_genes)}."
            )
        if extreme_sig_genes:
            description_parts.append(
                "Extremely significant genes (-log10(padj) > 50): "
                f"{', '.join(extreme_sig_genes)}."
            )

        interpretations.append(Interpretation(
            category='de',
            level='informative',
            title='Volcano Plot Pattern',
            description=' '.join(description_parts),
            evidence={
                'significant_genes': n_sig,
                'upregulated': n_up,
                'downregulated': n_down,
                'pattern': pattern,
                'outlier_lfc_genes': outlier_genes,
                'extreme_significance_genes': extreme_sig_genes
            },
            recommendations=[
                "Inspect extreme genes for potential technical artifacts",
                "Review balance of up/down pathways for biological context"
            ]
        ))

        self.interpretations.extend(interpretations)
        return interpretations

    def interpret_pca(
        self,
        pca_coordinates: pd.DataFrame,
        sample_metadata: Dict[str, Dict[str, Any]],
        explained_variance: Optional[Any] = None
    ) -> List[Interpretation]:
        """
        Interpret PCA separation and outliers.
        """
        interpretations = []

        if pca_coordinates is None or pca_coordinates.empty:
            interpretations.append(Interpretation(
                category='qc',
                level='informative',
                title='PCA Results Unavailable',
                description="No PCA coordinates were provided for interpretation.",
                evidence={},
                recommendations=["Ensure PCA coordinates are generated"]
            ))
            self.interpretations.extend(interpretations)
            return interpretations

        if 'PC1' not in pca_coordinates.columns or 'PC2' not in pca_coordinates.columns:
            interpretations.append(Interpretation(
                category='qc',
                level='informative',
                title='PCA Results Unavailable',
                description="PCA coordinates must include PC1 and PC2 columns.",
                evidence={'available_columns': list(pca_coordinates.columns)},
                recommendations=["Provide PCA output with PC1 and PC2 columns"]
            ))
            self.interpretations.extend(interpretations)
            return interpretations

        coords = pca_coordinates[['PC1', 'PC2']].copy()
        coords['sample'] = coords.index.astype(str)
        if sample_metadata:
            coords['condition'] = coords['sample'].map(
                lambda s: sample_metadata.get(s, {}).get('condition')
            )
        else:
            coords['condition'] = None

        valid_coords = coords.dropna(subset=['PC1', 'PC2'])
        if valid_coords.empty:
            interpretations.append(Interpretation(
                category='qc',
                level='informative',
                title='PCA Results Unavailable',
                description="No valid PCA coordinates were available after filtering.",
                evidence={},
                recommendations=["Check PCA output for missing values"]
            ))
            self.interpretations.extend(interpretations)
            return interpretations

        condition_groups = valid_coords.dropna(subset=['condition'])
        conditions = condition_groups['condition'].unique().tolist()

        centroids = None
        within_pc1 = np.nan
        within_pc2 = np.nan
        between_pc1 = np.nan
        between_pc2 = np.nan
        outliers = []

        if len(conditions) >= 1:
            centroids = condition_groups.groupby('condition')[['PC1', 'PC2']].mean()
            within_pc1 = condition_groups.groupby('condition')['PC1'].std(ddof=0).mean()
            within_pc2 = condition_groups.groupby('condition')['PC2'].std(ddof=0).mean()
            if len(conditions) >= 2:
                between_pc1 = centroids['PC1'].std(ddof=0)
                between_pc2 = centroids['PC2'].std(ddof=0)

            for condition, group in condition_groups.groupby('condition'):
                center = centroids.loc[condition]
                distances = np.sqrt(((group[['PC1', 'PC2']] - center) ** 2).sum(axis=1))
                if distances.std(ddof=0) > 0:
                    z_scores = (distances - distances.mean()) / distances.std(ddof=0)
                    outliers.extend(group.loc[z_scores > 2, 'sample'].tolist())

        def safe_ratio(between: float, within: float) -> float:
            if within is None or np.isnan(within):
                return np.nan
            if within == 0:
                return np.inf if between is not None and not np.isnan(between) and between > 0 else np.nan
            if between is None or np.isnan(between):
                return np.nan
            return between / within

        pc1_ratio = safe_ratio(between_pc1, within_pc1)
        pc2_ratio = safe_ratio(between_pc2, within_pc2)

        pc1_separates = bool(pc1_ratio > 1.0) if not np.isnan(pc1_ratio) else False
        pc2_separates = bool(pc2_ratio > 1.0) if not np.isnan(pc2_ratio) else False

        var_pc1 = None
        var_pc2 = None
        if explained_variance is not None:
            if isinstance(explained_variance, dict):
                var_pc1 = explained_variance.get('PC1') or explained_variance.get('pc1')
                var_pc2 = explained_variance.get('PC2') or explained_variance.get('pc2')
            elif isinstance(explained_variance, (list, tuple, np.ndarray, pd.Series)):
                if len(explained_variance) >= 2:
                    var_pc1 = explained_variance[0]
                    var_pc2 = explained_variance[1]

            if var_pc1 is not None and var_pc1 <= 1:
                var_pc1 *= 100
            if var_pc2 is not None and var_pc2 <= 1:
                var_pc2 *= 100

        variance_text = ""
        if var_pc1 is not None:
            variance_text = f" ({var_pc1:.1f}% variance)"

        if pc1_separates or pc2_separates:
            if pc1_separates and not pc2_separates:
                separation_text = f"Conditions separate clearly on PC1{variance_text}."
            elif pc2_separates and not pc1_separates:
                pc2_text = f" ({var_pc2:.1f}% variance)" if var_pc2 is not None else ""
                separation_text = f"Conditions separate clearly on PC2{pc2_text}."
            else:
                pc1_text = f"PC1{variance_text}" if variance_text else "PC1"
                pc2_text = f"PC2 ({var_pc2:.1f}% variance)" if var_pc2 is not None else "PC2"
                separation_text = f"Conditions separate across {pc1_text} and {pc2_text}."
            level = 'informative'
        else:
            separation_text = "No clear separation — biological effect may be subtle."
            level = 'important'

        description_parts = [separation_text]
        if outliers:
            description_parts.append(f"Outlier samples: {', '.join(sorted(set(outliers)))}.")
        if not conditions:
            description_parts.append("Condition labels were unavailable for group separation analysis.")

        interpretations.append(Interpretation(
            category='qc',
            level=level,
            title='PCA Separation Assessment',
            description=' '.join(description_parts),
            evidence={
                'conditions': conditions,
                'pc1_ratio': None if np.isnan(pc1_ratio) else float(pc1_ratio),
                'pc2_ratio': None if np.isnan(pc2_ratio) else float(pc2_ratio),
                'outlier_samples': sorted(set(outliers)),
                'explained_variance': {'PC1': var_pc1, 'PC2': var_pc2}
            },
            recommendations=[
                "Review outlier samples for QC issues",
                "Confirm metadata assignments if separation is poor"
            ]
        ))

        self.interpretations.extend(interpretations)
        return interpretations

    def interpret_gsea_results(self, gsea_df: pd.DataFrame) -> List[Interpretation]:
        """
        Interpret GSEA results with directionality.
        """
        interpretations = []

        if gsea_df is None or gsea_df.empty:
            interpretations.append(Interpretation(
                category='enrichment',
                level='important',
                title='No Significant GSEA Results',
                description="GSEA did not return significant pathways for interpretation.",
                evidence={},
                recommendations=["Check ranked gene list and gene set selection"]
            ))
            self.interpretations.extend(interpretations)
            return interpretations

        df = gsea_df.copy()
        term_col = 'Term' if 'Term' in df.columns else df.columns[0]
        nes_col = next((c for c in df.columns if c.lower() == 'nes' or 'nes' in c.lower()), None)
        fdr_col = next(
            (c for c in df.columns if 'fdr' in c.lower() or 'q-val' in c.lower() or 'q value' in c.lower()),
            None
        )

        if nes_col is None or fdr_col is None:
            interpretations.append(Interpretation(
                category='enrichment',
                level='important',
                title='GSEA Results Incomplete',
                description="GSEA results are missing NES or FDR columns needed for interpretation.",
                evidence={'available_columns': list(df.columns)},
                recommendations=["Ensure GSEA results include NES and FDR q-value"]
            ))
            self.interpretations.extend(interpretations)
            return interpretations

        df = df.dropna(subset=[term_col, nes_col, fdr_col])
        sig_df = df[df[fdr_col] < 0.25]

        activated = sig_df[sig_df[nes_col] > 0].sort_values(nes_col, ascending=False)
        suppressed = sig_df[sig_df[nes_col] < 0].sort_values(nes_col, ascending=True)

        top_activated = activated[term_col].astype(str).tolist()[:3]
        top_suppressed = suppressed[term_col].astype(str).tolist()[:3]

        description_parts = []
        if top_activated:
            description_parts.append(f"Activated pathways: {', '.join(top_activated)}.")
        if top_suppressed:
            description_parts.append(f"Suppressed pathways: {', '.join(top_suppressed)}.")
        if not description_parts:
            description_parts.append("No pathways passed FDR < 0.25 for directional interpretation.")

        concordant_themes = []
        if self.interpretations:
            ora_theme_map = {
                'immune_keywords': 'immune',
                'proliferation_keywords': 'cell cycle',
                'metabolic_keywords': 'metabolic',
                'ecm_keywords': 'ecm'
            }
            ora_themes = set()
            for interp in self.interpretations:
                if interp.category == 'enrichment':
                    for key, theme in ora_theme_map.items():
                        if key in interp.evidence:
                            ora_themes.add(theme)

            if ora_themes:
                gsea_terms_text = ' '.join(df[term_col].astype(str).tolist()).lower()
                theme_keywords = {
                    'immune': ['immune', 'interferon', 'cytokine', 't cell', 'b cell'],
                    'cell cycle': ['cell cycle', 'mitotic', 'proliferation', 'checkpoint'],
                    'metabolic': ['metabolic', 'metabolism', 'glycolysis', 'oxidative'],
                    'ecm': ['extracellular matrix', 'collagen', 'ecm', 'focal adhesion']
                }
                for theme in ora_themes:
                    if any(kw in gsea_terms_text for kw in theme_keywords.get(theme, [])):
                        concordant_themes.append(theme)

        if concordant_themes:
            description_parts.append(
                "GSEA themes align with ORA findings: " + ', '.join(sorted(set(concordant_themes))) + "."
            )

        interpretations.append(Interpretation(
            category='enrichment',
            level='important',
            title='GSEA Directional Pathways',
            description=' '.join(description_parts),
            evidence={
                'n_activated': int(activated.shape[0]),
                'n_suppressed': int(suppressed.shape[0]),
                'top_activated': top_activated,
                'top_suppressed': top_suppressed,
                'fdr_threshold': 0.25,
                'concordant_themes': sorted(set(concordant_themes))
            },
            recommendations=[
                "Inspect leading-edge genes for key drivers",
                "Compare with ORA results for consistent pathway signals"
            ]
        ))

        self.interpretations.extend(interpretations)
        return interpretations

    def interpret_deconvolution(
        self,
        deconv_results: pd.DataFrame,
        sample_metadata: Optional[Dict[str, Dict[str, Any]]] = None
    ) -> List[Interpretation]:
        """
        Interpret cell type deconvolution results.
        """
        interpretations = []

        if deconv_results is None or deconv_results.empty:
            interpretations.append(Interpretation(
                category='deconvolution',
                level='informative',
                title='Deconvolution Results Unavailable',
                description="No deconvolution results were provided for interpretation.",
                evidence={},
                recommendations=["Run cell type deconvolution before interpretation"]
            ))
            self.interpretations.extend(interpretations)
            return interpretations

        df = deconv_results.copy()
        df = df.dropna(how='all')
        if df.empty:
            interpretations.append(Interpretation(
                category='deconvolution',
                level='informative',
                title='Deconvolution Results Unavailable',
                description="Deconvolution results contained only missing values.",
                evidence={},
                recommendations=["Check deconvolution output for validity"]
            ))
            self.interpretations.extend(interpretations)
            return interpretations

        if sample_metadata:
            sample_conditions = {
                str(sample): meta.get('condition')
                for sample, meta in sample_metadata.items()
            }
            df['condition'] = df.index.astype(str).map(sample_conditions).fillna('Unknown')
        else:
            df['condition'] = 'All samples'

        dominant_cells = {}
        condition_means = {}
        for condition, group in df.groupby('condition'):
            cell_df = group.drop(columns=['condition'], errors='ignore')
            if cell_df.empty:
                continue
            mean_props = cell_df.mean(axis=0)
            if mean_props.empty:
                continue
            top_cells = mean_props.sort_values(ascending=False).head(2)
            dominant_cells[str(condition)] = top_cells.index.astype(str).tolist()
            condition_means[str(condition)] = mean_props

        condition_comparison = []
        if len(condition_means) >= 2:
            mean_df = pd.DataFrame(condition_means).T
            for cell_type in mean_df.columns:
                max_cond = mean_df[cell_type].idxmax()
                min_cond = mean_df[cell_type].idxmin()
                diff = mean_df.loc[max_cond, cell_type] - mean_df.loc[min_cond, cell_type]
                condition_comparison.append((cell_type, max_cond, min_cond, diff))
            condition_comparison.sort(key=lambda x: abs(x[3]), reverse=True)
            condition_comparison = condition_comparison[:3]

        immune_keywords = [
            't cell', 'b cell', 'nk', 'macrophage', 'monocyte', 'neutrophil',
            'dendritic', 'immune', 'lymphocyte', 'mast', 'microglia'
        ]
        immune_cols = [
            col for col in df.columns
            if any(kw in col.lower() for kw in immune_keywords)
        ]
        immune_flags = []
        if immune_cols:
            immune_fraction = df[immune_cols].sum(axis=1)
            immune_flags = immune_fraction[immune_fraction > 0.3].index.astype(str).tolist()

        description_parts = []
        if dominant_cells:
            dominant_text = '; '.join(
                f"{cond}: {', '.join(cells)}" for cond, cells in dominant_cells.items()
            )
            description_parts.append(f"Dominant cell types by condition — {dominant_text}.")
        if condition_comparison:
            diff_text = '; '.join(
                f"{cell} higher in {max_c} vs {min_c} (Δ={diff:.2f})"
                for cell, max_c, min_c, diff in condition_comparison
            )
            description_parts.append(f"Largest composition shifts: {diff_text}.")
        if immune_flags:
            description_parts.append(
                "High immune fractions (>30%) detected in samples: "
                f"{', '.join(immune_flags)}."
            )

        interpretations.append(Interpretation(
            category='deconvolution',
            level='informative',
            title='Cell Type Composition Summary',
            description=' '.join(description_parts) if description_parts else
                "Deconvolution results were provided but no dominant patterns were identified.",
            evidence={
                'dominant_cells': dominant_cells,
                'condition_comparison': [
                    {
                        'cell_type': cell,
                        'higher_in': max_c,
                        'lower_in': min_c,
                        'delta': float(diff)
                    }
                    for cell, max_c, min_c, diff in condition_comparison
                ],
                'immune_high_samples': immune_flags
            },
            recommendations=[
                "Cross-check dominant cell types with histology",
                "Review immune fraction if tissue is expected to be non-immune"
            ]
        ))

        self.interpretations.extend(interpretations)
        return interpretations

    def generate_misinterpretation_warnings(
        self,
        de_results_df: pd.DataFrame,
        enrichment_df: Optional[pd.DataFrame] = None
    ) -> List[Interpretation]:
        """
        Generate contextual warnings to avoid misinterpretation.
        """
        warnings = []

        def add_warning(title: str, description: str, evidence: Optional[Dict[str, Any]] = None):
            warnings.append(Interpretation(
                category='warning',
                level='important',
                title=title,
                description=description,
                evidence=evidence or {},
                recommendations=[]
            ))

        df = None
        if de_results_df is not None and not de_results_df.empty:
            df = de_results_df.copy()
            if 'padj' in df.columns:
                df = df.dropna(subset=['padj'])

        if df is not None and not df.empty and 'padj' in df.columns:
            sig_mask = df['padj'] < 0.05
            n_sig = int(sig_mask.sum())

            borderline_mask = df['padj'].between(0.045, 0.055)
            borderline_count = int(borderline_mask.sum())
            if borderline_count >= 20 or (n_sig >= 10 and borderline_count / max(n_sig, 1) >= 0.3):
                add_warning(
                    title='Borderline Significance',
                    description=(
                        f"A notable fraction of genes fall near the padj=0.05 threshold "
                        f"(n={borderline_count}). Interpret marginal hits cautiously."
                    ),
                    evidence={'borderline_count': borderline_count, 'padj_range': [0.045, 0.055]}
                )

            if n_sig > 500:
                add_warning(
                    title='Large DE List',
                    description=(
                        f"More than 500 genes are significant (n={n_sig}). "
                        "Large DE lists can reduce enrichment specificity."
                    ),
                    evidence={'n_significant': n_sig}
                )

            if 'log2FoldChange' in df.columns and n_sig > 0:
                median_abs_lfc = df.loc[sig_mask, 'log2FoldChange'].abs().median()
                if pd.notna(median_abs_lfc) and median_abs_lfc < 0.5:
                    add_warning(
                        title='Small Effect Sizes',
                        description=(
                            f"Significant genes have small effect sizes (median |LFC|={median_abs_lfc:.2f}). "
                            "Statistical significance may not translate to strong biological impact."
                        ),
                        evidence={'median_abs_lfc': float(median_abs_lfc)}
                    )

            replicate_counts = None
            attrs = getattr(de_results_df, 'attrs', {}) if de_results_df is not None else {}
            if isinstance(attrs, dict):
                for key in ['replicate_counts', 'sample_counts', 'condition_counts']:
                    if isinstance(attrs.get(key), dict):
                        replicate_counts = attrs.get(key)
                        break
                if replicate_counts is None and isinstance(attrs.get('sample_metadata'), dict):
                    replicate_counts = {}
                    for meta in attrs['sample_metadata'].values():
                        condition = meta.get('condition') if isinstance(meta, dict) else None
                        if condition:
                            replicate_counts[condition] = replicate_counts.get(condition, 0) + 1

            if replicate_counts:
                low_rep = {cond: cnt for cond, cnt in replicate_counts.items() if cnt < 3}
                if low_rep:
                    add_warning(
                        title='Low Replicate Count',
                        description=(
                            "Some conditions have fewer than 3 replicates, reducing statistical power: "
                            + ', '.join(f"{cond} (n={cnt})" for cond, cnt in low_rep.items())
                        ),
                        evidence={'low_replicate_conditions': low_rep}
                    )

        if enrichment_df is not None and not enrichment_df.empty:
            adj_col = next(
                (
                    c for c in enrichment_df.columns
                    if 'adjusted' in c.lower() or 'fdr' in c.lower() or 'q-val' in c.lower()
                    or c.lower() in {'padj', 'qvalue', 'q-value'}
                ),
                None
            )
            if adj_col is not None:
                n_sig_terms = int((enrichment_df[adj_col] < 0.05).sum())
                if n_sig_terms > 100:
                    add_warning(
                        title='Over-enrichment Risk',
                        description=(
                            f"More than 100 pathways are significant (n={n_sig_terms}). "
                            "Consider stricter thresholds to improve specificity."
                        ),
                        evidence={'n_significant_terms': n_sig_terms, 'threshold': 0.05}
                    )

        add_warning(
            title='Correlation ≠ Causation',
            description="Transcriptomic associations do not establish causal relationships."
        )
        add_warning(
            title='Orthogonal Validation Required',
            description="Validate key findings with orthogonal methods (qPCR, protein assays, or functional tests)."
        )

        self.interpretations.extend(warnings)
        return warnings
    
    def get_all_interpretations(self) -> List[Interpretation]:
        """Return all accumulated interpretations."""
        return self.interpretations
    
    def clear(self):
        """Clear all accumulated interpretations."""
        self.interpretations = []
