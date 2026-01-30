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
    
    def get_all_interpretations(self) -> List[Interpretation]:
        """Return all accumulated interpretations."""
        return self.interpretations
    
    def clear(self):
        """Clear all accumulated interpretations."""
        self.interpretations = []
