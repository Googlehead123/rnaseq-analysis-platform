# Bulk RNA-Seq Analysis Platform: Complete Technical Specification
## A PhD-Grade, Publication-Ready Platform Guide

---

# Table of Contents

1. [Platform Overview & Goals](#1-platform-overview--goals)
2. [Core Analysis Modules](#2-core-analysis-modules)
3. [Python Implementation Details](#3-python-implementation-details)
4. [Visualization Specifications](#4-visualization-specifications)
5. [Interpretation Engine](#5-interpretation-engine)
6. [UI Architecture](#6-ui-architecture)
7. [Pipeline Orchestration](#7-pipeline-orchestration)
8. [Project Structure](#8-project-structure)

---

# 1. Platform Overview & Goals

## 1.1 What Bulk RNA-Seq Measures

Bulk RNA-seq quantifies the **transcriptome** of a biological sample—the complete set of RNA transcripts. It captures an **averaged expression profile** across thousands to millions of cells.

### Core Information Extractable

| Data Type | Description | Use Case |
|-----------|-------------|----------|
| **Gene Expression Quantification** | Raw/normalized counts (TPM, FPKM, CPM) for every detected gene | Compare expression levels across samples |
| **Differential Expression** | Genes with significant expression changes between conditions | Identify treatment effects, disease signatures |
| **Pathway/Functional Enrichment** | Over-represented biological processes among DE genes | Understand mechanisms, generate hypotheses |
| **Sample Relationships** | Clustering, batch effects, outliers | QC validation, subtype discovery |
| **Cell Type Deconvolution** | Estimated proportions of cell types | Immune infiltration, tissue composition |

## 1.2 Platform Architecture

```
┌─────────────────────────────────────────────────────────────────┐
│                      USER INTERFACE LAYER                        │
│   ┌──────────────┐  ┌──────────────┐  ┌────────────────────┐   │
│   │ Web UI       │  │ CLI          │  │ Interactive Report │   │
│   │ (Streamlit)  │  │ (Click)      │  │ (Quarto/Jupyter)   │   │
│   └──────────────┘  └──────────────┘  └────────────────────┘   │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│                    ORCHESTRATION LAYER                           │
│   ┌─────────────────────────────────────────────────────────┐   │
│   │ Snakemake / Nextflow                                    │   │
│   │ • DAG-based execution  • Containerization              │   │
│   │ • Checkpointing        • Cloud/HPC support             │   │
│   └─────────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│                      ANALYSIS MODULES                            │
│  ┌────────────┐ ┌────────────┐ ┌────────────┐ ┌────────────┐   │
│  │ QC &       │ │ Quanti-    │ │ Differ-    │ │ Enrichment │   │
│  │ Preproc    │ │ fication   │ │ ential DE  │ │            │   │
│  └────────────┘ └────────────┘ └────────────┘ └────────────┘   │
│  ┌────────────┐ ┌────────────┐ ┌────────────┐ ┌────────────┐   │
│  │ Clustering │ │ Deconvo-   │ │ Visual-    │ │ Interpre-  │   │
│  │            │ │ lution     │ │ ization    │ │ tation     │   │
│  └────────────┘ └────────────┘ └────────────┘ └────────────┘   │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│                        DATA LAYER                                │
│  ┌────────────┐ ┌────────────┐ ┌────────────┐ ┌────────────┐   │
│  │ Raw FASTQ  │ │ Count      │ │ Results    │ │ Reference  │   │
│  │ Storage    │ │ Matrices   │ │ Database   │ │ Annotations│   │
│  └────────────┘ └────────────┘ └────────────┘ └────────────┘   │
└─────────────────────────────────────────────────────────────────┘
```

## 1.3 Technology Stack

| Component | Tool | Rationale |
|-----------|------|-----------|
| **Workflow** | Snakemake | Academic standard, reproducibility |
| **Backend** | FastAPI/Python | Modern, fast, typed |
| **Web UI** | Streamlit | Rapid development, Python-native |
| **Visualization** | Plotly + Matplotlib + Seaborn | Interactive + publication quality |
| **Statistics** | PyDESeq2 + scipy | Python-native DE analysis |
| **Enrichment** | gseapy | GSEA/ORA in Python |
| **Reports** | Quarto/Jupyter | Reproducible documents |

---

# 2. Core Analysis Modules

## 2.1 Quality Control & Preprocessing

### Tools & Functions

| Step | Tool | Key Outputs |
|------|------|-------------|
| Raw read QC | FastQC + MultiQC | Per-sample quality reports |
| Trimming | fastp | Adapter-free, quality-filtered reads |
| Post-alignment QC | RSeQC, Picard | Mapping rates, gene body coverage |

### Metrics to Assess

- **Per-base quality scores** (should be >Q30)
- **GC content** (species-specific, ~40-60% for human)
- **Adapter contamination** (<1% ideal)
- **Sequence duplication** (indicates PCR bias)
- **Mapping rate** (>70% for good libraries)

### Example Implementation

```python
class RawReadQC:
    def run_fastqc(self, fastq_files: List[Path]) -> Path:
        """Run FastQC on FASTQ files."""
        cmd = ["fastqc", "--outdir", str(self.output_dir), 
               "--threads", str(self.threads)] + [str(f) for f in fastq_files]
        subprocess.run(cmd, check=True)
        return self.output_dir
    
    def run_multiqc(self, input_dir: Path) -> Path:
        """Aggregate QC reports."""
        cmd = ["multiqc", str(input_dir), "-o", str(self.output_dir)]
        subprocess.run(cmd, check=True)
        return self.output_dir / "multiqc_report.html"
```

## 2.2 Gene Expression Quantification

### Two Main Approaches

#### Approach A: Pseudo-alignment (Salmon/Kallisto)
**Pros**: Fast, accurate, no BAM files needed
**Cons**: Cannot do variant calling or IGV visualization

```python
class SalmonQuantifier:
    def quantify_sample(self, r1: Path, r2: Path, sample_name: str) -> Path:
        """Quantify with Salmon."""
        cmd = [
            "salmon", "quant",
            "-i", str(self.index_path),
            "-l", "A",  # Auto-detect library type
            "-1", str(r1), "-2", str(r2),
            "-o", str(self.output_dir / sample_name),
            "--gcBias", "--seqBias",  # Bias correction
            "--validateMappings"
        ]
        subprocess.run(cmd, check=True)
        return self.output_dir / sample_name
```

#### Approach B: Alignment-based (STAR + featureCounts)
**Pros**: BAM files for visualization, variant calling
**Cons**: Slower, more disk space

### Normalization Methods

| Method | Description | When to Use |
|--------|-------------|-------------|
| **TPM** | Transcripts per million | Cross-sample comparison |
| **CPM** | Counts per million | Simple library size normalization |
| **TMM** | Trimmed mean of M-values | edgeR's method |
| **RLE** | Relative log expression | DESeq2's method |
| **VST** | Variance stabilizing transform | Visualization, clustering |

## 2.3 Differential Expression Analysis

### Statistical Framework

1. **Model**: Negative binomial distribution (handles overdispersion)
2. **Design**: Generalized linear models for complex designs
3. **Shrinkage**: Empirical Bayes for variance estimation
4. **Correction**: Benjamini-Hochberg FDR

### Key Parameters

```python
@dataclass
class DEConfig:
    padj_threshold: float = 0.05     # Adjusted p-value cutoff
    lfc_threshold: float = 1.0       # |log2FC| cutoff
    min_counts: int = 10             # Minimum total counts per gene
    min_samples: int = 3             # Minimum samples with counts
    shrinkage: bool = True           # Apply LFC shrinkage
    shrinkage_method: str = "apeglm" # apeglm, ashr, or normal
```

### Implementation with PyDESeq2

```python
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

class DifferentialExpression:
    def run_analysis(
        self,
        counts: pd.DataFrame,
        metadata: pd.DataFrame,
        design: str,
        contrast: Tuple[str, str, str]  # (variable, numerator, denominator)
    ) -> DEResult:
        """Run differential expression analysis."""
        
        # Filter low-count genes
        counts_filtered = self._filter_genes(counts)
        
        # Create DESeq dataset
        dds = DeseqDataSet(
            counts=counts_filtered.T,
            metadata=metadata,
            design_factors=design.replace("~", "").strip()
        )
        
        # Run DESeq2 pipeline
        dds.deseq2()
        
        # Get statistics
        stat_res = DeseqStats(dds, contrast=list(contrast))
        stat_res.summary()
        
        # Apply shrinkage
        if self.shrinkage:
            stat_res.lfc_shrink(coeff=f"{contrast[0]}_{contrast[1]}_vs_{contrast[2]}")
        
        return DEResult(
            results_table=stat_res.results_df,
            normalized_counts=dds.layers['normed_counts'],
            ...
        )
```

### Output Columns

| Column | Description |
|--------|-------------|
| `baseMean` | Average normalized count across samples |
| `log2FoldChange` | Effect size (treatment/control) |
| `lfcSE` | Standard error of LFC |
| `stat` | Wald test statistic |
| `pvalue` | Raw p-value |
| `padj` | BH-adjusted p-value |

## 2.4 Pathway & Functional Enrichment

### Methods

#### Over-Representation Analysis (ORA)
- Input: List of DE genes
- Test: Fisher's exact / hypergeometric
- Output: Enriched terms with p-values

#### Gene Set Enrichment Analysis (GSEA)
- Input: All genes ranked by expression change
- Test: Running sum enrichment score
- Output: Normalized enrichment score (NES)

### Databases

| Database | Content | Use Case |
|----------|---------|----------|
| **GO BP** | Biological processes | Mechanism discovery |
| **GO MF** | Molecular functions | Protein activities |
| **GO CC** | Cellular components | Localization |
| **KEGG** | Metabolic pathways | Pathway mapping |
| **Reactome** | Curated pathways | Detailed mechanisms |
| **MSigDB Hallmark** | Cancer hallmarks | Cancer biology |

### Implementation

```python
class OverRepresentationAnalysis:
    def run(
        self,
        query_genes: List[str],
        libraries: List[str],
        background_genes: Optional[Set[str]] = None
    ) -> pd.DataFrame:
        """Run ORA analysis."""
        
        results = []
        for library in libraries:
            gene_sets = self.load_gene_sets(library)
            
            for term_name, term_genes in gene_sets.items():
                # Fisher's exact test
                overlap = set(query_genes) & term_genes
                
                a = len(overlap)
                b = len(query_genes) - a
                c = len(term_genes) - a
                d = len(background) - a - b - c
                
                odds_ratio, p_value = stats.fisher_exact(
                    [[a, b], [c, d]], alternative='greater'
                )
                
                results.append({
                    'term_name': term_name,
                    'p_value': p_value,
                    'odds_ratio': odds_ratio,
                    'overlap_genes': list(overlap)
                })
        
        # Multiple testing correction
        results_df = pd.DataFrame(results)
        results_df['adjusted_p_value'] = multipletests(
            results_df['p_value'], method='fdr_bh'
        )[1]
        
        return results_df.sort_values('adjusted_p_value')
```

## 2.5 Sample Relationship Analysis

### Methods

| Analysis | Purpose | Tool |
|----------|---------|------|
| **PCA** | Dimensionality reduction, QC | scikit-learn |
| **UMAP** | Non-linear embedding | umap-learn |
| **Hierarchical clustering** | Sample grouping | scipy |
| **Distance matrices** | Sample similarity | scipy, sklearn |

### Implementation

```python
class SampleRelationshipAnalyzer:
    def run_pca(self, n_components: int = 50) -> PCAResult:
        """Principal Component Analysis."""
        
        # Use top variable genes
        gene_vars = self.vst_counts.var(axis=1)
        top_genes = gene_vars.nlargest(2000).index
        data = self.vst_counts.loc[top_genes].T
        
        # Fit PCA
        pca = PCA(n_components=n_components)
        embedding = pca.fit_transform(data)
        
        return PCAResult(
            embedding=pd.DataFrame(embedding, index=data.index),
            explained_variance=pca.explained_variance_ratio_,
            loadings=pd.DataFrame(pca.components_.T, index=data.columns)
        )
    
    def detect_outliers(self, pca_result, threshold_sd: float = 3.0) -> List[str]:
        """Detect outlier samples in PC1-PC2 space."""
        pc_data = pca_result.embedding[['PC1', 'PC2']]
        
        centroid = pc_data.mean()
        distances = np.sqrt(((pc_data - centroid) ** 2).sum(axis=1))
        
        threshold = distances.mean() + threshold_sd * distances.std()
        return distances[distances > threshold].index.tolist()
```

## 2.6 Cell Type Deconvolution

### Methods

| Method | Approach | Reference Required |
|--------|----------|-------------------|
| **NNLS** | Non-negative least squares | Signature matrix |
| **CIBERSORTx** | Support vector regression | LM22 signature |
| **EPIC** | Constrained optimization | EPIC signatures |
| **xCell** | Gene signature enrichment | xCell signatures |
| **MuSiC** | Weighted estimation | scRNA-seq reference |

### Implementation

```python
class CellTypeDeconvolution:
    def run_nnls(
        self,
        expression: pd.DataFrame,  # TPM matrix
        signature: pd.DataFrame    # Reference signatures
    ) -> pd.DataFrame:
        """NNLS deconvolution."""
        
        # Find common genes
        common = expression.index.intersection(signature.index)
        expr = expression.loc[common]
        sig = signature.loc[common]
        
        # NNLS for each sample
        proportions = []
        for sample in expr.columns:
            coef, _ = scipy.optimize.nnls(sig.values, expr[sample].values)
            coef_norm = coef / coef.sum()  # Normalize to 1
            proportions.append(coef_norm)
        
        return pd.DataFrame(
            proportions,
            index=expr.columns,
            columns=sig.columns
        )
```

---

# 3. Python Implementation Details

## 3.1 Key Dependencies

```python
# requirements.txt
numpy>=1.21.0
pandas>=1.3.0
scipy>=1.7.0
scikit-learn>=1.0.0

# Statistical analysis
pydeseq2>=0.4.0
statsmodels>=0.13.0
gseapy>=1.0.0

# Visualization
matplotlib>=3.5.0
seaborn>=0.11.0
plotly>=5.0.0

# Dimensionality reduction
umap-learn>=0.5.0

# Web application
streamlit>=1.20.0

# Workflow
snakemake>=7.0.0
```

## 3.2 Data Classes for Results

```python
from dataclasses import dataclass
from typing import Dict, List, Optional
import pandas as pd
import numpy as np

@dataclass
class DEResult:
    """Container for DE results."""
    results_table: pd.DataFrame
    normalized_counts: pd.DataFrame
    size_factors: pd.Series
    design_matrix: pd.DataFrame
    contrasts: Dict[str, str]
    dispersions: pd.Series
    
    def get_significant(
        self,
        padj: float = 0.05,
        lfc: float = 1.0
    ) -> pd.DataFrame:
        """Get significant genes."""
        mask = (
            (self.results_table['padj'] < padj) &
            (self.results_table['log2FoldChange'].abs() > lfc)
        )
        return self.results_table[mask].sort_values('padj')

@dataclass
class EnrichmentResult:
    """Container for enrichment results."""
    results_table: pd.DataFrame
    query_genes: List[str]
    background_size: int
    databases_used: List[str]

@dataclass
class PCAResult:
    """Container for PCA results."""
    embedding: pd.DataFrame
    explained_variance: np.ndarray
    loadings: pd.DataFrame
    n_components: int
```

---

# 4. Visualization Specifications

## 4.1 Publication-Quality Settings

```python
import matplotlib.pyplot as plt

# Global settings for publication
plt.rcParams.update({
    'font.family': 'Arial',
    'font.size': 10,
    'axes.labelsize': 12,
    'axes.titlesize': 14,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 10,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'axes.spines.top': False,
    'axes.spines.right': False,
})
```

## 4.2 Core Visualizations

### Volcano Plot

```python
def plot_volcano(
    de_results: pd.DataFrame,
    lfc_threshold: float = 1.0,
    pval_threshold: float = 0.05,
    label_top_n: int = 10,
    interactive: bool = False
):
    """Create volcano plot for DE results."""
    
    df = de_results.copy()
    df['-log10(padj)'] = -np.log10(df['padj'].clip(lower=1e-300))
    
    # Classify genes
    df['status'] = 'NS'
    df.loc[(df['log2FoldChange'] > lfc_threshold) & 
           (df['padj'] < pval_threshold), 'status'] = 'Up'
    df.loc[(df['log2FoldChange'] < -lfc_threshold) & 
           (df['padj'] < pval_threshold), 'status'] = 'Down'
    
    if interactive:
        import plotly.express as px
        fig = px.scatter(
            df.reset_index(),
            x='log2FoldChange',
            y='-log10(padj)',
            color='status',
            color_discrete_map={'Up': '#e74c3c', 'Down': '#3498db', 'NS': '#95a5a6'},
            hover_data=['index']
        )
        fig.add_vline(x=lfc_threshold, line_dash="dash")
        fig.add_vline(x=-lfc_threshold, line_dash="dash")
        fig.add_hline(y=-np.log10(pval_threshold), line_dash="dash")
        return fig
    
    else:
        fig, ax = plt.subplots(figsize=(8, 6))
        
        colors = {'Up': '#e74c3c', 'Down': '#3498db', 'NS': '#95a5a6'}
        for status, color in colors.items():
            mask = df['status'] == status
            ax.scatter(
                df.loc[mask, 'log2FoldChange'],
                df.loc[mask, '-log10(padj)'],
                c=color, alpha=0.6, s=20, label=f"{status} ({mask.sum()})"
            )
        
        ax.axvline(lfc_threshold, ls='--', c='gray', alpha=0.5)
        ax.axvline(-lfc_threshold, ls='--', c='gray', alpha=0.5)
        ax.axhline(-np.log10(pval_threshold), ls='--', c='gray', alpha=0.5)
        
        ax.set_xlabel('log₂(Fold Change)')
        ax.set_ylabel('-log₁₀(Adjusted P-value)')
        ax.legend()
        
        return fig
```

### PCA Plot

```python
def plot_pca(
    pca_result,
    metadata: pd.DataFrame,
    color_by: str,
    pc_x: int = 1,
    pc_y: int = 2,
    interactive: bool = False
):
    """Create PCA scatter plot."""
    
    embedding = pca_result.embedding
    var = pca_result.explained_variance
    
    x_label = f'PC{pc_x} ({var[pc_x-1]*100:.1f}%)'
    y_label = f'PC{pc_y} ({var[pc_y-1]*100:.1f}%)'
    
    if interactive:
        import plotly.express as px
        df = embedding.join(metadata)
        fig = px.scatter(
            df.reset_index(),
            x=f'PC{pc_x}', y=f'PC{pc_y}',
            color=color_by,
            hover_data=['index'],
            labels={f'PC{pc_x}': x_label, f'PC{pc_y}': y_label}
        )
        return fig
    
    else:
        fig, ax = plt.subplots(figsize=(8, 6))
        
        for group in metadata[color_by].unique():
            mask = metadata[color_by] == group
            samples = metadata[mask].index
            ax.scatter(
                embedding.loc[samples, f'PC{pc_x}'],
                embedding.loc[samples, f'PC{pc_y}'],
                label=group, s=80, alpha=0.7
            )
        
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)
        ax.legend(title=color_by)
        
        return fig
```

### Heatmap with Annotations

```python
def plot_heatmap(
    expression: pd.DataFrame,
    genes: List[str],
    metadata: pd.DataFrame,
    annotation_cols: List[str],
    scale: str = 'row',
    figsize: Tuple = (12, 10)
):
    """Create clustered heatmap with annotations."""
    
    import seaborn as sns
    
    # Subset data
    data = expression.loc[genes, metadata.index]
    
    # Scale
    if scale == 'row':
        data = data.sub(data.mean(axis=1), axis=0).div(data.std(axis=1), axis=0)
    
    # Create annotation colors
    col_colors = []
    for col in annotation_cols:
        values = metadata[col]
        palette = sns.color_palette("Set2", values.nunique())
        color_map = dict(zip(values.unique(), palette))
        col_colors.append([color_map[v] for v in values])
    
    # Create clustermap
    g = sns.clustermap(
        data,
        col_colors=col_colors,
        cmap='RdBu_r',
        center=0,
        figsize=figsize,
        xticklabels=True,
        yticklabels=len(genes) <= 50
    )
    
    return g.fig
```

### Enrichment Dot Plot

```python
def plot_enrichment_dotplot(
    enrichment_results: pd.DataFrame,
    top_n: int = 20,
    interactive: bool = False
):
    """Create enrichment dot plot."""
    
    df = enrichment_results.head(top_n).copy()
    df['-log10(padj)'] = -np.log10(df['adjusted_p_value'].clip(1e-300))
    
    if interactive:
        import plotly.express as px
        fig = px.scatter(
            df,
            x='odds_ratio',
            y='term_name',
            size='overlap_size',
            color='-log10(padj)',
            color_continuous_scale='Viridis'
        )
        return fig
    
    else:
        fig, ax = plt.subplots(figsize=(8, top_n * 0.3 + 2))
        
        scatter = ax.scatter(
            df['odds_ratio'],
            range(len(df)),
            s=df['overlap_size'] * 10,
            c=df['-log10(padj)'],
            cmap='viridis'
        )
        
        ax.set_yticks(range(len(df)))
        ax.set_yticklabels(df['term_name'])
        ax.set_xlabel('Odds Ratio')
        plt.colorbar(scatter, label='-log₁₀(Adj. P)')
        
        return fig
```

## 4.3 Visualization Summary Table

| Plot Type | Purpose | Static (Matplotlib) | Interactive (Plotly) |
|-----------|---------|---------------------|---------------------|
| **Volcano** | DE overview | ✓ Publication | ✓ Hover for gene names |
| **MA** | Expression vs FC | ✓ | Optional |
| **PCA** | Sample relationships | ✓ | ✓ Sample labels |
| **Heatmap** | Expression patterns | ✓ Clustermap | Limited |
| **Enrichment dot** | Pathway results | ✓ | ✓ |
| **Box/Violin** | Gene expression | ✓ | ✓ |
| **Bar (stacked)** | Cell proportions | ✓ | ✓ |

---

# 5. Interpretation Engine

## 5.1 Automated Interpretation Framework

```python
@dataclass
class Interpretation:
    """Single interpretation or insight."""
    category: str      # 'de', 'enrichment', 'qc', 'deconvolution'
    level: str         # 'critical', 'important', 'informative'
    title: str
    description: str
    evidence: Dict
    recommendations: List[str]

class InterpretationEngine:
    """Generate publication-ready interpretations."""
    
    def interpret_de_results(
        self,
        de_results: pd.DataFrame,
        contrast_name: str,
        padj_threshold: float = 0.05,
        lfc_threshold: float = 1.0
    ) -> List[Interpretation]:
        """Generate DE interpretations."""
        
        interpretations = []
        
        # Count significant genes
        n_genes = len(de_results)
        n_up = ((de_results['padj'] < padj_threshold) & 
                (de_results['log2FoldChange'] > lfc_threshold)).sum()
        n_down = ((de_results['padj'] < padj_threshold) & 
                  (de_results['log2FoldChange'] < -lfc_threshold)).sum()
        
        # Summary interpretation
        interpretations.append(Interpretation(
            category='de',
            level='important',
            title='Differential Expression Summary',
            description=(
                f"Analysis of {n_genes:,} genes revealed {n_up + n_down:,} "
                f"significantly differentially expressed genes in {contrast_name}. "
                f"Of these, {n_up:,} were upregulated and {n_down:,} downregulated."
            ),
            evidence={
                'total_genes': n_genes,
                'upregulated': n_up,
                'downregulated': n_down
            },
            recommendations=[
                "Examine top DE genes for biological relevance",
                "Perform pathway enrichment analysis",
                "Validate key findings with qPCR"
            ]
        ))
        
        return interpretations
    
    def interpret_enrichment(
        self,
        enrichment_results: pd.DataFrame,
        top_n: int = 10
    ) -> List[Interpretation]:
        """Generate enrichment interpretations."""
        
        interpretations = []
        
        if len(enrichment_results) == 0:
            interpretations.append(Interpretation(
                category='enrichment',
                level='important',
                title='No Significant Enrichment',
                description="No significantly enriched pathways identified.",
                evidence={},
                recommendations=["Try different databases", "Use GSEA"]
            ))
            return interpretations
        
        top_terms = enrichment_results.head(top_n)
        
        # Check for immune enrichment
        immune_keywords = ['immune', 'inflammatory', 'cytokine', 'T cell', 'B cell']
        immune_terms = enrichment_results[
            enrichment_results['term_name'].str.contains(
                '|'.join(immune_keywords), case=False, na=False
            )
        ]
        
        if len(immune_terms) > 0:
            interpretations.append(Interpretation(
                category='enrichment',
                level='important',
                title='Immune-Related Enrichment Detected',
                description=(
                    f"Strong enrichment of {len(immune_terms)} immune-related pathways, "
                    f"suggesting significant immune involvement."
                ),
                evidence={'immune_terms': immune_terms['term_name'].tolist()[:5]},
                recommendations=[
                    "Consider cell type deconvolution",
                    "Investigate cytokine/chemokine profiles"
                ]
            ))
        
        return interpretations
    
    def generate_methods_text(self, parameters: Dict) -> str:
        """Generate methods section text."""
        
        methods = """
## Methods

### RNA Sequencing and Data Processing

Reads were quantified using {quant_tool} with a {reference} reference. 
Gene-level counts were aggregated using tximport.

### Differential Expression Analysis

Differential expression analysis was performed using {de_tool} with the 
design formula: {design}. Genes with adjusted p-value < {padj} and 
|log₂ fold change| > {lfc} were considered significant.

### Pathway Enrichment Analysis

{enrichment_method} was performed using gene sets from {databases}. 
P-values were adjusted using the Benjamini-Hochberg method.
""".format(**parameters)
        
        return methods
```

---

# 6. UI Architecture

## 6.1 Streamlit Application Structure

```python
# app.py - Main Streamlit application
import streamlit as st
import pandas as pd

# Page configuration
st.set_page_config(
    page_title="RNA-seq Analysis Platform",
    page_icon="🧬",
    layout="wide"
)

# Sidebar navigation
page = st.sidebar.radio(
    "Navigation",
    ["🏠 Home", "📊 Upload Data", "🔬 QC", "📈 DE Analysis",
     "🔗 Enrichment", "🎯 Clustering", "📝 Report"]
)

# Session state for data persistence
if 'data' not in st.session_state:
    st.session_state.data = {}
if 'results' not in st.session_state:
    st.session_state.results = {}

# Route to pages
if page == "🏠 Home":
    show_home()
elif page == "📊 Upload Data":
    show_upload()
elif page == "📈 DE Analysis":
    show_de_analysis()
# ... etc
```

### DE Analysis Page Example

```python
def show_de_analysis():
    st.header("📈 Differential Expression Analysis")
    
    if st.session_state.data.get('counts') is None:
        st.warning("Please upload data first")
        return
    
    # Configuration sidebar
    col1, col2 = st.columns(2)
    
    with col1:
        design_var = st.selectbox("Comparison variable", 
                                   metadata.columns.tolist())
        numerator = st.selectbox("Treatment group", levels)
        denominator = st.selectbox("Control group", 
                                    [l for l in levels if l != numerator])
    
    with col2:
        padj = st.slider("Adjusted p-value threshold", 0.001, 0.1, 0.05)
        lfc = st.slider("Log2 fold change threshold", 0.0, 3.0, 1.0)
    
    # Run analysis button
    if st.button("🚀 Run Analysis", type="primary"):
        with st.spinner("Running..."):
            result = run_de_analysis(...)
            st.session_state.results['de'] = result
    
    # Display results
    if st.session_state.results.get('de'):
        # Tabs for different views
        tab1, tab2, tab3 = st.tabs(["Volcano", "Table", "Interpretations"])
        
        with tab1:
            fig = plot_volcano(results, interactive=True)
            st.plotly_chart(fig, use_container_width=True)
        
        with tab2:
            st.dataframe(results.sort_values('padj').head(100))
            st.download_button("Download", results.to_csv(), "de_results.csv")
        
        with tab3:
            interpretations = interpret_de_results(results)
            for interp in interpretations:
                st.info(f"**{interp.title}**\n\n{interp.description}")
```

---

# 7. Pipeline Orchestration

## 7.1 Snakemake Workflow

```python
# Snakefile
configfile: "config.yaml"

SAMPLES = config["samples"]

rule all:
    input:
        "results/multiqc/multiqc_report.html",
        "results/de/de_results.csv",
        "results/enrichment/enrichment_results.csv",
        "results/report/analysis_report.html"

rule fastqc:
    input:
        r1 = "data/raw/{sample}_R1.fastq.gz",
        r2 = "data/raw/{sample}_R2.fastq.gz"
    output:
        "results/qc/fastqc/{sample}_R1_fastqc.html"
    shell:
        "fastqc -o results/qc/fastqc {input.r1} {input.r2}"

rule salmon_quant:
    input:
        r1 = "results/trimmed/{sample}_R1.fastq.gz",
        r2 = "results/trimmed/{sample}_R2.fastq.gz",
        index = config["salmon_index"]
    output:
        directory("results/salmon/{sample}")
    threads: 8
    shell:
        """
        salmon quant -i {input.index} -l A \
            -1 {input.r1} -2 {input.r2} \
            -o {output} --gcBias --seqBias
        """

rule differential_expression:
    input:
        counts = "results/counts/gene_counts.csv",
        metadata = config["metadata"]
    output:
        "results/de/de_results.csv"
    script:
        "scripts/run_de.py"
```

## 7.2 Configuration File

```yaml
# config.yaml
samples:
  - sample1
  - sample2
  - sample3

salmon_index: "references/salmon_index"
metadata: "data/sample_metadata.csv"

design_formula: "~ condition"
contrast:
  variable: "condition"
  numerator: "treatment"
  denominator: "control"

de_params:
  padj_threshold: 0.05
  lfc_threshold: 1.0
  shrinkage: true

enrichment_databases:
  - "GO_Biological_Process_2023"
  - "KEGG_2021_Human"
```

---

# 8. Project Structure

```
bulk_rnaseq_platform/
│
├── README.md
├── pyproject.toml
├── requirements.txt
│
├── config/
│   ├── config.yaml
│   └── default_params.yaml
│
├── workflow/
│   ├── Snakefile
│   ├── rules/
│   │   ├── qc.smk
│   │   ├── quantification.smk
│   │   ├── de_analysis.smk
│   │   └── enrichment.smk
│   └── scripts/
│       ├── run_de.py
│       └── generate_report.py
│
├── src/rnaseq_platform/
│   ├── __init__.py
│   ├── qc/
│   │   └── quality_control.py
│   ├── quantification/
│   │   ├── salmon.py
│   │   └── aggregation.py
│   ├── analysis/
│   │   ├── differential_expression.py
│   │   ├── enrichment.py
│   │   ├── clustering.py
│   │   └── deconvolution.py
│   ├── visualization/
│   │   ├── plots.py
│   │   └── interactive.py
│   └── interpretation/
│       └── engine.py
│
├── app/
│   ├── main.py              # Streamlit entry point
│   └── pages/
│       ├── home.py
│       ├── upload.py
│       ├── de_analysis.py
│       └── report.py
│
├── tests/
│   ├── test_de.py
│   └── test_enrichment.py
│
└── examples/
    ├── demo_data/
    └── notebooks/
        └── quickstart.ipynb
```

---

# Summary: Key Capabilities Checklist

## Analysis Capabilities

- [x] **Gene Expression Quantification**: Salmon/STAR + featureCounts
- [x] **Differential Expression**: PyDESeq2/DESeq2 with shrinkage
- [x] **Pathway Enrichment**: ORA + GSEA via gseapy
- [x] **Sample Clustering**: PCA, UMAP, hierarchical
- [x] **Cell Type Deconvolution**: NNLS, EPIC, CIBERSORTx support

## Visualization Capabilities

- [x] Volcano plots (static + interactive)
- [x] MA plots
- [x] PCA plots with metadata coloring
- [x] Clustered heatmaps with annotations
- [x] Enrichment dot plots
- [x] Cell type proportion bar charts
- [x] Sample correlation heatmaps

## Platform Features

- [x] Reproducible Snakemake pipelines
- [x] Interactive Streamlit web UI
- [x] Automated interpretation engine
- [x] Publication-ready figures (300 DPI)
- [x] Methods text generation
- [x] Downloadable results (CSV, PDF)
- [x] Session state persistence

---

*This specification provides a complete blueprint for building a PhD-grade, publication-ready bulk RNA-seq analysis platform.*
