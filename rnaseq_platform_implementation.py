"""
Bulk RNA-Seq Analysis Platform - Complete Python Implementation
================================================================

This module contains all core analysis components for a publication-ready
RNA-seq analysis platform.

Author: Generated for PhD-level research
"""

# ============================================================================
# IMPORTS
# ============================================================================

import numpy as np
import pandas as pd
from dataclasses import dataclass, field
from typing import List, Dict, Optional, Tuple, Set, Union, Any
from pathlib import Path
import subprocess
from scipy import stats
from scipy.optimize import nnls
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import pdist, squareform
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from statsmodels.stats.multitest import multipletests
import warnings

# Visualization imports
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns

# Set publication defaults
plt.rcParams.update({
    'font.family': 'Arial',
    'font.size': 10,
    'axes.labelsize': 12,
    'axes.titlesize': 14,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'axes.spines.top': False,
    'axes.spines.right': False,
})

# ============================================================================
# DATA CLASSES
# ============================================================================

@dataclass
class DEResult:
    """Container for differential expression results."""
    results_table: pd.DataFrame
    normalized_counts: pd.DataFrame
    size_factors: pd.Series
    design_matrix: pd.DataFrame
    contrasts: Dict[str, str]
    dispersions: Optional[pd.Series] = None
    
    def get_significant(
        self,
        padj_threshold: float = 0.05,
        lfc_threshold: float = 1.0
    ) -> pd.DataFrame:
        """Get significantly DE genes."""
        mask = (
            (self.results_table['padj'] < padj_threshold) &
            (self.results_table['log2FoldChange'].abs() > lfc_threshold)
        )
        return self.results_table[mask].sort_values('padj')
    
    def get_upregulated(self, padj: float = 0.05, lfc: float = 1.0) -> pd.DataFrame:
        """Get upregulated genes."""
        mask = (self.results_table['padj'] < padj) & (self.results_table['log2FoldChange'] > lfc)
        return self.results_table[mask].sort_values('padj')
    
    def get_downregulated(self, padj: float = 0.05, lfc: float = 1.0) -> pd.DataFrame:
        """Get downregulated genes."""
        mask = (self.results_table['padj'] < padj) & (self.results_table['log2FoldChange'] < -lfc)
        return self.results_table[mask].sort_values('padj')


@dataclass
class PCAResult:
    """Container for PCA results."""
    embedding: pd.DataFrame
    explained_variance: np.ndarray
    loadings: pd.DataFrame
    n_components: int
    genes_used: List[str]


@dataclass
class EnrichmentResult:
    """Container for enrichment results."""
    results_table: pd.DataFrame
    query_genes: List[str]
    background_size: int
    databases_used: List[str]


@dataclass
class DeconvolutionResult:
    """Container for deconvolution results."""
    proportions: pd.DataFrame
    method: str
    goodness_of_fit: Optional[pd.Series] = None


@dataclass
class Interpretation:
    """Single interpretation or insight."""
    category: str
    level: str  # 'critical', 'important', 'informative'
    title: str
    description: str
    evidence: Dict[str, Any]
    recommendations: List[str] = field(default_factory=list)


# ============================================================================
# MODULE 1: QUALITY CONTROL
# ============================================================================

class QualityControl:
    """
    Quality control module for RNA-seq data.
    
    Wraps FastQC, MultiQC, and fastp for comprehensive QC.
    """
    
    def __init__(self, output_dir: Path, threads: int = 4):
        self.output_dir = Path(output_dir)
        self.threads = threads
        self.output_dir.mkdir(parents=True, exist_ok=True)
    
    def run_fastqc(self, fastq_files: List[Path]) -> Path:
        """
        Run FastQC on FASTQ files.
        
        Parameters
        ----------
        fastq_files : List[Path]
            List of FASTQ files (can be gzipped)
        
        Returns
        -------
        Path
            Directory containing FastQC outputs
        """
        fastqc_dir = self.output_dir / "fastqc"
        fastqc_dir.mkdir(exist_ok=True)
        
        cmd = [
            "fastqc",
            "--outdir", str(fastqc_dir),
            "--threads", str(self.threads),
            "--quiet"
        ] + [str(f) for f in fastq_files]
        
        subprocess.run(cmd, check=True)
        return fastqc_dir
    
    def run_multiqc(self, input_dir: Path, title: str = "RNA-seq QC") -> Path:
        """Aggregate QC reports with MultiQC."""
        multiqc_dir = self.output_dir / "multiqc"
        
        cmd = [
            "multiqc", str(input_dir),
            "--outdir", str(multiqc_dir),
            "--title", title,
            "--force"
        ]
        
        subprocess.run(cmd, check=True)
        return multiqc_dir / "multiqc_report.html"
    
    def trim_reads(
        self,
        r1: Path,
        r2: Optional[Path],
        sample_name: str,
        min_length: int = 36,
        quality: int = 20
    ) -> Tuple[Path, Path, Path]:
        """
        Trim reads with fastp.
        
        Returns
        -------
        Tuple[Path, Path, Path]
            Trimmed R1, trimmed R2 (or None), JSON report
        """
        trim_dir = self.output_dir / "trimmed"
        trim_dir.mkdir(exist_ok=True)
        
        r1_out = trim_dir / f"{sample_name}_R1_trimmed.fastq.gz"
        json_out = trim_dir / f"{sample_name}_fastp.json"
        html_out = trim_dir / f"{sample_name}_fastp.html"
        
        cmd = [
            "fastp",
            "--in1", str(r1),
            "--out1", str(r1_out),
            "--json", str(json_out),
            "--html", str(html_out),
            "--thread", str(self.threads),
            "--length_required", str(min_length),
            "--qualified_quality_phred", str(quality),
            "--detect_adapter_for_pe",
            "--trim_poly_g"
        ]
        
        if r2:
            r2_out = trim_dir / f"{sample_name}_R2_trimmed.fastq.gz"
            cmd.extend(["--in2", str(r2), "--out2", str(r2_out)])
        else:
            r2_out = None
        
        subprocess.run(cmd, check=True)
        return r1_out, r2_out, json_out


# ============================================================================
# MODULE 2: QUANTIFICATION
# ============================================================================

class SalmonQuantifier:
    """
    Gene expression quantification using Salmon.
    
    Salmon provides fast, accurate transcript-level quantification
    without full alignment.
    """
    
    def __init__(self, index_path: Path, output_dir: Path, threads: int = 8):
        self.index_path = Path(index_path)
        self.output_dir = Path(output_dir)
        self.threads = threads
        self.output_dir.mkdir(parents=True, exist_ok=True)
    
    def quantify_sample(
        self,
        r1: Path,
        r2: Optional[Path],
        sample_name: str,
        library_type: str = "A"
    ) -> Path:
        """
        Quantify a single sample.
        
        Parameters
        ----------
        r1 : Path
            R1 FASTQ file
        r2 : Optional[Path]
            R2 FASTQ file (None for single-end)
        sample_name : str
            Sample identifier
        library_type : str
            Library type (A=auto-detect)
        
        Returns
        -------
        Path
            Output directory with quant.sf
        """
        sample_out = self.output_dir / sample_name
        
        cmd = [
            "salmon", "quant",
            "-i", str(self.index_path),
            "-l", library_type,
            "-o", str(sample_out),
            "-p", str(self.threads),
            "--validateMappings",
            "--gcBias",
            "--seqBias"
        ]
        
        if r2:
            cmd.extend(["-1", str(r1), "-2", str(r2)])
        else:
            cmd.extend(["-r", str(r1)])
        
        subprocess.run(cmd, check=True)
        return sample_out
    
    def aggregate_counts(
        self,
        sample_dirs: List[Path],
        tx2gene: pd.DataFrame
    ) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """
        Aggregate transcript counts to gene level.
        
        Parameters
        ----------
        sample_dirs : List[Path]
            Salmon output directories
        tx2gene : pd.DataFrame
            Transcript to gene mapping (transcript_id, gene_id, gene_name)
        
        Returns
        -------
        Tuple[pd.DataFrame, pd.DataFrame]
            Gene-level counts and TPM matrices
        """
        counts_list = []
        tpm_list = []
        
        for sample_dir in sample_dirs:
            quant_file = sample_dir / "quant.sf"
            sample_name = sample_dir.name
            
            df = pd.read_csv(quant_file, sep='\t')
            
            # Merge with gene mapping
            df = df.merge(tx2gene, left_on='Name', right_on='transcript_id')
            
            # Aggregate to gene level
            gene_counts = df.groupby('gene_id')['NumReads'].sum()
            gene_tpm = df.groupby('gene_id')['TPM'].sum()
            
            gene_counts.name = sample_name
            gene_tpm.name = sample_name
            
            counts_list.append(gene_counts)
            tpm_list.append(gene_tpm)
        
        counts_matrix = pd.concat(counts_list, axis=1)
        tpm_matrix = pd.concat(tpm_list, axis=1)
        
        return counts_matrix, tpm_matrix


# ============================================================================
# MODULE 3: DIFFERENTIAL EXPRESSION
# ============================================================================

class DifferentialExpression:
    """
    Differential expression analysis using PyDESeq2.
    
    Implements the DESeq2 statistical framework:
    - Negative binomial distribution
    - Empirical Bayes shrinkage
    - Wald test for significance
    """
    
    def __init__(self, min_counts: int = 10, min_samples: int = 3):
        self.min_counts = min_counts
        self.min_samples = min_samples
    
    def run_analysis(
        self,
        counts: pd.DataFrame,
        metadata: pd.DataFrame,
        design_formula: str,
        contrast: Tuple[str, str, str],
        shrinkage: bool = True
    ) -> DEResult:
        """
        Run differential expression analysis.
        
        Parameters
        ----------
        counts : pd.DataFrame
            Raw count matrix (genes x samples)
        metadata : pd.DataFrame
            Sample metadata with design variables
        design_formula : str
            Design formula (e.g., "~ condition")
        contrast : Tuple[str, str, str]
            (variable, numerator, denominator)
        shrinkage : bool
            Apply LFC shrinkage
        
        Returns
        -------
        DEResult
            Differential expression results
        """
        try:
            from pydeseq2.dds import DeseqDataSet
            from pydeseq2.ds import DeseqStats
        except ImportError:
            raise ImportError("PyDESeq2 required. Install with: pip install pydeseq2")
        
        # Validate inputs
        self._validate_inputs(counts, metadata, contrast)
        
        # Filter low-count genes
        counts_filtered = self._filter_genes(counts)
        
        # Align samples
        common_samples = counts_filtered.columns.intersection(metadata.index)
        counts_aligned = counts_filtered[common_samples]
        metadata_aligned = metadata.loc[common_samples]
        
        # Create DESeq dataset
        design_factor = design_formula.replace("~", "").strip()
        
        dds = DeseqDataSet(
            counts=counts_aligned.T,  # PyDESeq2 expects samples x genes
            metadata=metadata_aligned,
            design_factors=design_factor
        )
        
        # Run DESeq2
        dds.deseq2()
        
        # Get statistics for contrast
        stat_res = DeseqStats(dds, contrast=list(contrast))
        stat_res.summary()
        
        # Apply shrinkage if requested
        if shrinkage:
            coef_name = f"{contrast[0]}_{contrast[1]}_vs_{contrast[2]}"
            try:
                stat_res.lfc_shrink(coeff=coef_name)
            except Exception:
                warnings.warn("LFC shrinkage failed, using unshrunk estimates")
        
        # Build results
        results_df = stat_res.results_df.copy()
        
        # Get normalized counts
        normalized = pd.DataFrame(
            dds.layers['normed_counts'],
            index=dds.obs_names,
            columns=dds.var_names
        ).T
        
        return DEResult(
            results_table=results_df,
            normalized_counts=normalized,
            size_factors=pd.Series(dds.obsm['size_factors'], index=dds.obs_names),
            design_matrix=metadata_aligned,
            contrasts={f"{contrast[1]}_vs_{contrast[2]}": str(contrast)},
            dispersions=pd.Series(dds.varm['dispersions'], index=dds.var_names)
        )
    
    def _validate_inputs(
        self,
        counts: pd.DataFrame,
        metadata: pd.DataFrame,
        contrast: Tuple[str, str, str]
    ):
        """Validate input data."""
        # Check contrast variable exists
        if contrast[0] not in metadata.columns:
            raise ValueError(f"Contrast variable '{contrast[0]}' not in metadata")
        
        # Check levels exist
        levels = metadata[contrast[0]].unique()
        for level in contrast[1:]:
            if level not in levels:
                raise ValueError(f"Level '{level}' not found in '{contrast[0]}'")
    
    def _filter_genes(self, counts: pd.DataFrame) -> pd.DataFrame:
        """Filter lowly expressed genes."""
        # Total counts filter
        total = counts.sum(axis=1)
        mask_total = total >= self.min_counts
        
        # Detection filter
        detected = (counts > 0).sum(axis=1)
        mask_detected = detected >= self.min_samples
        
        filtered = counts[mask_total & mask_detected]
        
        n_removed = len(counts) - len(filtered)
        print(f"Filtered {n_removed} low-count genes ({len(filtered)} remaining)")
        
        return filtered


# ============================================================================
# MODULE 4: PATHWAY ENRICHMENT
# ============================================================================

class PathwayEnrichment:
    """
    Pathway and functional enrichment analysis.
    
    Supports:
    - Over-representation analysis (ORA)
    - Gene Set Enrichment Analysis (GSEA)
    """
    
    ENRICHR_DATABASES = [
        'GO_Biological_Process_2023',
        'GO_Molecular_Function_2023',
        'GO_Cellular_Component_2023',
        'KEGG_2021_Human',
        'Reactome_2022',
        'WikiPathway_2023_Human',
        'MSigDB_Hallmark_2020'
    ]
    
    def __init__(self):
        self.gene_sets: Dict[str, Dict[str, Set[str]]] = {}
    
    def run_ora(
        self,
        query_genes: List[str],
        databases: List[str],
        background_genes: Optional[Set[str]] = None,
        min_size: int = 10,
        max_size: int = 500,
        padj_threshold: float = 0.05
    ) -> pd.DataFrame:
        """
        Run over-representation analysis.
        
        Parameters
        ----------
        query_genes : List[str]
            Genes to test (e.g., DE genes)
        databases : List[str]
            Gene set databases to use
        background_genes : Optional[Set[str]]
            Background gene universe
        min_size : int
            Minimum gene set size
        max_size : int
            Maximum gene set size
        padj_threshold : float
            Adjusted p-value threshold
        
        Returns
        -------
        pd.DataFrame
            Enrichment results
        """
        try:
            import gseapy as gp
        except ImportError:
            raise ImportError("gseapy required. Install with: pip install gseapy")
        
        query_set = set(query_genes)
        results = []
        
        for db in databases:
            # Use gseapy's enrichr
            try:
                enr = gp.enrichr(
                    gene_list=list(query_genes),
                    gene_sets=db,
                    organism='human',
                    cutoff=0.5  # Initial cutoff, filter later
                )
                
                if enr.results is not None and len(enr.results) > 0:
                    df = enr.results.copy()
                    df['source'] = db
                    results.append(df)
            except Exception as e:
                warnings.warn(f"Error with {db}: {e}")
                continue
        
        if not results:
            return pd.DataFrame()
        
        combined = pd.concat(results, ignore_index=True)
        
        # Standardize column names
        combined = combined.rename(columns={
            'Term': 'term_name',
            'P-value': 'p_value',
            'Adjusted P-value': 'adjusted_p_value',
            'Odds Ratio': 'odds_ratio',
            'Combined Score': 'combined_score',
            'Genes': 'overlap_genes'
        })
        
        # Filter
        combined = combined[combined['adjusted_p_value'] < padj_threshold]
        combined = combined.sort_values('combined_score', ascending=False)
        
        return combined
    
    def run_gsea(
        self,
        gene_ranks: pd.Series,
        database: str,
        min_size: int = 15,
        max_size: int = 500,
        permutations: int = 1000
    ) -> pd.DataFrame:
        """
        Run Gene Set Enrichment Analysis.
        
        Parameters
        ----------
        gene_ranks : pd.Series
            Genes ranked by expression change (e.g., -log10(p) * sign(lfc))
        database : str
            Gene set database
        
        Returns
        -------
        pd.DataFrame
            GSEA results
        """
        try:
            import gseapy as gp
        except ImportError:
            raise ImportError("gseapy required")
        
        pre_res = gp.prerank(
            rnk=gene_ranks,
            gene_sets=database,
            min_size=min_size,
            max_size=max_size,
            permutation_num=permutations,
            seed=42,
            verbose=False
        )
        
        results = pre_res.res2d.copy()
        results = results.rename(columns={
            'Term': 'term_name',
            'NES': 'normalized_enrichment_score',
            'NOM p-val': 'p_value',
            'FDR q-val': 'fdr_q_value'
        })
        
        return results
    
    @staticmethod
    def create_ranking_metric(de_results: pd.DataFrame) -> pd.Series:
        """Create gene ranking from DE results."""
        pval = de_results['pvalue'].clip(lower=1e-300)
        ranking = -np.log10(pval) * np.sign(de_results['log2FoldChange'])
        return ranking.dropna().sort_values(ascending=False)


# ============================================================================
# MODULE 5: SAMPLE CLUSTERING
# ============================================================================

class SampleClustering:
    """
    Sample relationship analysis including PCA, UMAP, and clustering.
    """
    
    def __init__(self, counts: pd.DataFrame, metadata: pd.DataFrame):
        """
        Initialize with normalized counts.
        
        Parameters
        ----------
        counts : pd.DataFrame
            Normalized counts (genes x samples)
        metadata : pd.DataFrame
            Sample metadata
        """
        self.counts = counts
        self.metadata = metadata.loc[counts.columns]
        
        # Variance stabilizing transformation (simplified)
        self.vst = self._variance_stabilize(counts)
    
    def _variance_stabilize(self, counts: pd.DataFrame) -> pd.DataFrame:
        """Apply variance stabilizing transformation."""
        log_counts = np.log2(counts + 1)
        scaler = StandardScaler()
        scaled = scaler.fit_transform(log_counts.T).T
        return pd.DataFrame(scaled, index=counts.index, columns=counts.columns)
    
    def run_pca(
        self,
        n_components: int = 50,
        n_top_genes: int = 2000
    ) -> PCAResult:
        """
        Run Principal Component Analysis.
        
        Parameters
        ----------
        n_components : int
            Number of components
        n_top_genes : int
            Number of highly variable genes to use
        
        Returns
        -------
        PCAResult
            PCA results
        """
        # Select top variable genes
        gene_vars = self.vst.var(axis=1)
        top_genes = gene_vars.nlargest(n_top_genes).index.tolist()
        data = self.vst.loc[top_genes].T
        
        # Fit PCA
        n_components = min(n_components, min(data.shape) - 1)
        pca = PCA(n_components=n_components)
        embedding = pca.fit_transform(data)
        
        # Create results
        pc_names = [f"PC{i+1}" for i in range(n_components)]
        
        return PCAResult(
            embedding=pd.DataFrame(embedding, index=data.index, columns=pc_names),
            explained_variance=pca.explained_variance_ratio_,
            loadings=pd.DataFrame(pca.components_.T, index=data.columns, columns=pc_names),
            n_components=n_components,
            genes_used=top_genes
        )
    
    def compute_distances(
        self,
        method: str = 'correlation',
        n_top_genes: int = 2000
    ) -> pd.DataFrame:
        """Compute sample distance matrix."""
        gene_vars = self.vst.var(axis=1)
        top_genes = gene_vars.nlargest(n_top_genes).index
        data = self.vst.loc[top_genes].T
        
        if method == 'correlation':
            corr = data.T.corr()
            return 1 - corr
        else:
            dist = pdist(data, metric=method)
            return pd.DataFrame(
                squareform(dist),
                index=data.index,
                columns=data.index
            )
    
    def detect_outliers(
        self,
        pca_result: PCAResult,
        threshold_sd: float = 3.0
    ) -> List[str]:
        """Detect outlier samples."""
        pc_data = pca_result.embedding[['PC1', 'PC2']]
        
        centroid = pc_data.mean()
        distances = np.sqrt(((pc_data - centroid) ** 2).sum(axis=1))
        
        threshold = distances.mean() + threshold_sd * distances.std()
        return distances[distances > threshold].index.tolist()
    
    def assess_batch_effects(
        self,
        pca_result: PCAResult,
        batch_variable: str
    ) -> Dict[str, float]:
        """Assess variance explained by batch."""
        results = {}
        batch_labels = self.metadata[batch_variable]
        
        for pc in ['PC1', 'PC2', 'PC3']:
            if pc not in pca_result.embedding.columns:
                continue
            
            pc_values = pca_result.embedding[pc]
            
            # Calculate variance explained by batch
            groups = [pc_values[batch_labels == b].values for b in batch_labels.unique()]
            ss_between = sum(len(g) * (g.mean() - pc_values.mean())**2 for g in groups)
            ss_total = ((pc_values - pc_values.mean())**2).sum()
            
            results[pc] = ss_between / ss_total
        
        return results


# ============================================================================
# MODULE 6: CELL TYPE DECONVOLUTION
# ============================================================================

class CellTypeDeconvolution:
    """
    Cell type deconvolution from bulk RNA-seq.
    
    Estimates cell type proportions using signature-based methods.
    """
    
    def __init__(self):
        self.signature = None
    
    def load_signature(self, signature_df: pd.DataFrame):
        """Load cell type signature matrix."""
        self.signature = signature_df
    
    def run_nnls(
        self,
        expression: pd.DataFrame,
        signature: Optional[pd.DataFrame] = None
    ) -> DeconvolutionResult:
        """
        Run NNLS deconvolution.
        
        Parameters
        ----------
        expression : pd.DataFrame
            Normalized expression (TPM preferred), genes x samples
        signature : Optional[pd.DataFrame]
            Signature matrix (genes x cell types)
        
        Returns
        -------
        DeconvolutionResult
            Cell type proportions
        """
        if signature is None:
            signature = self.signature
        
        if signature is None:
            raise ValueError("No signature matrix provided")
        
        # Find common genes
        common = expression.index.intersection(signature.index)
        
        if len(common) < 50:
            raise ValueError(f"Only {len(common)} genes overlap")
        
        expr = expression.loc[common]
        sig = signature.loc[common]
        
        # NNLS for each sample
        proportions = []
        residuals = []
        
        for sample in expr.columns:
            coef, residual = nnls(sig.values, expr[sample].values)
            
            # Normalize to sum to 1
            coef_norm = coef / coef.sum() if coef.sum() > 0 else coef
            
            proportions.append(coef_norm)
            residuals.append(residual)
        
        return DeconvolutionResult(
            proportions=pd.DataFrame(
                proportions,
                index=expr.columns,
                columns=sig.columns
            ),
            method='NNLS',
            goodness_of_fit=pd.Series(residuals, index=expr.columns)
        )


# ============================================================================
# MODULE 7: VISUALIZATION
# ============================================================================

class RNASeqVisualizer:
    """
    Publication-quality visualization suite.
    """
    
    def __init__(self, style: str = 'ticks', context: str = 'paper'):
        sns.set_style(style)
        sns.set_context(context, font_scale=1.2)
    
    def plot_volcano(
        self,
        de_results: pd.DataFrame,
        lfc_threshold: float = 1.0,
        pval_threshold: float = 0.05,
        label_top_n: int = 10,
        figsize: Tuple[float, float] = (8, 6)
    ) -> plt.Figure:
        """Create volcano plot."""
        
        df = de_results.copy()
        df['-log10(padj)'] = -np.log10(df['padj'].clip(lower=1e-300))
        
        # Classify
        df['status'] = 'NS'
        df.loc[(df['log2FoldChange'] > lfc_threshold) & 
               (df['padj'] < pval_threshold), 'status'] = 'Up'
        df.loc[(df['log2FoldChange'] < -lfc_threshold) & 
               (df['padj'] < pval_threshold), 'status'] = 'Down'
        
        fig, ax = plt.subplots(figsize=figsize)
        
        colors = {'Up': '#e74c3c', 'Down': '#3498db', 'NS': '#95a5a6'}
        
        for status in ['NS', 'Up', 'Down']:
            mask = df['status'] == status
            n = mask.sum()
            ax.scatter(
                df.loc[mask, 'log2FoldChange'],
                df.loc[mask, '-log10(padj)'],
                c=colors[status],
                alpha=0.6,
                s=20,
                label=f'{status} ({n})'
            )
        
        # Threshold lines
        ax.axvline(lfc_threshold, ls='--', c='gray', alpha=0.5)
        ax.axvline(-lfc_threshold, ls='--', c='gray', alpha=0.5)
        ax.axhline(-np.log10(pval_threshold), ls='--', c='gray', alpha=0.5)
        
        # Label top genes
        if label_top_n > 0:
            top = df[df['status'] != 'NS'].nlargest(label_top_n, '-log10(padj)')
            for gene, row in top.iterrows():
                ax.annotate(gene, (row['log2FoldChange'], row['-log10(padj)']),
                           fontsize=8, alpha=0.8)
        
        ax.set_xlabel('log₂(Fold Change)')
        ax.set_ylabel('-log₁₀(Adjusted P-value)')
        ax.legend(loc='upper right')
        
        plt.tight_layout()
        return fig
    
    def plot_pca(
        self,
        pca_result: PCAResult,
        metadata: pd.DataFrame,
        color_by: str,
        pc_x: int = 1,
        pc_y: int = 2,
        figsize: Tuple[float, float] = (8, 6)
    ) -> plt.Figure:
        """Create PCA plot."""
        
        embedding = pca_result.embedding
        var = pca_result.explained_variance
        
        fig, ax = plt.subplots(figsize=figsize)
        
        groups = metadata[color_by].unique()
        colors = sns.color_palette('Set2', len(groups))
        
        for group, color in zip(groups, colors):
            mask = metadata[color_by] == group
            samples = metadata[mask].index
            samples = [s for s in samples if s in embedding.index]
            
            ax.scatter(
                embedding.loc[samples, f'PC{pc_x}'],
                embedding.loc[samples, f'PC{pc_y}'],
                c=[color],
                label=group,
                s=80,
                alpha=0.7
            )
        
        ax.set_xlabel(f'PC{pc_x} ({var[pc_x-1]*100:.1f}%)')
        ax.set_ylabel(f'PC{pc_y} ({var[pc_y-1]*100:.1f}%)')
        ax.legend(title=color_by)
        
        plt.tight_layout()
        return fig
    
    def plot_heatmap(
        self,
        expression: pd.DataFrame,
        genes: List[str],
        metadata: pd.DataFrame,
        annotation_cols: List[str],
        scale: str = 'row',
        figsize: Tuple[float, float] = (12, 10)
    ) -> plt.Figure:
        """Create clustered heatmap."""
        
        data = expression.loc[genes, metadata.index]
        
        if scale == 'row':
            data = data.sub(data.mean(axis=1), axis=0).div(data.std(axis=1), axis=0)
        
        # Annotation colors
        col_colors = []
        for col in annotation_cols:
            values = metadata[col]
            palette = sns.color_palette('Set2', values.nunique())
            color_map = dict(zip(values.unique(), palette))
            col_colors.append([color_map[v] for v in values])
        
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
    
    def plot_enrichment_dotplot(
        self,
        enrichment_results: pd.DataFrame,
        top_n: int = 20,
        figsize: Tuple[float, float] = (8, None)
    ) -> plt.Figure:
        """Create enrichment dot plot."""
        
        df = enrichment_results.head(top_n).copy()
        df['-log10(padj)'] = -np.log10(df['adjusted_p_value'].clip(1e-300))
        
        height = max(6, top_n * 0.3)
        fig, ax = plt.subplots(figsize=(figsize[0], height))
        
        scatter = ax.scatter(
            df['odds_ratio'] if 'odds_ratio' in df.columns else df['combined_score'],
            range(len(df)),
            s=100,
            c=df['-log10(padj)'],
            cmap='viridis'
        )
        
        ax.set_yticks(range(len(df)))
        ax.set_yticklabels(df['term_name'])
        ax.set_xlabel('Odds Ratio' if 'odds_ratio' in df.columns else 'Combined Score')
        
        plt.colorbar(scatter, label='-log₁₀(Adj. P)')
        plt.tight_layout()
        
        return fig
    
    def plot_cell_proportions(
        self,
        deconv_result: DeconvolutionResult,
        figsize: Tuple[float, float] = (12, 6)
    ) -> plt.Figure:
        """Create stacked bar chart of cell proportions."""
        
        fig, ax = plt.subplots(figsize=figsize)
        
        deconv_result.proportions.plot(
            kind='bar',
            stacked=True,
            ax=ax,
            colormap='tab20'
        )
        
        ax.set_xlabel('Sample')
        ax.set_ylabel('Proportion')
        ax.legend(title='Cell Type', bbox_to_anchor=(1.05, 1))
        
        plt.tight_layout()
        return fig


# ============================================================================
# MODULE 8: INTERPRETATION ENGINE
# ============================================================================

class InterpretationEngine:
    """
    Automated interpretation of RNA-seq results.
    """
    
    def __init__(self):
        self.interpretations: List[Interpretation] = []
    
    def interpret_de_results(
        self,
        de_results: pd.DataFrame,
        contrast_name: str,
        padj_threshold: float = 0.05,
        lfc_threshold: float = 1.0
    ) -> List[Interpretation]:
        """Generate DE interpretations."""
        
        n_genes = len(de_results)
        n_up = ((de_results['padj'] < padj_threshold) & 
                (de_results['log2FoldChange'] > lfc_threshold)).sum()
        n_down = ((de_results['padj'] < padj_threshold) & 
                  (de_results['log2FoldChange'] < -lfc_threshold)).sum()
        
        interpretations = []
        
        # Summary
        interpretations.append(Interpretation(
            category='de',
            level='important',
            title='Differential Expression Summary',
            description=(
                f"Analysis of {n_genes:,} genes revealed {n_up + n_down:,} "
                f"significantly differentially expressed genes in {contrast_name}. "
                f"{n_up:,} were upregulated and {n_down:,} downregulated "
                f"(|log₂FC| > {lfc_threshold}, padj < {padj_threshold})."
            ),
            evidence={
                'total_genes': n_genes,
                'upregulated': n_up,
                'downregulated': n_down,
                'thresholds': {'padj': padj_threshold, 'lfc': lfc_threshold}
            },
            recommendations=[
                "Validate top DE genes with qPCR",
                "Perform pathway enrichment analysis",
                "Examine known markers for your system"
            ]
        ))
        
        # Effect size check
        median_lfc = de_results.loc[
            de_results['padj'] < padj_threshold, 'log2FoldChange'
        ].abs().median()
        
        if median_lfc < 0.5:
            interpretations.append(Interpretation(
                category='de',
                level='informative',
                title='Modest Effect Sizes',
                description=(
                    f"Median |log₂FC| among significant genes is {median_lfc:.2f}, "
                    f"indicating relatively modest effect sizes. Consider GSEA "
                    f"to detect coordinated pathway changes."
                ),
                evidence={'median_abs_lfc': median_lfc},
                recommendations=["Run GSEA with full ranked gene list"]
            ))
        
        self.interpretations.extend(interpretations)
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
                description="No significantly enriched pathways were identified.",
                evidence={},
                recommendations=["Try GSEA", "Use different databases"]
            ))
            return interpretations
        
        # Check for immune terms
        immune_kw = ['immune', 'inflammatory', 'cytokine', 'T cell', 'B cell']
        immune_terms = enrichment_results[
            enrichment_results['term_name'].str.contains(
                '|'.join(immune_kw), case=False, na=False
            )
        ]
        
        if len(immune_terms) > 0:
            interpretations.append(Interpretation(
                category='enrichment',
                level='important',
                title='Immune-Related Enrichment',
                description=(
                    f"Detected {len(immune_terms)} immune-related pathways, "
                    f"suggesting significant immune involvement."
                ),
                evidence={'immune_terms': immune_terms['term_name'].tolist()[:5]},
                recommendations=[
                    "Consider cell type deconvolution",
                    "Investigate cytokine profiles"
                ]
            ))
        
        self.interpretations.extend(interpretations)
        return interpretations
    
    def generate_methods_text(self, params: Dict) -> str:
        """Generate methods section."""
        
        return f"""
## Methods

### RNA Sequencing Analysis

Gene expression quantification was performed using {params.get('quantification', 'Salmon')}.
Differential expression analysis was conducted using {params.get('de_tool', 'PyDESeq2')} 
with the design formula: {params.get('design', '~ condition')}.
Genes with adjusted p-value < {params.get('padj', 0.05)} and 
|log₂ fold change| > {params.get('lfc', 1.0)} were considered significant.

Pathway enrichment analysis was performed using {params.get('enrichment_method', 'over-representation analysis')}
with gene sets from {', '.join(params.get('databases', ['GO', 'KEGG']))}.
"""


# ============================================================================
# MAIN PIPELINE
# ============================================================================

class RNASeqPipeline:
    """
    Main pipeline orchestrator.
    """
    
    def __init__(self, output_dir: Path):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        self.results = {}
    
    def run_full_analysis(
        self,
        counts: pd.DataFrame,
        metadata: pd.DataFrame,
        design_formula: str,
        contrast: Tuple[str, str, str],
        enrichment_databases: List[str] = None
    ) -> Dict:
        """
        Run complete RNA-seq analysis.
        
        Parameters
        ----------
        counts : pd.DataFrame
            Raw count matrix
        metadata : pd.DataFrame
            Sample metadata
        design_formula : str
            Design formula
        contrast : Tuple[str, str, str]
            Contrast specification
        enrichment_databases : List[str]
            Databases for enrichment
        
        Returns
        -------
        Dict
            All analysis results
        """
        if enrichment_databases is None:
            enrichment_databases = ['GO_Biological_Process_2023', 'KEGG_2021_Human']
        
        print("=" * 60)
        print("BULK RNA-SEQ ANALYSIS PIPELINE")
        print("=" * 60)
        
        # 1. Sample Clustering
        print("\n[1/4] Running sample clustering...")
        clustering = SampleClustering(counts, metadata)
        pca_result = clustering.run_pca()
        outliers = clustering.detect_outliers(pca_result)
        
        self.results['pca'] = pca_result
        self.results['outliers'] = outliers
        
        if outliers:
            print(f"    ⚠ Potential outliers: {outliers}")
        else:
            print("    ✓ No outliers detected")
        
        # 2. Differential Expression
        print("\n[2/4] Running differential expression...")
        de = DifferentialExpression()
        de_result = de.run_analysis(counts, metadata, design_formula, contrast)
        
        n_sig = de_result.get_significant().shape[0]
        print(f"    ✓ Found {n_sig} significant genes")
        
        self.results['de'] = de_result
        
        # 3. Enrichment Analysis
        print("\n[3/4] Running pathway enrichment...")
        enrichment = PathwayEnrichment()
        
        sig_genes = de_result.get_significant().index.tolist()
        
        if len(sig_genes) > 0:
            enrichment_results = enrichment.run_ora(sig_genes, enrichment_databases)
            self.results['enrichment'] = enrichment_results
            print(f"    ✓ Found {len(enrichment_results)} enriched terms")
        else:
            print("    ⚠ No significant genes for enrichment")
        
        # 4. Interpretation
        print("\n[4/4] Generating interpretations...")
        interpreter = InterpretationEngine()
        interpreter.interpret_de_results(
            de_result.results_table,
            f"{contrast[1]} vs {contrast[2]}"
        )
        
        if 'enrichment' in self.results:
            interpreter.interpret_enrichment(self.results['enrichment'])
        
        self.results['interpretations'] = interpreter.interpretations
        
        print(f"    ✓ Generated {len(interpreter.interpretations)} interpretations")
        
        print("\n" + "=" * 60)
        print("ANALYSIS COMPLETE")
        print("=" * 60)
        
        return self.results
    
    def generate_visualizations(self, metadata: pd.DataFrame, color_by: str) -> Dict[str, plt.Figure]:
        """Generate all visualizations."""
        
        viz = RNASeqVisualizer()
        figures = {}
        
        # Volcano
        if 'de' in self.results:
            figures['volcano'] = viz.plot_volcano(
                self.results['de'].results_table
            )
        
        # PCA
        if 'pca' in self.results:
            figures['pca'] = viz.plot_pca(
                self.results['pca'],
                metadata,
                color_by
            )
        
        # Enrichment
        if 'enrichment' in self.results and len(self.results['enrichment']) > 0:
            figures['enrichment'] = viz.plot_enrichment_dotplot(
                self.results['enrichment']
            )
        
        return figures
    
    def save_results(self):
        """Save all results to files."""
        
        if 'de' in self.results:
            self.results['de'].results_table.to_csv(
                self.output_dir / 'de_results.csv'
            )
            self.results['de'].normalized_counts.to_csv(
                self.output_dir / 'normalized_counts.csv'
            )
        
        if 'enrichment' in self.results:
            self.results['enrichment'].to_csv(
                self.output_dir / 'enrichment_results.csv',
                index=False
            )
        
        if 'pca' in self.results:
            self.results['pca'].embedding.to_csv(
                self.output_dir / 'pca_coordinates.csv'
            )
        
        print(f"Results saved to {self.output_dir}")


# ============================================================================
# EXAMPLE USAGE
# ============================================================================

if __name__ == "__main__":
    print("""
    ╔═══════════════════════════════════════════════════════════════╗
    ║         BULK RNA-SEQ ANALYSIS PLATFORM                        ║
    ║         Publication-Ready Analysis Pipeline                    ║
    ╚═══════════════════════════════════════════════════════════════╝
    
    Example usage:
    
    ```python
    from rnaseq_platform import RNASeqPipeline
    import pandas as pd
    
    # Load data
    counts = pd.read_csv('counts.csv', index_col=0)
    metadata = pd.read_csv('metadata.csv', index_col=0)
    
    # Run pipeline
    pipeline = RNASeqPipeline(output_dir='results')
    results = pipeline.run_full_analysis(
        counts=counts,
        metadata=metadata,
        design_formula='~ condition',
        contrast=('condition', 'treatment', 'control')
    )
    
    # Generate figures
    figures = pipeline.generate_visualizations(metadata, color_by='condition')
    
    # Save results
    pipeline.save_results()
    ```
    """)
