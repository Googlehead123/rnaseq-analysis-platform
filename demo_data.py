"""
Demo dataset generator for RNA-seq analysis platform.

Generates realistic RNA-seq data for dermatology/cosmetics experiments
with built-in differential expression patterns.
"""

from typing import Tuple
import pandas as pd
import numpy as np


def load_demo_dataset() -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Generate realistic RNA-seq demo dataset for dermatology research.

    Returns:
        Tuple of (counts_df, metadata_df):
        - counts_df: samples × genes (index=sample names, columns=gene names)
          Integer counts from realistic Poisson distribution
        - metadata_df: index=sample names, column "condition" with "Control" or "Treatment"

    Dataset characteristics:
    - 6 samples: Control_Rep1-3, Treatment_Rep1-3
    - ~200 genes including dermatology panel genes
    - Built-in differential expression:
      * Upregulated in Treatment: VEGFA, COL1A1, FN1, MMP1 (wound healing)
      * Downregulated in Treatment: FLG, LOR, CLDN1 (barrier genes)
    - Realistic count distributions (Poisson with log-normal base means)
    - Reproducible with np.random.seed(42)
    """
    np.random.seed(42)

    # Sample names and conditions
    sample_names = [
        "Control_Rep1",
        "Control_Rep2",
        "Control_Rep3",
        "Treatment_Rep1",
        "Treatment_Rep2",
        "Treatment_Rep3",
    ]

    # Dermatology gene panels (from config/gene_panels.yaml)
    dermatology_genes = {
        "Anti-aging": [
            "COL1A1",
            "COL3A1",
            "ELN",
            "FBN1",
            "MMP1",
            "MMP2",
            "MMP9",
            "TIMP1",
        ],
        "Skin Barrier": [
            "FLG",
            "LOR",
            "IVL",
            "CLDN1",
            "CLDN4",
            "TJP1",
            "AQP3",
            "HAS2",
            "HAS3",
        ],
        "Anti-inflammation": [
            "IL1A",
            "IL1B",
            "IL6",
            "IL8",
            "TNFA",
            "CXCL1",
            "CCL2",
            "PTGS2",
        ],
        "Whitening/Melanogenesis": [
            "MITF",
            "TYR",
            "TYRP1",
            "DCT",
            "MC1R",
            "PMEL",
            "OCA2",
        ],
        "Sebum Regulation": ["SREBF1", "PPARG", "FASN", "SCD", "DGAT2", "ELOVL6"],
        "Wound Healing": [
            "VEGFA",
            "FGF2",
            "TGFB1",
            "EGFR",
            "KRT6A",
            "KRT16",
            "KRT17",
            "MMP3",
        ],
    }

    # Flatten all dermatology genes
    all_dermatology_genes = []
    for genes in dermatology_genes.values():
        all_dermatology_genes.extend(genes)
    all_dermatology_genes = list(set(all_dermatology_genes))  # Remove duplicates

    # Add additional common genes to reach ~200 total
    additional_genes = [
        f"GENE_{i:03d}" for i in range(1, 201 - len(all_dermatology_genes))
    ]
    all_genes = all_dermatology_genes + additional_genes

    # Define differential expression patterns
    # Upregulated in Treatment (wound healing response)
    upregulated_genes = {"VEGFA", "COL1A1", "FN1", "MMP1", "FGF2", "TGFB1"}
    # Downregulated in Treatment (barrier disruption)
    downregulated_genes = {"FLG", "LOR", "CLDN1", "IVL", "TIMP1"}

    # Generate base means (log-normal distribution for realistic counts)
    # Most genes have low expression, some have high
    base_means = np.random.lognormal(mean=3, sigma=1.5, size=len(all_genes))

    # Create count matrix
    counts_data = np.zeros((len(sample_names), len(all_genes)), dtype=int)

    for gene_idx, gene in enumerate(all_genes):
        base_mean = base_means[gene_idx]

        for sample_idx, sample in enumerate(sample_names):
            is_treatment = "Treatment" in sample

            # Apply differential expression
            if gene in upregulated_genes and is_treatment:
                # Upregulate: 2-3x fold change
                mean = base_mean * np.random.uniform(2.0, 3.0)
            elif gene in downregulated_genes and is_treatment:
                # Downregulate: 0.3-0.5x fold change
                mean = base_mean * np.random.uniform(0.3, 0.5)
            else:
                mean = base_mean

            # Add biological noise (replicate variation)
            mean = mean * np.random.normal(1.0, 0.1)

            # Generate count from Poisson distribution
            count = np.random.poisson(max(1, mean))
            counts_data[sample_idx, gene_idx] = count

    # Create DataFrames
    counts_df = pd.DataFrame(
        counts_data, index=sample_names, columns=all_genes, dtype=int
    )

    # Create metadata
    metadata_df = pd.DataFrame(
        {
            "condition": [
                "Control" if "Control" in s else "Treatment" for s in sample_names
            ]
        },
        index=sample_names,
    )

    return counts_df, metadata_df


def get_demo_description() -> str:
    """
    Get markdown description of the demo dataset.

    Returns:
        Markdown string describing dataset characteristics, design, and DE patterns
    """
    description = """# RNA-seq Demo Dataset

## Overview
Realistic RNA-seq dataset for dermatology/cosmetics research with built-in differential expression patterns.

## Experimental Design
- **Samples**: 6 total (3 Control replicates, 3 Treatment replicates)
- **Genes**: ~200 genes including curated dermatology panels
- **Conditions**: Control vs Treatment

### Sample Names
- Control: `Control_Rep1`, `Control_Rep2`, `Control_Rep3`
- Treatment: `Treatment_Rep1`, `Treatment_Rep2`, `Treatment_Rep3`

## Gene Coverage
Includes genes from all dermatology panels:
- **Anti-aging**: COL1A1, COL3A1, ELN, FBN1, MMP1, MMP2, MMP9, TIMP1
- **Skin Barrier**: FLG, LOR, IVL, CLDN1, CLDN4, TJP1, AQP3, HAS2, HAS3
- **Anti-inflammation**: IL1A, IL1B, IL6, IL8, TNFA, CXCL1, CCL2, PTGS2
- **Whitening/Melanogenesis**: MITF, TYR, TYRP1, DCT, MC1R, PMEL, OCA2
- **Sebum Regulation**: SREBF1, PPARG, FASN, SCD, DGAT2, ELOVL6
- **Wound Healing**: VEGFA, FGF2, TGFB1, EGFR, KRT6A, KRT16, KRT17, MMP3
- **Additional**: ~100 background genes for realistic analysis

## Differential Expression Patterns

### Upregulated in Treatment (Wound Healing Response)
- **VEGFA**: Vascular endothelial growth factor (angiogenesis)
- **COL1A1**: Type I collagen (ECM remodeling)
- **FN1**: Fibronectin (cell adhesion, ECM)
- **MMP1**: Matrix metalloproteinase 1 (ECM degradation)
- **FGF2**: Fibroblast growth factor 2 (proliferation)
- **TGFB1**: TGF-beta 1 (wound healing, fibrosis)

**Fold Change**: 2-3x upregulation

### Downregulated in Treatment (Barrier Disruption)
- **FLG**: Filaggrin (stratum corneum protein)
- **LOR**: Loricrin (stratum corneum protein)
- **CLDN1**: Claudin-1 (tight junctions)
- **IVL**: Involucrin (stratum corneum protein)
- **TIMP1**: TIMP metalloproteinase inhibitor 1

**Fold Change**: 0.3-0.5x downregulation

## Data Characteristics
- **Count Distribution**: Poisson with log-normal base means
- **Realistic Variation**: Biological noise added to replicate variation
- **Reproducibility**: Generated with `np.random.seed(42)`
- **Data Type**: Integer counts (raw, not normalized)

## Usage
```python
from demo_data import load_demo_dataset, get_demo_description

# Load dataset
counts_df, metadata_df = load_demo_dataset()

# Get description
description = get_demo_description()
print(description)

# Explore
print(f"Counts shape: {counts_df.shape}")  # (6, ~200)
print(f"Metadata shape: {metadata_df.shape}")  # (6, 1)
print(counts_df.head())
print(metadata_df)
```

## Expected Analysis Results
When running differential expression analysis (Treatment vs Control):
- **Upregulated genes** should show positive log2FoldChange (~1-1.5)
- **Downregulated genes** should show negative log2FoldChange (~-0.5 to -1)
- **Background genes** should show minimal change (log2FC ≈ 0)
- **Adjusted p-values** should be significant for DE genes (padj < 0.05)

## Notes
- This is synthetic data for demonstration and testing purposes
- Real biological experiments will have different patterns
- Use for validating analysis pipelines and UI functionality
"""
    return description
