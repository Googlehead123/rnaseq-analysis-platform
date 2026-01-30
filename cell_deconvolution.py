"""
Cell Type Deconvolution Module for RNA-seq Analysis Platform.

Estimates cell type proportions from bulk RNA-seq data using
Non-Negative Least Squares (NNLS) deconvolution with reference signature matrices.
"""

from dataclasses import dataclass
from typing import Optional, Dict
import pandas as pd
import numpy as np
from scipy.optimize import nnls


@dataclass
class DeconvolutionResult:
    """Container for deconvolution results."""
    proportions: pd.DataFrame  # samples × cell types, values sum to ~1 per sample
    method: str
    goodness_of_fit: Optional[pd.Series] = None  # Residual per sample


class CellTypeDeconvolution:
    """
    Cell type deconvolution from bulk RNA-seq data.
    
    Uses NNLS (Non-Negative Least Squares) to estimate the proportion of
    different cell types in each sample, given a reference signature matrix.
    """
    
    def run_nnls(
        self,
        expression: pd.DataFrame,
        signature: pd.DataFrame
    ) -> DeconvolutionResult:
        """
        Run NNLS deconvolution.
        
        Parameters
        ----------
        expression : pd.DataFrame
            Normalized expression matrix (genes × samples).
            Gene names as index, sample names as columns.
        signature : pd.DataFrame
            Reference signature matrix (genes × cell types).
            Gene names as index, cell type names as columns.
        
        Returns
        -------
        DeconvolutionResult
            Cell type proportions (samples × cell types)
        
        Raises
        ------
        ValueError
            If fewer than 50 genes overlap between expression and signature
        """
        # Find common genes
        common = expression.index.intersection(signature.index)
        
        if len(common) < 50:
            raise ValueError(
                f"Insufficient gene overlap for deconvolution: only {len(common)} genes "
                f"found in common between expression data ({len(expression)}) and "
                f"signature matrix ({len(signature)}). Need at least 50."
            )
        
        expr = expression.loc[common]
        sig = signature.loc[common]
        
        # NNLS for each sample
        proportions = []
        residuals = []
        
        for sample in expr.columns:
            coef, residual = nnls(sig.values, expr[sample].values)
            
            # Normalize to sum to 1
            total = coef.sum()
            if total > 0:
                coef_norm = coef / total
            else:
                coef_norm = coef
            
            proportions.append(coef_norm)
            residuals.append(residual)
        
        return DeconvolutionResult(
            proportions=pd.DataFrame(
                proportions,
                index=expr.columns,
                columns=sig.columns
            ),
            method='NNLS',
            goodness_of_fit=pd.Series(residuals, index=expr.columns, name='residual')
        )
    
    @staticmethod
    def get_default_skin_signature() -> pd.DataFrame:
        """
        Return a built-in signature matrix for skin tissue cell types.
        
        Cell types:
        - Keratinocyte
        - Fibroblast  
        - Melanocyte
        - T cell (Immune)
        - Endothelial
        
        Returns a genes × cell_types DataFrame with realistic relative expression values.
        Values represent relative expression (arbitrary units scaled like TPM).
        """
        # Marker genes with relative expression patterns
        # Each row: gene → expression across cell types
        # Values represent relative enrichment in each cell type
        signature_data = {
            # Gene: [Keratinocyte, Fibroblast, Melanocyte, T_cell, Endothelial]
            # Keratinocyte markers
            'KRT14': [100, 2, 1, 0, 1],
            'KRT10': [95, 1, 1, 0, 0],
            'KRT5': [90, 3, 1, 0, 1],
            'KRT1': [85, 1, 0, 0, 0],
            'IVL': [80, 1, 0, 0, 0],
            'LOR': [75, 0, 0, 0, 0],
            'FLG': [70, 0, 0, 0, 0],
            'DSG1': [65, 2, 0, 0, 1],
            'SPRR1A': [60, 0, 0, 0, 0],
            'CDSN': [55, 0, 0, 0, 0],
            # Fibroblast markers
            'COL1A1': [5, 100, 2, 0, 5],
            'COL1A2': [3, 95, 1, 0, 3],
            'COL3A1': [2, 90, 1, 0, 4],
            'VIM': [10, 85, 15, 20, 30],
            'FN1': [3, 80, 2, 1, 8],
            'DCN': [1, 75, 0, 0, 2],
            'LUM': [1, 70, 0, 0, 1],
            'FAP': [0, 65, 0, 0, 1],
            'PDGFRA': [1, 60, 1, 0, 3],
            'THY1': [0, 55, 0, 5, 2],
            # Melanocyte markers
            'MITF': [1, 0, 100, 0, 0],
            'TYR': [0, 0, 95, 0, 0],
            'TYRP1': [0, 0, 90, 0, 0],
            'DCT': [0, 0, 85, 0, 0],
            'PMEL': [0, 0, 80, 0, 0],
            'MLANA': [0, 0, 75, 0, 0],
            'SOX10': [0, 0, 70, 0, 0],
            'SLC45A2': [0, 0, 65, 0, 0],
            'EDNRB': [0, 1, 60, 0, 5],
            'KIT': [0, 1, 55, 5, 2],
            # T cell / Immune markers
            'CD3D': [0, 0, 0, 100, 0],
            'CD3E': [0, 0, 0, 95, 0],
            'CD2': [0, 0, 0, 90, 0],
            'CD8A': [0, 0, 0, 85, 0],
            'CD4': [0, 0, 0, 80, 0],
            'IL7R': [0, 0, 0, 75, 0],
            'PTPRC': [0, 0, 0, 70, 0],
            'LCK': [0, 0, 0, 65, 0],
            'GZMB': [0, 0, 0, 60, 0],
            'PRF1': [0, 0, 0, 55, 0],
            # Endothelial markers
            'PECAM1': [0, 2, 0, 1, 100],
            'CDH5': [0, 1, 0, 0, 95],
            'VWF': [0, 1, 0, 0, 90],
            'KDR': [0, 1, 0, 0, 85],
            'FLT1': [0, 0, 0, 0, 80],
            'TEK': [0, 0, 0, 0, 75],
            'ENG': [0, 2, 0, 0, 70],
            'EMCN': [0, 0, 0, 0, 65],
            'CLDN5': [0, 0, 0, 0, 60],
            'ERG': [0, 1, 0, 0, 55],
        }
        
        cell_types = ['Keratinocyte', 'Fibroblast', 'Melanocyte', 'T_cell', 'Endothelial']
        
        df = pd.DataFrame.from_dict(signature_data, orient='index', columns=cell_types)
        df = df.astype(float)
        
        return df
