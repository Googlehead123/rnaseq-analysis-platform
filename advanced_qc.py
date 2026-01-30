"""
Advanced Quality Control Module for RNA-seq Analysis Platform.

Provides outlier detection in PCA space and batch effect assessment.
All plots use Plotly (no matplotlib).
"""

from typing import List, Dict, Tuple, Optional
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.express as px


class AdvancedQC:
    """
    Advanced QC analyses for RNA-seq data.
    
    Features:
    - Outlier detection in PCA space (Euclidean distance from centroid)
    - Batch effect assessment (variance decomposition per PC)
    """
    
    def detect_outliers(
        self,
        pca_coords: pd.DataFrame,
        threshold_sd: float = 3.0
    ) -> Tuple[List[str], pd.Series]:
        """
        Detect outlier samples in PCA space.
        
        Calculates Euclidean distance of each sample from the centroid
        in PC1-PC2 space. Samples beyond threshold_sd standard deviations
        from the mean distance are flagged as outliers.
        
        Parameters
        ----------
        pca_coords : pd.DataFrame
            DataFrame with 'PC1' and 'PC2' columns, sample names as index
        threshold_sd : float
            Number of standard deviations for outlier threshold (default: 3.0)
        
        Returns
        -------
        Tuple[List[str], pd.Series]
            (list of outlier sample names, Series of distances for all samples)
        """
        pc_data = pca_coords[['PC1', 'PC2']]
        
        centroid = pc_data.mean()
        distances = np.sqrt(((pc_data - centroid) ** 2).sum(axis=1))
        
        threshold = distances.median() + threshold_sd * distances.std()
        outliers = distances[distances > threshold].index.tolist()
        
        return outliers, distances
    
    def assess_batch_effects(
        self,
        pca_coords: pd.DataFrame,
        batch_labels: pd.Series
    ) -> Dict[str, float]:
        """
        Assess variance explained by a batch variable for each PC.
        
        Uses sum-of-squares decomposition (one-way ANOVA-like) to calculate
        what fraction of each PC's variance is explained by the batch variable.
        
        Parameters
        ----------
        pca_coords : pd.DataFrame
            DataFrame with PC columns (PC1, PC2, ...), sample names as index
        batch_labels : pd.Series
            Batch assignment for each sample (index must match pca_coords)
        
        Returns
        -------
        Dict[str, float]
            Mapping of PC name to fraction of variance explained by batch (0-1)
        """
        results = {}
        
        # Align batch labels to PCA coords
        common = pca_coords.index.intersection(batch_labels.index)
        batch = batch_labels.loc[common]
        
        for pc in pca_coords.columns:
            if not pc.startswith('PC'):
                continue
            
            pc_values = pca_coords.loc[common, pc]
            grand_mean = pc_values.mean()
            
            # Between-group sum of squares
            ss_between = 0.0
            for group_name in batch.unique():
                group_mask = batch == group_name
                group_vals = pc_values[group_mask]
                ss_between += len(group_vals) * (group_vals.mean() - grand_mean) ** 2
            
            # Total sum of squares
            ss_total = ((pc_values - grand_mean) ** 2).sum()
            
            if ss_total > 0:
                results[pc] = float(ss_between / ss_total)
            else:
                results[pc] = 0.0
        
        return results
    
    def create_outlier_plot(
        self,
        pca_coords: pd.DataFrame,
        sample_conditions: Dict[str, str],
        outliers: List[str],
        distances: pd.Series
    ) -> go.Figure:
        """
        Create PCA scatter plot highlighting outlier samples.
        
        Parameters
        ----------
        pca_coords : pd.DataFrame
            PCA coordinates with PC1, PC2 columns
        sample_conditions : Dict[str, str]
            Mapping of sample name to condition
        outliers : List[str]
            List of outlier sample names
        distances : pd.Series
            Distance from centroid for each sample
        
        Returns
        -------
        go.Figure
            Plotly figure with outliers highlighted
        """
        df = pca_coords[['PC1', 'PC2']].copy()
        df['condition'] = [sample_conditions.get(s, 'Unknown') for s in df.index]
        df['distance'] = distances
        df['is_outlier'] = [s in outliers for s in df.index]
        df['sample'] = df.index
        
        # Plot normal samples by condition
        fig = px.scatter(
            df[~df['is_outlier']],
            x='PC1', y='PC2',
            color='condition',
            hover_name='sample',
            hover_data={'distance': ':.2f', 'condition': True, 'PC1': ':.2f', 'PC2': ':.2f'},
            title='PCA Outlier Detection'
        )
        
        # Overlay outliers with distinct markers
        if outliers:
            outlier_df = df[df['is_outlier']]
            fig.add_trace(go.Scatter(
                x=outlier_df['PC1'],
                y=outlier_df['PC2'],
                mode='markers+text',
                marker=dict(
                    size=14,
                    color='red',
                    symbol='x',
                    line=dict(width=2, color='darkred')
                ),
                text=outlier_df.index,
                textposition='top center',
                textfont=dict(size=10, color='red'),
                name='Outlier',
                hovertemplate=(
                    '<b>%{text}</b><br>'
                    'PC1: %{x:.2f}<br>'
                    'PC2: %{y:.2f}<br>'
                    'Distance: %{customdata:.2f}<extra>Outlier</extra>'
                ),
                customdata=outlier_df['distance'].values
            ))
        
        fig.update_layout(
            xaxis_title='PC1',
            yaxis_title='PC2',
            showlegend=True
        )
        
        return fig
    
    def create_batch_effect_plot(
        self,
        batch_results: Dict[str, float]
    ) -> go.Figure:
        """
        Create bar chart showing variance explained by batch per PC.
        
        Parameters
        ----------
        batch_results : Dict[str, float]
            Mapping of PC name to fraction of variance explained by batch
        
        Returns
        -------
        go.Figure
            Plotly bar chart
        """
        pcs = list(batch_results.keys())
        values = [batch_results[pc] * 100 for pc in pcs]  # Convert to percentage
        
        # Color bars: red if > 30% variance explained by batch (concerning)
        colors = ['#e74c3c' if v > 30 else '#3498db' for v in values]
        
        fig = go.Figure(data=[
            go.Bar(
                x=pcs,
                y=values,
                marker_color=colors,
                text=[f'{v:.1f}%' for v in values],
                textposition='auto',
                hovertemplate='%{x}: %{y:.1f}% variance explained by batch<extra></extra>'
            )
        ])
        
        # Add warning threshold line
        fig.add_hline(
            y=30, line_dash="dash", line_color="red",
            annotation_text="Concerning threshold (30%)",
            annotation_position="top right"
        )
        
        fig.update_layout(
            title='Batch Effect Assessment: Variance Explained per PC',
            xaxis_title='Principal Component',
            yaxis_title='Variance Explained by Batch (%)',
            yaxis=dict(range=[0, 100]),
            showlegend=False
        )
        
        return fig
