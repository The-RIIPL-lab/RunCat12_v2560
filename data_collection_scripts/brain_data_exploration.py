#!/usr/bin/env python3
"""
Data Exploration and Visualization Script for CAT12 Hypothalamus Data

This script provides cross-sectional visualization and analysis of:
- Volume measurements
- DTI metrics
- NODDI metrics
For hypothalamus ROIs and brain tissues

For longitudinal subjects (multiple timepoints), values are averaged.

Author: Data exploration for simplified hypothalamus analysis
Date: December 2025
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import argparse
import warnings
warnings.filterwarnings('ignore')


# Set plotting style
sns.set_style("whitegrid")
plt.rcParams['figure.figsize'] = (12, 8)
plt.rcParams['font.size'] = 10


def load_and_prepare_data(csv_path, average_longitudinal=True):
    """
    Load CSV data and prepare for cross-sectional analysis.
    
    Parameters:
    -----------
    csv_path : str
        Path to the CSV file
    average_longitudinal : bool
        If True, average values for subjects with multiple timepoints
    
    Returns:
    --------
    df : pandas.DataFrame
        Prepared dataframe
    """
    df = pd.read_csv(csv_path)
    
    print(f"Loaded {len(df)} rows from {csv_path}")
    
    # Remove failed/skipped subjects
    if 'processing_status' in df.columns:
        df_valid = df[df['processing_status'] == 'completed'].copy()
        print(f"Valid subjects: {len(df_valid)} (removed {len(df) - len(df_valid)} failed/skipped)")
        df = df_valid
    
    if average_longitudinal and len(df) > 0:
        # Extract base subject ID (remove timepoint suffixes like _v1, _v2, _baseline, etc.)
        df['base_subject_id'] = df['subject_id'].str.replace(r'_v\d+$', '', regex=True)
        df['base_subject_id'] = df['base_subject_id'].str.replace(r'_(baseline|followup|fu\d+)$', '', regex=True)
        
        # Count how many subjects have multiple timepoints
        subject_counts = df['base_subject_id'].value_counts()
        longitudinal_subjects = subject_counts[subject_counts > 1]
        
        if len(longitudinal_subjects) > 0:
            print(f"\nLongitudinal subjects detected: {len(longitudinal_subjects)} subjects with multiple timepoints")
            print(f"Total scans from these subjects: {longitudinal_subjects.sum()}")
            
            # Average numeric columns for each base subject
            numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()
            numeric_cols = [c for c in numeric_cols if c not in ['base_subject_id']]
            
            # Group by base subject and average
            df_averaged = df.groupby('base_subject_id')[numeric_cols].mean().reset_index()
            df_averaged.rename(columns={'base_subject_id': 'subject_id'}, inplace=True)
            
            print(f"After averaging: {len(df_averaged)} unique subjects")
            df = df_averaged
        else:
            df = df.drop(columns=['base_subject_id'])
            print("No longitudinal subjects detected")
    
    return df


def get_metric_columns(df, metric_type='volume'):
    """Get columns for specific metric type."""
    if metric_type == 'volume':
        return [c for c in df.columns if c.endswith('_volume_mm3')]
    elif metric_type == 'native_dti':
        return [c for c in df.columns if c.startswith('native_DTI_')]
    elif metric_type == 'native_noddi':
        return [c for c in df.columns if c.startswith('native_NODDI_')]
    elif metric_type == 'normalized_dti':
        return [c for c in df.columns if c.startswith('normalized_DTI_')]
    elif metric_type == 'normalized_noddi':
        return [c for c in df.columns if c.startswith('normalized_NODDI_')]
    return []


def print_summary_statistics(df):
    """Print summary statistics for all metrics."""
    print("\n" + "=" * 80)
    print("SUMMARY STATISTICS")
    print("=" * 80)
    
    metric_groups = {
        'Volumes (mm³)': get_metric_columns(df, 'volume'),
        'Native DTI Metrics': get_metric_columns(df, 'native_dti'),
        'Native NODDI Metrics': get_metric_columns(df, 'native_noddi'),
        'Normalized DTI Metrics': get_metric_columns(df, 'normalized_dti'),
        'Normalized NODDI Metrics': get_metric_columns(df, 'normalized_noddi'),
    }
    
    for group_name, cols in metric_groups.items():
        if not cols:
            continue
        
        print(f"\n{group_name}")
        print("-" * 80)
        
        for col in cols:
            if col in df.columns:
                data = df[col].dropna()
                if len(data) > 0:
                    print(f"{col:55s}: "
                          f"n={len(data):4d}, "
                          f"mean={np.mean(data):10.4f}, "
                          f"std={np.std(data):10.4f}, "
                          f"min={np.min(data):10.4f}, "
                          f"max={np.max(data):10.4f}")


def plot_volume_distributions(df, output_dir):
    """Create distribution plots for volume measurements."""
    volume_cols = get_metric_columns(df, 'volume')
    if not volume_cols:
        print("No volume columns found")
        return
    
    # Create figure
    n_cols = len(volume_cols)
    n_rows = (n_cols + 2) // 3  # 3 columns per row
    fig, axes = plt.subplots(n_rows, 3, figsize=(15, 5*n_rows))
    axes = axes.flatten() if n_rows > 1 else [axes]
    
    for idx, col in enumerate(volume_cols):
        data = df[col].dropna()
        if len(data) > 0:
            axes[idx].hist(data, bins=30, edgecolor='black', alpha=0.7)
            axes[idx].set_xlabel('Volume (mm³)')
            axes[idx].set_ylabel('Frequency')
            axes[idx].set_title(col.replace('_volume_mm3', '').replace('_', ' ').title())
            axes[idx].grid(True, alpha=0.3)
    
    # Hide unused subplots
    for idx in range(len(volume_cols), len(axes)):
        axes[idx].set_visible(False)
    
    plt.tight_layout()
    output_path = Path(output_dir) / 'volume_distributions.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_path}")
    plt.close()


def plot_hypothalamus_roi_comparison(df, output_dir):
    """Compare hypothalamus ROI 1 vs ROI 2 across all metrics."""
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    
    # Volume comparison
    vol_data = []
    for roi in [1, 2]:
        col = f'hypothalamus_roi{roi}_volume_mm3'
        if col in df.columns:
            vol_data.append(df[col].dropna())
    
    if len(vol_data) == 2:
        axes[0, 0].boxplot(vol_data, labels=['ROI 1', 'ROI 2'])
        axes[0, 0].set_ylabel('Volume (mm³)')
        axes[0, 0].set_title('Hypothalamus Volume')
        axes[0, 0].grid(True, alpha=0.3)
    
    # Native DTI FA comparison
    dti_data = []
    for roi in [1, 2]:
        col = f'native_DTI_FA_hyp_roi{roi}_mean'
        if col in df.columns:
            dti_data.append(df[col].dropna())
    
    if len(dti_data) == 2:
        axes[0, 1].boxplot(dti_data, labels=['ROI 1', 'ROI 2'])
        axes[0, 1].set_ylabel('FA')
        axes[0, 1].set_title('Native DTI FA')
        axes[0, 1].grid(True, alpha=0.3)
    
    # Native NODDI ICVF comparison
    noddi_data = []
    for roi in [1, 2]:
        col = f'native_NODDI_ICVF_hyp_roi{roi}_mean'
        if col in df.columns:
            noddi_data.append(df[col].dropna())
    
    if len(noddi_data) == 2:
        axes[0, 2].boxplot(noddi_data, labels=['ROI 1', 'ROI 2'])
        axes[0, 2].set_ylabel('ICVF')
        axes[0, 2].set_title('Native NODDI ICVF')
        axes[0, 2].grid(True, alpha=0.3)
    
    # Normalized DTI FA comparison
    norm_dti_data = []
    for roi in [1, 2]:
        col = f'normalized_DTI_FA_hyp_roi{roi}_mean'
        if col in df.columns:
            norm_dti_data.append(df[col].dropna())
    
    if len(norm_dti_data) == 2:
        axes[1, 0].boxplot(norm_dti_data, labels=['ROI 1', 'ROI 2'])
        axes[1, 0].set_ylabel('FA')
        axes[1, 0].set_title('Normalized DTI FA')
        axes[1, 0].grid(True, alpha=0.3)
    
    # Normalized NODDI ICVF comparison
    norm_noddi_data = []
    for roi in [1, 2]:
        col = f'normalized_NODDI_ICVF_hyp_roi{roi}_mean'
        if col in df.columns:
            norm_noddi_data.append(df[col].dropna())
    
    if len(norm_noddi_data) == 2:
        axes[1, 1].boxplot(norm_noddi_data, labels=['ROI 1', 'ROI 2'])
        axes[1, 1].set_ylabel('ICVF')
        axes[1, 1].set_title('Normalized NODDI ICVF')
        axes[1, 1].grid(True, alpha=0.3)
    
    # ROI 1 vs ROI 2 volume scatter
    if 'hypothalamus_roi1_volume_mm3' in df.columns and 'hypothalamus_roi2_volume_mm3' in df.columns:
        roi1 = df['hypothalamus_roi1_volume_mm3'].dropna()
        roi2 = df['hypothalamus_roi2_volume_mm3'].dropna()
        common_idx = roi1.index.intersection(roi2.index)
        axes[1, 2].scatter(roi1[common_idx], roi2[common_idx], alpha=0.5)
        axes[1, 2].set_xlabel('ROI 1 Volume (mm³)')
        axes[1, 2].set_ylabel('ROI 2 Volume (mm³)')
        axes[1, 2].set_title('ROI 1 vs ROI 2 Volume')
        axes[1, 2].grid(True, alpha=0.3)
        
        # Add diagonal line
        max_val = max(roi1[common_idx].max(), roi2[common_idx].max())
        axes[1, 2].plot([0, max_val], [0, max_val], 'r--', alpha=0.5, label='y=x')
        axes[1, 2].legend()
    
    plt.tight_layout()
    output_path = Path(output_dir) / 'hypothalamus_roi_comparison.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_path}")
    plt.close()


def plot_correlation_matrix(df, output_dir, metric_type='native_dti'):
    """Create correlation matrix for specified metric type."""
    cols = get_metric_columns(df, metric_type)
    if len(cols) < 2:
        print(f"Not enough columns for {metric_type} correlation matrix")
        return
    
    # Filter to only ROI 1 and 2
    cols = [c for c in cols if 'roi1' in c or 'roi2' in c]
    
    # Calculate correlation matrix
    corr_data = df[cols].dropna()
    if len(corr_data) < 3:
        print(f"Not enough data for {metric_type} correlation matrix")
        return
    
    corr = corr_data.corr()
    
    # Simplify labels
    labels = [c.replace(f'{metric_type}_', '').replace('_mean', '').replace('_hyp_', ' ') 
              for c in cols]
    
    # Create heatmap
    plt.figure(figsize=(12, 10))
    sns.heatmap(corr, annot=True, fmt='.2f', cmap='coolwarm', center=0,
                xticklabels=labels, yticklabels=labels,
                cbar_kws={'label': 'Correlation'})
    plt.title(f'Correlation Matrix: {metric_type.replace("_", " ").title()}')
    plt.tight_layout()
    
    output_path = Path(output_dir) / f'correlation_matrix_{metric_type}.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_path}")
    plt.close()


def plot_tissue_volume_comparison(df, output_dir):
    """Create comparison plot for tissue volumes."""
    tissue_cols = ['gm_volume_mm3', 'wm_volume_mm3', 'csf_volume_mm3']
    tissue_cols = [c for c in tissue_cols if c in df.columns]
    
    if len(tissue_cols) < 2:
        print("Not enough tissue columns for comparison")
        return
    
    # Prepare data
    tissue_data = []
    tissue_labels = []
    for col in tissue_cols:
        data = df[col].dropna()
        if len(data) > 0:
            tissue_data.append(data)
            tissue_labels.append(col.replace('_volume_mm3', '').upper())
    
    # Create box plot
    fig, axes = plt.subplots(1, 2, figsize=(15, 6))
    
    # Box plot
    axes[0].boxplot(tissue_data, labels=tissue_labels)
    axes[0].set_ylabel('Volume (mm³)')
    axes[0].set_title('Brain Tissue Volume Distribution')
    axes[0].grid(True, alpha=0.3)
    
    # Stacked bar chart showing proportions
    if len(tissue_data) >= 2:
        # Calculate total for each subject
        tissue_df = df[tissue_cols].copy()
        tissue_df['total'] = tissue_df.sum(axis=1)
        
        # Calculate proportions
        for col in tissue_cols:
            tissue_df[f'{col}_prop'] = tissue_df[col] / tissue_df['total'] * 100
        
        # Plot mean proportions
        prop_cols = [f'{col}_prop' for col in tissue_cols]
        means = tissue_df[prop_cols].mean()
        stds = tissue_df[prop_cols].std()
        
        x = np.arange(len(tissue_labels))
        axes[1].bar(x, means, yerr=stds, capsize=5, alpha=0.7)
        axes[1].set_xticks(x)
        axes[1].set_xticklabels(tissue_labels)
        axes[1].set_ylabel('Proportion (%)')
        axes[1].set_title('Mean Tissue Proportions')
        axes[1].grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    output_path = Path(output_dir) / 'tissue_volume_comparison.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_path}")
    plt.close()


def plot_dti_metrics_overview(df, output_dir, space='native'):
    """Create overview of DTI metrics for both ROIs."""
    prefix = f'{space}_DTI_'
    dti_metrics = ['FA', 'MD', 'L1', 'L2', 'L3']
    
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    axes = axes.flatten()
    
    for idx, metric in enumerate(dti_metrics):
        data_roi1 = []
        data_roi2 = []
        
        col1 = f'{prefix}{metric}_hyp_roi1_mean'
        col2 = f'{prefix}{metric}_hyp_roi2_mean'
        
        if col1 in df.columns:
            data_roi1 = df[col1].dropna().values
        if col2 in df.columns:
            data_roi2 = df[col2].dropna().values
        
        if len(data_roi1) > 0 and len(data_roi2) > 0:
            axes[idx].boxplot([data_roi1, data_roi2], labels=['ROI 1', 'ROI 2'])
            axes[idx].set_title(f'{space.title()} DTI {metric}')
            axes[idx].set_ylabel(metric)
            axes[idx].grid(True, alpha=0.3)
    
    # Hide unused subplot
    axes[-1].set_visible(False)
    
    plt.suptitle(f'{space.title()} Space DTI Metrics Comparison', fontsize=14, y=1.00)
    plt.tight_layout()
    
    output_path = Path(output_dir) / f'dti_metrics_overview_{space}.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_path}")
    plt.close()


def plot_noddi_metrics_overview(df, output_dir, space='native'):
    """Create overview of NODDI metrics for both ROIs."""
    prefix = f'{space}_NODDI_'
    noddi_metrics = ['ICVF', 'ISOVF', 'OD']
    
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    
    for idx, metric in enumerate(noddi_metrics):
        data_roi1 = []
        data_roi2 = []
        
        col1 = f'{prefix}{metric}_hyp_roi1_mean'
        col2 = f'{prefix}{metric}_hyp_roi2_mean'
        
        if col1 in df.columns:
            data_roi1 = df[col1].dropna().values
        if col2 in df.columns:
            data_roi2 = df[col2].dropna().values
        
        if len(data_roi1) > 0 and len(data_roi2) > 0:
            axes[idx].boxplot([data_roi1, data_roi2], labels=['ROI 1', 'ROI 2'])
            axes[idx].set_title(f'{space.title()} NODDI {metric}')
            axes[idx].set_ylabel(metric)
            axes[idx].grid(True, alpha=0.3)
    
    plt.suptitle(f'{space.title()} Space NODDI Metrics Comparison', fontsize=14, y=1.00)
    plt.tight_layout()
    
    output_path = Path(output_dir) / f'noddi_metrics_overview_{space}.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_path}")
    plt.close()


def create_summary_report(df, output_dir):
    """Create a text summary report."""
    output_path = Path(output_dir) / 'summary_report.txt'
    
    with open(output_path, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("CAT12 HYPOTHALAMUS DATA - SUMMARY REPORT\n")
        f.write("=" * 80 + "\n\n")
        
        f.write(f"Total subjects: {len(df)}\n\n")
        
        # Volume summary
        f.write("VOLUME MEASUREMENTS (mm³)\n")
        f.write("-" * 80 + "\n")
        volume_cols = get_metric_columns(df, 'volume')
        for col in volume_cols:
            if col in df.columns:
                data = df[col].dropna()
                if len(data) > 0:
                    f.write(f"{col:50s}: "
                           f"n={len(data):4d}, "
                           f"mean={np.mean(data):10.2f} ± {np.std(data):8.2f}, "
                           f"range=[{np.min(data):10.2f}, {np.max(data):10.2f}]\n")
        
        # Native DTI summary
        f.write("\n\nNATIVE DTI METRICS\n")
        f.write("-" * 80 + "\n")
        dti_cols = get_metric_columns(df, 'native_dti')
        for col in dti_cols:
            if col in df.columns:
                data = df[col].dropna()
                if len(data) > 0:
                    f.write(f"{col:50s}: "
                           f"n={len(data):4d}, "
                           f"mean={np.mean(data):10.4f} ± {np.std(data):8.4f}\n")
        
        # Native NODDI summary
        f.write("\n\nNATIVE NODDI METRICS\n")
        f.write("-" * 80 + "\n")
        noddi_cols = get_metric_columns(df, 'native_noddi')
        for col in noddi_cols:
            if col in df.columns:
                data = df[col].dropna()
                if len(data) > 0:
                    f.write(f"{col:50s}: "
                           f"n={len(data):4d}, "
                           f"mean={np.mean(data):10.4f} ± {np.std(data):8.4f}\n")
        
        # Normalized DTI summary
        f.write("\n\nNORMALIZED DTI METRICS\n")
        f.write("-" * 80 + "\n")
        norm_dti_cols = get_metric_columns(df, 'normalized_dti')
        for col in norm_dti_cols:
            if col in df.columns:
                data = df[col].dropna()
                if len(data) > 0:
                    f.write(f"{col:50s}: "
                           f"n={len(data):4d}, "
                           f"mean={np.mean(data):10.4f} ± {np.std(data):8.4f}\n")
        
        # Normalized NODDI summary
        f.write("\n\nNORMALIZED NODDI METRICS\n")
        f.write("-" * 80 + "\n")
        norm_noddi_cols = get_metric_columns(df, 'normalized_noddi')
        for col in norm_noddi_cols:
            if col in df.columns:
                data = df[col].dropna()
                if len(data) > 0:
                    f.write(f"{col:50s}: "
                           f"n={len(data):4d}, "
                           f"mean={np.mean(data):10.4f} ± {np.std(data):8.4f}\n")
    
    print(f"Saved: {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description='Explore and visualize CAT12 hypothalamus data (cross-sectional)'
    )
    parser.add_argument('csv_file', help='Path to the CSV file generated by extract_roi_values.py')
    parser.add_argument('-o', '--output-dir', default='./analysis_output', 
                       help='Output directory for plots and reports')
    parser.add_argument('--no-average', action='store_true',
                       help='Do not average longitudinal subjects (use all timepoints)')
    
    args = parser.parse_args()
    
    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("=" * 80)
    print("CAT12 HYPOTHALAMUS DATA EXPLORATION")
    print("=" * 80)
    
    # Load and prepare data
    df = load_and_prepare_data(args.csv_file, average_longitudinal=not args.no_average)
    
    if len(df) == 0:
        print("No valid data to analyze")
        return
    
    # Print summary statistics
    print_summary_statistics(df)
    
    print("\n" + "=" * 80)
    print("GENERATING VISUALIZATIONS")
    print("=" * 80)
    
    # Generate plots
    plot_volume_distributions(df, output_dir)
    plot_hypothalamus_roi_comparison(df, output_dir)
    plot_tissue_volume_comparison(df, output_dir)
    plot_dti_metrics_overview(df, output_dir, space='native')
    plot_dti_metrics_overview(df, output_dir, space='normalized')
    plot_noddi_metrics_overview(df, output_dir, space='native')
    plot_noddi_metrics_overview(df, output_dir, space='normalized')
    
    # Correlation matrices
    plot_correlation_matrix(df, output_dir, 'native_dti')
    plot_correlation_matrix(df, output_dir, 'native_noddi')
    plot_correlation_matrix(df, output_dir, 'normalized_dti')
    plot_correlation_matrix(df, output_dir, 'normalized_noddi')
    
    # Create summary report
    create_summary_report(df, output_dir)
    
    print("\n" + "=" * 80)
    print(f"Analysis complete! All outputs saved to: {output_dir}")
    print("=" * 80)


if __name__ == "__main__":
    main()