#!/usr/bin/env python3
"""
Brain Imaging Data Exploration Script
======================================
This script performs exploratory data analysis on neuroimaging data including:
- Volume normalization by brain mask volume
- Variability analysis of relative volumes and template-spaced metrics
- Longitudinal trend analysis
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime
import re
from pathlib import Path

# Set style for better-looking plots
sns.set_style("whitegrid")
plt.rcParams['figure.figsize'] = (12, 8)

def parse_subject_id(subject_id):
    """
    Parse subject_id to extract subject, project, and date information.
    
    Format: <subject>_<project>_<date> or <subject>_<project>
    """
    parts = subject_id.split('_')
    
    if len(parts) >= 3:
        # Check if last part is a date (8 digits)
        if parts[-1].isdigit() and len(parts[-1]) == 8:
            subject = '_'.join(parts[:-2])
            project = parts[-2]
            date_str = parts[-1]
            try:
                date = datetime.strptime(date_str, '%Y%m%d')
            except ValueError:
                date = None
            return subject, project, date_str, date
    
    # Baseline or non-standard format
    if len(parts) >= 2:
        subject = '_'.join(parts[:-1])
        project = parts[-1]
        return subject, project, 'baseline', None
    
    return subject_id, 'unknown', 'baseline', None

def load_and_prepare_data(csv_path):
    """Load CSV and prepare data with parsed subject information."""
    df = pd.read_csv(csv_path)
    
    # Parse subject IDs
    parsed = df['subject_id'].apply(parse_subject_id)
    df['subject'] = parsed.apply(lambda x: x[0])
    df['project'] = parsed.apply(lambda x: x[1])
    df['date_str'] = parsed.apply(lambda x: x[2])
    df['date'] = parsed.apply(lambda x: x[3])
    
    # Identify volume columns (need normalization)
    volume_cols = [col for col in df.columns if '_volume_mm3' in col and col != 'brain_mask_volume_mm3']
    
    # Normalize volumes by brain mask volume
    for col in volume_cols:
        normalized_col = col.replace('_volume_mm3', '_norm')
        df[normalized_col] = df[col] / df['brain_mask_volume_mm3']
    
    return df, volume_cols

def plot_volume_variability(df, volume_cols, output_dir='./plots'):
    """Plot variability of normalized volumes."""
    Path(output_dir).mkdir(exist_ok=True)
    
    # Get normalized column names
    norm_cols = [col.replace('_volume_mm3', '_norm') for col in volume_cols]
    
    # Create figure with subplots
    n_cols = len(norm_cols)
    fig, axes = plt.subplots(2, 1, figsize=(14, 10))
    
    # Box plot of normalized volumes
    df[norm_cols].boxplot(ax=axes[0], vert=True, patch_artist=True)
    axes[0].set_title('Normalized Volume Distributions (relative to brain mask volume)', fontsize=14, fontweight='bold')
    axes[0].set_ylabel('Normalized Volume (fraction of brain mask)')
    axes[0].tick_params(axis='x', rotation=45)
    axes[0].grid(True, alpha=0.3)
    
    # Violin plot for better distribution visualization
    norm_data = df[norm_cols].melt(var_name='Region', value_name='Normalized Volume')
    sns.violinplot(data=norm_data, x='Region', y='Normalized Volume', ax=axes[1])
    axes[1].set_title('Normalized Volume Distributions (Violin Plot)', fontsize=14, fontweight='bold')
    axes[1].tick_params(axis='x', rotation=45)
    axes[1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/normalized_volume_variability.png', dpi=300, bbox_inches='tight')
    print(f"Saved: {output_dir}/normalized_volume_variability.png")
    
    # Summary statistics
    print("\n" + "="*70)
    print("NORMALIZED VOLUME STATISTICS")
    print("="*70)
    print(df[norm_cols].describe())
    
    return fig

def plot_template_variability(df, output_dir='./plots'):
    """Plot variability of template-spaced metrics (DTI and NODDI)."""
    Path(output_dir).mkdir(exist_ok=True)
    
    # DTI FA metrics
    fa_cols = [col for col in df.columns if 'DTI_FA_' in col]
    # DTI MD metrics
    md_cols = [col for col in df.columns if 'DTI_MD_' in col]
    # NODDI metrics
    noddi_icvf_cols = [col for col in df.columns if 'NODDI_ICVF_' in col]
    noddi_isovf_cols = [col for col in df.columns if 'NODDI_ISOVF_' in col]
    noddi_od_cols = [col for col in df.columns if 'NODDI_OD_' in col]
    
    # Create comprehensive figure
    fig, axes = plt.subplots(3, 2, figsize=(16, 14))
    
    # DTI FA
    if fa_cols:
        fa_data = df[fa_cols].melt(var_name='Region', value_name='FA Value')
        sns.boxplot(data=fa_data, x='Region', y='FA Value', ax=axes[0, 0])
        axes[0, 0].set_title('DTI Fractional Anisotropy (FA)', fontsize=12, fontweight='bold')
        axes[0, 0].tick_params(axis='x', rotation=45)
        axes[0, 0].grid(True, alpha=0.3)
    
    # DTI MD
    if md_cols:
        md_data = df[md_cols].melt(var_name='Region', value_name='MD Value')
        sns.boxplot(data=md_data, x='Region', y='MD Value', ax=axes[0, 1])
        axes[0, 1].set_title('DTI Mean Diffusivity (MD)', fontsize=12, fontweight='bold')
        axes[0, 1].tick_params(axis='x', rotation=45)
        axes[0, 1].grid(True, alpha=0.3)
    
    # NODDI ICVF
    if noddi_icvf_cols:
        icvf_data = df[noddi_icvf_cols].melt(var_name='Region', value_name='ICVF Value')
        sns.boxplot(data=icvf_data, x='Region', y='ICVF Value', ax=axes[1, 0])
        axes[1, 0].set_title('NODDI Intracellular Volume Fraction (ICVF)', fontsize=12, fontweight='bold')
        axes[1, 0].tick_params(axis='x', rotation=45)
        axes[1, 0].grid(True, alpha=0.3)
    
    # NODDI ISOVF
    if noddi_isovf_cols:
        isovf_data = df[noddi_isovf_cols].melt(var_name='Region', value_name='ISOVF Value')
        sns.boxplot(data=isovf_data, x='Region', y='ISOVF Value', ax=axes[1, 1])
        axes[1, 1].set_title('NODDI Isotropic Volume Fraction (ISOVF)', fontsize=12, fontweight='bold')
        axes[1, 1].tick_params(axis='x', rotation=45)
        axes[1, 1].grid(True, alpha=0.3)
    
    # NODDI OD
    if noddi_od_cols:
        od_data = df[noddi_od_cols].melt(var_name='Region', value_name='OD Value')
        sns.boxplot(data=od_data, x='Region', y='OD Value', ax=axes[2, 0])
        axes[2, 0].set_title('NODDI Orientation Dispersion (OD)', fontsize=12, fontweight='bold')
        axes[2, 0].tick_params(axis='x', rotation=45)
        axes[2, 0].grid(True, alpha=0.3)
    
    # Gray matter probability
    if 'gm_mean_probability' in df.columns:
        axes[2, 1].hist(df['gm_mean_probability'].dropna(), bins=30, edgecolor='black', alpha=0.7)
        axes[2, 1].set_title('Gray Matter Mean Probability Distribution', fontsize=12, fontweight='bold')
        axes[2, 1].set_xlabel('Probability')
        axes[2, 1].set_ylabel('Frequency')
        axes[2, 1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/template_metrics_variability.png', dpi=300, bbox_inches='tight')
    print(f"Saved: {output_dir}/template_metrics_variability.png")
    
    # Print statistics
    print("\n" + "="*70)
    print("TEMPLATE-SPACED METRICS STATISTICS")
    print("="*70)
    
    if fa_cols:
        print("\nDTI FA Statistics:")
        print(df[fa_cols].describe())
    
    if md_cols:
        print("\nDTI MD Statistics:")
        print(df[md_cols].describe())
    
    return fig

def plot_longitudinal_trends(df, output_dir='./plots'):
    """Plot longitudinal trends for subjects with multiple timepoints."""
    Path(output_dir).mkdir(exist_ok=True)
    
    # Filter subjects with longitudinal data (multiple timepoints)
    longitudinal_subjects = df.groupby('subject').filter(lambda x: len(x) > 1 and x['date'].notna().any())
    
    if len(longitudinal_subjects) == 0:
        print("\nNo longitudinal data found (subjects with multiple dated timepoints).")
        return None
    
    # Sort by subject and date
    longitudinal_subjects = longitudinal_subjects.sort_values(['subject', 'date'])
    
    print(f"\n{'='*70}")
    print(f"LONGITUDINAL ANALYSIS")
    print(f"{'='*70}")
    print(f"Found {longitudinal_subjects['subject'].nunique()} subjects with longitudinal data")
    print(f"Total longitudinal scans: {len(longitudinal_subjects)}")
    
    # Calculate days from baseline for each subject
    def calc_days_from_baseline(group):
        if group['date'].notna().any():
            baseline_date = group['date'].min()
            group['days_from_baseline'] = (group['date'] - baseline_date).dt.days
        else:
            group['days_from_baseline'] = 0
        return group
    
    longitudinal_subjects = longitudinal_subjects.groupby('subject').apply(calc_days_from_baseline).reset_index(drop=True)
    
    # Select key metrics to plot
    volume_norm_cols = [col for col in longitudinal_subjects.columns if '_norm' in col and 'hypothalamus' in col]
    dti_fa_cols = [col for col in longitudinal_subjects.columns if 'DTI_FA_' in col][:3]  # First 3
    
    # Create longitudinal plot
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    
    # Plot 1: Normalized hypothalamus volumes over time
    if volume_norm_cols:
        ax = axes[0, 0]
        for subject in longitudinal_subjects['subject'].unique():
            subject_data = longitudinal_subjects[longitudinal_subjects['subject'] == subject]
            for col in volume_norm_cols:
                ax.plot(subject_data['days_from_baseline'], subject_data[col], 
                       marker='o', alpha=0.6, label=f"{subject}_{col.replace('_norm', '')}")
        ax.set_title('Normalized Hypothalamus Volumes Over Time', fontsize=12, fontweight='bold')
        ax.set_xlabel('Days from Baseline')
        ax.set_ylabel('Normalized Volume')
        ax.grid(True, alpha=0.3)
        # Only show legend if not too many lines
        if len(volume_norm_cols) * longitudinal_subjects['subject'].nunique() < 10:
            ax.legend(fontsize=8, loc='best')
    
    # Plot 2: DTI FA metrics over time
    if dti_fa_cols:
        ax = axes[0, 1]
        for subject in longitudinal_subjects['subject'].unique():
            subject_data = longitudinal_subjects[longitudinal_subjects['subject'] == subject]
            for col in dti_fa_cols:
                ax.plot(subject_data['days_from_baseline'], subject_data[col], 
                       marker='s', alpha=0.6, label=f"{subject}_{col}")
        ax.set_title('DTI FA Metrics Over Time', fontsize=12, fontweight='bold')
        ax.set_xlabel('Days from Baseline')
        ax.set_ylabel('FA Value')
        ax.grid(True, alpha=0.3)
        if len(dti_fa_cols) * longitudinal_subjects['subject'].nunique() < 10:
            ax.legend(fontsize=8, loc='best')
    
    # Plot 3: Brain mask volume over time
    ax = axes[1, 0]
    for subject in longitudinal_subjects['subject'].unique():
        subject_data = longitudinal_subjects[longitudinal_subjects['subject'] == subject]
        ax.plot(subject_data['days_from_baseline'], subject_data['brain_mask_volume_mm3'], 
               marker='o', alpha=0.7, linewidth=2, label=subject)
    ax.set_title('Brain Mask Volume Over Time', fontsize=12, fontweight='bold')
    ax.set_xlabel('Days from Baseline')
    ax.set_ylabel('Volume (mm³)')
    ax.grid(True, alpha=0.3)
    if longitudinal_subjects['subject'].nunique() < 15:
        ax.legend(fontsize=8, loc='best')
    
    # Plot 4: Individual subject trajectory heatmap
    ax = axes[1, 1]
    # Create a simple summary: number of scans per subject over time
    scan_counts = longitudinal_subjects.groupby(['subject', 'days_from_baseline']).size().reset_index(name='count')
    pivot_data = scan_counts.pivot(index='subject', columns='days_from_baseline', values='count').fillna(0)
    
    if len(pivot_data) > 0:
        sns.heatmap(pivot_data, cmap='YlOrRd', ax=ax, cbar_kws={'label': 'Scan Present'}, 
                   linewidths=0.5, linecolor='gray')
        ax.set_title('Longitudinal Scan Timeline', fontsize=12, fontweight='bold')
        ax.set_xlabel('Days from Baseline')
        ax.set_ylabel('Subject')
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/longitudinal_trends.png', dpi=300, bbox_inches='tight')
    print(f"Saved: {output_dir}/longitudinal_trends.png")
    
    # Summary statistics
    print("\nLongitudinal Follow-up Summary:")
    follow_up_stats = longitudinal_subjects.groupby('subject')['days_from_baseline'].agg(['count', 'min', 'max'])
    follow_up_stats['duration_days'] = follow_up_stats['max'] - follow_up_stats['min']
    print(follow_up_stats)
    
    return fig

def generate_summary_report(df, volume_cols, output_dir='./plots'):
    """Generate a text summary report."""
    Path(output_dir).mkdir(exist_ok=True)
    
    report_path = f'{output_dir}/summary_report.txt'
    
    with open(report_path, 'w') as f:
        f.write("="*70 + "\n")
        f.write("BRAIN IMAGING DATA EXPLORATION REPORT\n")
        f.write("="*70 + "\n\n")
        
        f.write(f"Total number of scans: {len(df)}\n")
        f.write(f"Number of unique subjects: {df['subject'].nunique()}\n")
        f.write(f"Number of projects: {df['project'].nunique()}\n")
        f.write(f"Projects: {', '.join(df['project'].unique())}\n\n")
        
        # Longitudinal info
        longitudinal_count = df.groupby('subject').filter(lambda x: len(x) > 1)['subject'].nunique()
        f.write(f"Subjects with multiple scans: {longitudinal_count}\n\n")
        
        f.write("="*70 + "\n")
        f.write("NORMALIZED VOLUME STATISTICS\n")
        f.write("="*70 + "\n")
        norm_cols = [col.replace('_volume_mm3', '_norm') for col in volume_cols]
        f.write(df[norm_cols].describe().to_string())
        f.write("\n\n")
        
        f.write("="*70 + "\n")
        f.write("BRAIN MASK VOLUME STATISTICS\n")
        f.write("="*70 + "\n")
        f.write(df['brain_mask_volume_mm3'].describe().to_string())
        f.write("\n\n")
        
        # DTI statistics
        fa_cols = [col for col in df.columns if 'DTI_FA_' in col]
        if fa_cols:
            f.write("="*70 + "\n")
            f.write("DTI FA STATISTICS\n")
            f.write("="*70 + "\n")
            f.write(df[fa_cols].describe().to_string())
            f.write("\n\n")
    
    print(f"\nSaved: {report_path}")
    return report_path

def main():
    """Main execution function."""
    print("="*70)
    print("BRAIN IMAGING DATA EXPLORATION")
    print("="*70)
    print("\nPlease provide the path to your CSV file:")
    csv_path = input("CSV path: ").strip()
    
    # Load and prepare data
    print("\nLoading and preparing data...")
    df, volume_cols = load_and_prepare_data(csv_path)
    print(f"Loaded {len(df)} scans from {df['subject'].nunique()} subjects")
    
    # Create output directory
    output_dir = './brain_data_plots'
    Path(output_dir).mkdir(exist_ok=True)
    print(f"\nOutput directory: {output_dir}")
    
    # Generate plots
    print("\n" + "="*70)
    print("GENERATING VISUALIZATIONS")
    print("="*70)
    
    print("\n1. Analyzing normalized volume variability...")
    plot_volume_variability(df, volume_cols, output_dir)
    
    print("\n2. Analyzing template-spaced metrics variability...")
    plot_template_variability(df, output_dir)
    
    print("\n3. Analyzing longitudinal trends...")
    plot_longitudinal_trends(df, output_dir)
    
    print("\n4. Generating summary report...")
    generate_summary_report(df, volume_cols, output_dir)
    
    print("\n" + "="*70)
    print("ANALYSIS COMPLETE!")
    print("="*70)
    print(f"\nAll outputs saved to: {output_dir}/")
    print("\nGenerated files:")
    print("  - normalized_volume_variability.png")
    print("  - template_metrics_variability.png")
    print("  - longitudinal_trends.png")
    print("  - summary_report.txt")

if __name__ == "__main__":
    main()