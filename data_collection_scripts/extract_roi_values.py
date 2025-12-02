#!/usr/bin/env python3
"""
Post-processing script for CAT12 workflow
Extracts volumes and template-space diffusion values from normalized ROIs

This script processes CAT12 outputs to extract:
1. Volumes from native space segmentations
2. Diffusion metrics (DTI, NODDI) from template space with template atlases

Author: Refactored for volume and template-space analysis
Date: November 2025
"""

import numpy as np
import nibabel as nib
import pandas as pd
import os
from pathlib import Path
import argparse


# Template atlas paths
# UPDATE THESE
HYPOTHALAMUS_ATLAS_PATH = "/isilon/datalake/riipl/original/DEMONco/Hellcat-12.9/resampled_template_ROIs/rhypothalamusAtlas.nii"
JHU_ATLAS_PATH = "/isilon/datalake/riipl/original/DEMONco/Hellcat-12.9/resampled_template_ROIs/rJHU.nii"


def load_nifti_data(filepath):
    """Load NIfTI data and return data array and affine matrix."""
    try:
        img = nib.load(filepath)
        return img.get_fdata(), img.affine, img.header
    except Exception as e:
        print(f"Error loading {filepath}: {e}")
        return None, None, None


def calculate_volume_mm3(data, header):
    """Calculate volume in mm³ from binary mask or segmentation."""
    try:
        # Get voxel dimensions from header
        voxel_dims = header.get_zooms()[:3]  # Get x, y, z dimensions
        voxel_volume_mm3 = np.prod(voxel_dims)
        
        # Count non-zero voxels
        n_voxels = np.sum(data > 0)
        
        # Calculate total volume
        volume_mm3 = n_voxels * voxel_volume_mm3
        
        return volume_mm3, n_voxels, voxel_volume_mm3
    except Exception as e:
        print(f"Error calculating volume: {e}")
        return np.nan, 0, np.nan


def calculate_segmentation_volumes(data, header):
    """Calculate volumes for each segment in an atlas/segmentation."""
    try:
        voxel_dims = header.get_zooms()[:3]
        voxel_volume_mm3 = np.prod(voxel_dims)
        
        # Get unique ROI IDs (excluding 0 which is background)
        unique_ids = np.unique(data)
        unique_ids = unique_ids[unique_ids > 0]
        
        volumes = {}
        for roi_id in unique_ids:
            roi_id_int = int(roi_id)
            n_voxels = np.sum(data == roi_id)
            volume_mm3 = n_voxels * voxel_volume_mm3
            volumes[roi_id_int] = {
                'volume_mm3': volume_mm3,
                'n_voxels': int(n_voxels)
            }
        
        return volumes
    except Exception as e:
        print(f"Error calculating segmentation volumes: {e}")
        return {}


def extract_roi_mean_nonzero(data_img, atlas_img, roi_id):
    """Extract mean of non-zero values from ROI."""
    try:
        mask_data = atlas_img == roi_id
        if np.sum(mask_data) == 0:
            return {'mean_nonzero': np.nan, 'n_voxels': 0, 'n_nonzero_voxels': 0}
        
        roi_values = data_img[mask_data]
        
        # Get non-zero values
        nonzero_values = roi_values[roi_values != 0]
        # Remove NaN values
        nonzero_values = nonzero_values[~np.isnan(nonzero_values)]
        
        if len(nonzero_values) == 0:
            return {'mean_nonzero': np.nan, 'n_voxels': int(np.sum(mask_data)), 'n_nonzero_voxels': 0}
        
        return {
            'mean_nonzero': float(np.mean(nonzero_values)),
            'n_voxels': int(np.sum(mask_data)),
            'n_nonzero_voxels': int(len(nonzero_values))
        }
    except Exception as e:
        print(f"Error extracting ROI values: {e}")
        return {'mean_nonzero': np.nan, 'n_voxels': 0, 'n_nonzero_voxels': 0}


def get_jhu_rois():
    """Define JHU atlas ROI IDs for regions of interest."""
    # These IDs should be verified against your specific JHU atlas version
    # Common JHU-ICBM-DTI-81 ROI IDs:
    jhu_rois = {
        'cc_genu': 3,           # Genu of corpus callosum
        'cc_body': 4,           # Body of corpus callosum  
        'cc_splenium': 5,       # Splenium of corpus callosum
        # Add more ROIs as needed
    }
    return jhu_rois


def get_hypothalamus_rois():
    """Define hypothalamic atlas ROI IDs."""
    # These IDs should be verified against your specific hypothalamus atlas
    hyp_rois = {}
    
    # If your atlas has specific ROI IDs, add them here
    # For example, if left hypothalamus = 1, right = 2:
    # hyp_rois = {
    #     'hypothalamus_left': 1,
    #     'hypothalamus_right': 2,
    # }
    
    # If you want all hypothalamus ROIs combined:
    # We'll extract all non-zero ROIs and name them dynamically
    
    return hyp_rois


def validate_subject(subject_dir):
    """Validate that a subject has the required CAT12 outputs."""
    subject_dir = Path(subject_dir)
    
    # Check for CAT12 directory
    cat12_dir = subject_dir / 'nifti' / 'cat12_v2560'
    if not cat12_dir.exists():
        return False, "CAT12 output directory not found"
    
    # Check for basic required files
    mri_dir = cat12_dir / 'mri'
    if not mri_dir.exists():
        return False, "MRI directory not found"
    
    return True, "Valid"


def process_volumes(subject_dir):
    """Extract volume measurements from native space files."""
    results = {}
    cat12_dir = subject_dir / 'nifti' / 'cat12_v2560'
    
    # 1. Brain mask volume
    brain_mask_file = cat12_dir / 'wbrainmask_T1.nii'
    if brain_mask_file.exists():
        try:
            data, affine, header = load_nifti_data(str(brain_mask_file))
            if data is not None:
                volume_mm3, n_voxels, voxel_vol = calculate_volume_mm3(data, header)
                results['brain_mask_volume_mm3'] = volume_mm3
                results['brain_mask_n_voxels'] = n_voxels
        except Exception as e:
            results['brain_mask_error'] = str(e)
    
    # 2. Hypothalamus volumes
    hyp_atlas_file = cat12_dir / 'whypothalamusAtlas.nii'
    if hyp_atlas_file.exists():
        try:
            data, affine, header = load_nifti_data(str(hyp_atlas_file))
            if data is not None:
                seg_volumes = calculate_segmentation_volumes(data, header)
                
                # Add each hypothalamus ROI volume
                for roi_id, vol_info in seg_volumes.items():
                    results[f'hypothalamus_roi{roi_id}_volume_mm3'] = vol_info['volume_mm3']
                    results[f'hypothalamus_roi{roi_id}_n_voxels'] = vol_info['n_voxels']
                
                # Total hypothalamus volume (all non-zero voxels)
                total_vol, total_voxels, _ = calculate_volume_mm3(data, header)
                results['hypothalamus_total_volume_mm3'] = total_vol
                results['hypothalamus_total_n_voxels'] = total_voxels
        except Exception as e:
            results['hypothalamus_error'] = str(e)
    
    # 3. Tissue volumes (GM, WM, CSF)
    mri_dir = cat12_dir / 'mri'
    tissue_types = {
        'gm': 'p1',  # Gray matter
        'wm': 'p2',  # White matter
        'csf': 'p3'  # CSF
    }
    
    for tissue_name, prefix in tissue_types.items():
        tissue_files = list(mri_dir.glob(f'{prefix}*-tfl3d116ns.nii'))
        if not tissue_files:
            # Try alternative pattern
            tissue_files = list(mri_dir.glob(f'{prefix}*.nii'))
        
        if tissue_files:
            try:
                data, affine, header = load_nifti_data(str(tissue_files[0]))
                if data is not None:
                    # For tissue probability maps, threshold at 0.5
                    binary_data = (data > 0.5).astype(int)
                    volume_mm3, n_voxels, voxel_vol = calculate_volume_mm3(binary_data, header)
                    results[f'{tissue_name}_volume_mm3'] = volume_mm3
                    results[f'{tissue_name}_n_voxels'] = n_voxels
                    
                    # Also store mean probability
                    results[f'{tissue_name}_mean_probability'] = float(np.mean(data))
            except Exception as e:
                results[f'{tissue_name}_error'] = str(e)
    
    return results


def process_template_diffusion(subject_dir, hyp_atlas_data, jhu_atlas_data):
    """Process template-space diffusion data (DTI and NODDI)."""
    results = {}
    cat12_dir = subject_dir / 'nifti' / 'cat12_v2560'
    mri_dir = cat12_dir / 'mri'
    
    # Define diffusion metrics to extract
    dti_metrics = ['FA', 'MD', 'L1', 'L2', 'L3', 'V1', 'V2', 'V3']
    noddi_metrics = ['dir', 'ICVF', 'ISOVF', 'OD']
    
    # Get ROI definitions
    jhu_rois = get_jhu_rois()
    
    # Process DTI metrics
    dti_dir = mri_dir / 'DTI'
    if dti_dir.exists():
        for metric in dti_metrics:
            # Look for normalized (template space) files with 'wr' prefix
            metric_files = list(dti_dir.glob(f'wr*_{metric}.nii')) + list(dti_dir.glob(f'wr*_ECC_{metric}.nii'))
            
            if metric_files:
                try:
                    data, _, _ = load_nifti_data(str(metric_files[0]))
                    if data is not None:
                        # Extract JHU ROI values
                        if jhu_atlas_data is not None:
                            for roi_name, roi_id in jhu_rois.items():
                                roi_values = extract_roi_mean_nonzero(data, jhu_atlas_data, roi_id)
                                results[f'DTI_{metric}_{roi_name}_mean'] = roi_values['mean_nonzero']
                                results[f'DTI_{metric}_{roi_name}_n_nonzero'] = roi_values['n_nonzero_voxels']
                        
                        # Extract hypothalamus ROI values
                        if hyp_atlas_data is not None:
                            # Get all unique ROI IDs in hypothalamus atlas
                            unique_hyp_ids = np.unique(hyp_atlas_data)
                            unique_hyp_ids = unique_hyp_ids[unique_hyp_ids > 0]
                            
                            for roi_id in unique_hyp_ids:
                                roi_id_int = int(roi_id)
                                roi_values = extract_roi_mean_nonzero(data, hyp_atlas_data, roi_id)
                                results[f'DTI_{metric}_hyp_roi{roi_id_int}_mean'] = roi_values['mean_nonzero']
                                results[f'DTI_{metric}_hyp_roi{roi_id_int}_n_nonzero'] = roi_values['n_nonzero_voxels']
                except Exception as e:
                    results[f'DTI_{metric}_error'] = str(e)
    
    # Process NODDI metrics
    noddi_dir = mri_dir / 'NODDI'
    if noddi_dir.exists():
        for metric in noddi_metrics:
            # Look for normalized (template space) files with 'wr' prefix
            # NODDI files might use different naming conventions
            if metric == 'dir':
                metric_files = list(noddi_dir.glob(f'wr*NODDI_dir.nii'))
            elif metric == 'ICVF':
                metric_files = list(noddi_dir.glob(f'wr*NODDI_ICVF.nii')) + list(noddi_dir.glob(f'wr*NODDI_ficvf.nii'))
            elif metric == 'ISOVF':
                metric_files = list(noddi_dir.glob(f'wr*NODDI_ISOVF.nii')) + list(noddi_dir.glob(f'wr*NODDI_fiso.nii'))
            elif metric == 'OD':
                metric_files = list(noddi_dir.glob(f'wr*NODDI_OD.nii')) + list(noddi_dir.glob(f'wr*NODDI_odi.nii'))
            
            if metric_files:
                try:
                    data, _, _ = load_nifti_data(str(metric_files[0]))
                    if data is not None:
                        # Extract JHU ROI values
                        if jhu_atlas_data is not None:
                            for roi_name, roi_id in jhu_rois.items():
                                roi_values = extract_roi_mean_nonzero(data, jhu_atlas_data, roi_id)
                                results[f'NODDI_{metric}_{roi_name}_mean'] = roi_values['mean_nonzero']
                                results[f'NODDI_{metric}_{roi_name}_n_nonzero'] = roi_values['n_nonzero_voxels']
                        
                        # Extract hypothalamus ROI values
                        if hyp_atlas_data is not None:
                            unique_hyp_ids = np.unique(hyp_atlas_data)
                            unique_hyp_ids = unique_hyp_ids[unique_hyp_ids > 0]
                            
                            for roi_id in unique_hyp_ids:
                                roi_id_int = int(roi_id)
                                roi_values = extract_roi_mean_nonzero(data, hyp_atlas_data, roi_id)
                                results[f'NODDI_{metric}_hyp_roi{roi_id_int}_mean'] = roi_values['mean_nonzero']
                                results[f'NODDI_{metric}_hyp_roi{roi_id_int}_n_nonzero'] = roi_values['n_nonzero_voxels']
                except Exception as e:
                    results[f'NODDI_{metric}_error'] = str(e)
    
    return results


def process_subject(subject_dir, hyp_atlas_data, jhu_atlas_data):
    """Process a single subject's data."""
    subject_dir = Path(subject_dir)
    
    # Initialize results
    results = {
        'subject_id': subject_dir.name,
        'processing_status': 'started'
    }
    
    try:
        # Validate subject
        is_valid, validation_msg = validate_subject(subject_dir)
        if not is_valid:
            results['processing_status'] = f'skipped: {validation_msg}'
            return results
        
        # Process volume measurements
        try:
            volume_results = process_volumes(subject_dir)
            results.update(volume_results)
            results['volumes_status'] = 'success'
        except Exception as e:
            results['volumes_status'] = f'failed: {str(e)}'
        
        # Process template-space diffusion data
        try:
            diffusion_results = process_template_diffusion(subject_dir, hyp_atlas_data, jhu_atlas_data)
            results.update(diffusion_results)
            results['diffusion_status'] = 'success'
        except Exception as e:
            results['diffusion_status'] = f'failed: {str(e)}'
        
        results['processing_status'] = 'completed'
        
    except Exception as e:
        results['processing_status'] = f'failed: {str(e)}'
    
    return results


def main():
    parser = argparse.ArgumentParser(
        description='Extract volumes and template-space diffusion values from CAT12 data'
    )
    parser.add_argument('input_path', help='Path to subject directory or parent directory with multiple subjects')
    parser.add_argument('-o', '--output', default='cat12_volumes_diffusion.csv', help='Output CSV filename')
    parser.add_argument('--batch', action='store_true', help='Process multiple subjects in batch mode')
    parser.add_argument('--hyp-atlas', default=HYPOTHALAMUS_ATLAS_PATH, 
                       help='Path to hypothalamus atlas template')
    parser.add_argument('--jhu-atlas', default=JHU_ATLAS_PATH,
                       help='Path to JHU atlas template')
    parser.add_argument('--verbose', action='store_true', help='Print detailed processing information')
    
    args = parser.parse_args()
    
    # Load template atlases
    print("Loading template atlases...")
    hyp_atlas_data, jhu_atlas_data = None, None
    
    if Path(args.hyp_atlas).exists():
        hyp_atlas_data, _, _ = load_nifti_data(args.hyp_atlas)
        if hyp_atlas_data is not None:
            print(f"✓ Loaded hypothalamus atlas: {args.hyp_atlas}")
            print(f"  Unique ROIs: {np.unique(hyp_atlas_data[hyp_atlas_data > 0])}")
        else:
            print(f"✗ Failed to load hypothalamus atlas")
    else:
        print(f"✗ Hypothalamus atlas not found: {args.hyp_atlas}")
    
    if Path(args.jhu_atlas).exists():
        jhu_atlas_data, _, _ = load_nifti_data(args.jhu_atlas)
        if jhu_atlas_data is not None:
            print(f"✓ Loaded JHU atlas: {args.jhu_atlas}")
            print(f"  Unique ROIs: {len(np.unique(jhu_atlas_data[jhu_atlas_data > 0]))}")
        else:
            print(f"✗ Failed to load JHU atlas")
    else:
        print(f"✗ JHU atlas not found: {args.jhu_atlas}")
    
    if hyp_atlas_data is None and jhu_atlas_data is None:
        print("\nWARNING: No atlases loaded. Diffusion analysis will be skipped.")
        print("Continuing with volume measurements only...\n")
    
    input_path = Path(args.input_path)
    all_results = []
    
    processing_stats = {
        'total_folders': 0,
        'valid_subjects': 0,
        'completed': 0,
        'skipped': 0,
        'failed': 0
    }
    
    if args.batch:
        # Process multiple subjects
        subject_dirs = [d for d in input_path.iterdir() if d.is_dir()]
        processing_stats['total_folders'] = len(subject_dirs)
        
        print(f"\nStarting batch processing of {len(subject_dirs)} folders...")
        print("=" * 60)
        
        for i, subject_dir in enumerate(subject_dirs, 1):
            try:
                # Quick pre-check for CAT12 output
                if not (subject_dir / 'nifti' / 'cat12_v2560').exists():
                    processing_stats['skipped'] += 1
                    if args.verbose:
                        print(f"[{i:4d}/{len(subject_dirs):4d}] SKIP: {subject_dir.name} - No CAT12 output")
                    continue
                
                processing_stats['valid_subjects'] += 1
                
                if args.verbose:
                    print(f"[{i:4d}/{len(subject_dirs):4d}] Processing: {subject_dir.name}")
                
                results = process_subject(subject_dir, hyp_atlas_data, jhu_atlas_data)
                all_results.append(results)
                
                # Update statistics
                status = results.get('processing_status', 'unknown')
                if status == 'completed':
                    processing_stats['completed'] += 1
                elif status.startswith('failed'):
                    processing_stats['failed'] += 1
                else:
                    processing_stats['skipped'] += 1
                
                if args.verbose:
                    print(f"[{i:4d}/{len(subject_dirs):4d}] Result: {status}")
                
                # Progress update every 50 subjects
                if not args.verbose and i % 50 == 0:
                    print(f"Processed {i}/{len(subject_dirs)} folders...")
                    
            except Exception as e:
                processing_stats['failed'] += 1
                if args.verbose:
                    print(f"[{i:4d}/{len(subject_dirs):4d}] ERROR: {subject_dir.name} - {str(e)}")
                all_results.append({
                    'subject_id': subject_dir.name,
                    'processing_status': f'exception: {str(e)}'
                })
    else:
        # Process single subject
        processing_stats['total_folders'] = 1
        print(f"Processing subject: {input_path.name}")
        
        results = process_subject(input_path, hyp_atlas_data, jhu_atlas_data)
        all_results.append(results)
        processing_stats['valid_subjects'] = 1
        if results.get('processing_status') == 'completed':
            processing_stats['completed'] = 1
    
    # Save results to CSV
    if all_results:
        df = pd.DataFrame(all_results)
        df.to_csv(args.output, index=False)
        
        # Print summary
        print("\n" + "=" * 60)
        print("PROCESSING SUMMARY")
        print("=" * 60)
        print(f"Total folders examined: {processing_stats['total_folders']}")
        print(f"Valid CAT12 subjects: {processing_stats['valid_subjects']}")
        print(f"Successfully completed: {processing_stats['completed']}")
        print(f"Skipped: {processing_stats['skipped']}")
        print(f"Failed/errors: {processing_stats['failed']}")
        print(f"\nResults saved to: {args.output}")
        print(f"Total rows in output: {len(df)}")
        print(f"Total columns: {len(df.columns)}")
        
        # Show completion rate
        if processing_stats['valid_subjects'] > 0:
            completion_rate = (processing_stats['completed'] / processing_stats['valid_subjects']) * 100
            print(f"Success rate: {completion_rate:.1f}%")
        
        # Show sample columns
        print(f"\nSample columns in output:")
        for col in list(df.columns)[:20]:
            print(f"  - {col}")
        if len(df.columns) > 20:
            print(f"  ... and {len(df.columns) - 20} more columns")
    else:
        print("No subjects found for processing")
        print(f"Total folders examined: {processing_stats['total_folders']}")


if __name__ == "__main__":
    main()