#!/usr/bin/env python3
"""
Post-processing script for CAT12 workflow
Extracts volumes and diffusion metrics from hypothalamus ROIs

This script processes CAT12 outputs to extract:
1. Native space volumes (brain mask, hypothalamus, tissues)
2. Native space hypothalamus diffusion metrics (DTI, NODDI)
3. Normalized space hypothalamus diffusion metrics (DTI, NODDI)

Author: Simplified for hypothalamus-focused analysis
Date: December 2025
"""

import numpy as np
import nibabel as nib
import pandas as pd
from pathlib import Path
import argparse


# Template atlas path for normalized diffusion analysis
# UPDATE THIS
HYPOTHALAMUS_ATLAS_PATH = "/isilon/datalake/riipl/original/DEMONco/Hellcat-12.9/resampled_template_ROIs/rhypothalamusAtlas_template_v2.nii"


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
        voxel_dims = header.get_zooms()[:3]
        voxel_volume_mm3 = np.prod(voxel_dims)
        n_voxels = np.sum(data > 0)
        volume_mm3 = n_voxels * voxel_volume_mm3
        return volume_mm3
    except Exception as e:
        print(f"Error calculating volume: {e}")
        return np.nan


def calculate_roi_volume(data, header, roi_id):
    """Calculate volume for specific ROI ID."""
    try:
        voxel_dims = header.get_zooms()[:3]
        voxel_volume_mm3 = np.prod(voxel_dims)
        n_voxels = np.sum(data == roi_id)
        volume_mm3 = n_voxels * voxel_volume_mm3
        return volume_mm3
    except Exception as e:
        print(f"Error calculating ROI volume: {e}")
        return np.nan


def extract_roi_mean_nonzero(data_img, atlas_img, roi_id):
    """Extract mean of non-zero values from ROI."""
    try:
        mask_data = atlas_img == roi_id
        if np.sum(mask_data) == 0:
            return np.nan
        
        roi_values = data_img[mask_data]
        nonzero_values = roi_values[roi_values != 0]
        nonzero_values = nonzero_values[~np.isnan(nonzero_values)]
        
        if len(nonzero_values) == 0:
            return np.nan
        
        return float(np.mean(nonzero_values))
    except Exception as e:
        print(f"Error extracting ROI values: {e}")
        return np.nan


def validate_subject(subject_dir):
    """Validate that a subject has the required CAT12 outputs."""
    subject_dir = Path(subject_dir)
    cat12_dir = subject_dir / 'nifti' / 'cat12_v2560'
    if not cat12_dir.exists():
        return False, "CAT12 output directory not found"
    
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
                volume_mm3 = calculate_volume_mm3(data, header)
                results['brain_mask_volume_mm3'] = volume_mm3
        except Exception as e:
            results['brain_mask_volume_mm3'] = np.nan
    else:
        results['brain_mask_volume_mm3'] = np.nan
    
    # 2. Hypothalamus volumes (ROI 1 and 2)
    hyp_atlas_file = cat12_dir / 'whypothalamusAtlas.nii'
    if hyp_atlas_file.exists():
        try:
            data, affine, header = load_nifti_data(str(hyp_atlas_file))
            if data is not None:
                results['hypothalamus_roi1_volume_mm3'] = calculate_roi_volume(data, header, 1)
                results['hypothalamus_roi2_volume_mm3'] = calculate_roi_volume(data, header, 2)
        except Exception as e:
            results['hypothalamus_roi1_volume_mm3'] = np.nan
            results['hypothalamus_roi2_volume_mm3'] = np.nan
    else:
        results['hypothalamus_roi1_volume_mm3'] = np.nan
        results['hypothalamus_roi2_volume_mm3'] = np.nan
    
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
            tissue_files = list(mri_dir.glob(f'{prefix}*.nii'))
        
        if tissue_files:
            try:
                data, affine, header = load_nifti_data(str(tissue_files[0]))
                if data is not None:
                    # For tissue probability maps, threshold at 0.5
                    binary_data = (data > 0.5).astype(int)
                    volume_mm3 = calculate_volume_mm3(binary_data, header)
                    results[f'{tissue_name}_volume_mm3'] = volume_mm3
            except Exception as e:
                results[f'{tissue_name}_volume_mm3'] = np.nan
        else:
            results[f'{tissue_name}_volume_mm3'] = np.nan
    
    return results


def process_native_diffusion(subject_dir):
    """Process native-space diffusion data from hypothalamus ROIs."""
    results = {}
    cat12_dir = subject_dir / 'nifti' / 'cat12_v2560'
    mri_dir = cat12_dir / 'mri'
    
    # Load native hypothalamus atlas
    hyp_atlas_file = cat12_dir / 'whypothalamusAtlas.nii'
    if not hyp_atlas_file.exists():
        return results
    
    hyp_atlas_data, _, _ = load_nifti_data(str(hyp_atlas_file))
    if hyp_atlas_data is None:
        return results
    
    # Define diffusion metrics
    dti_metrics = ['FA', 'MD', 'L1', 'L2', 'L3']
    noddi_metrics = ['ICVF', 'ISOVF', 'OD']
    
    # Process DTI metrics
    dti_dir = mri_dir / 'DTI'
    if dti_dir.exists():
        for metric in dti_metrics:
            # Look for native space files (without 'w' or 'wr' prefix)
            metric_files = list(dti_dir.glob(f'*_{metric}.nii'))
            # Filter out warped files
            metric_files = [f for f in metric_files if not f.name.startswith('w')]
            
            if metric_files:
                try:
                    data, _, _ = load_nifti_data(str(metric_files[0]))
                    if data is not None:
                        # Extract for hypothalamus ROI 1 and 2
                        for roi_id in [1, 2]:
                            mean_val = extract_roi_mean_nonzero(data, hyp_atlas_data, roi_id)
                            results[f'native_DTI_{metric}_hyp_roi{roi_id}_mean'] = mean_val
                except Exception as e:
                    for roi_id in [1, 2]:
                        results[f'native_DTI_{metric}_hyp_roi{roi_id}_mean'] = np.nan
    
    # Process NODDI metrics
    noddi_dir = mri_dir / 'NODDI'
    if noddi_dir.exists():
        for metric in noddi_metrics:
            # Look for native space files
            if metric == 'ICVF':
                metric_files = list(noddi_dir.glob('*NODDI_ICVF.nii')) + list(noddi_dir.glob('*NODDI_ficvf.nii'))
            elif metric == 'ISOVF':
                metric_files = list(noddi_dir.glob('*NODDI_ISOVF.nii')) + list(noddi_dir.glob('*NODDI_fiso.nii'))
            elif metric == 'OD':
                metric_files = list(noddi_dir.glob('*NODDI_OD.nii')) + list(noddi_dir.glob('*NODDI_odi.nii'))
            
            # Filter out warped files
            metric_files = [f for f in metric_files if not f.name.startswith('w')]
            
            if metric_files:
                try:
                    data, _, _ = load_nifti_data(str(metric_files[0]))
                    if data is not None:
                        # Extract for hypothalamus ROI 1 and 2
                        for roi_id in [1, 2]:
                            mean_val = extract_roi_mean_nonzero(data, hyp_atlas_data, roi_id)
                            results[f'native_NODDI_{metric}_hyp_roi{roi_id}_mean'] = mean_val
                except Exception as e:
                    for roi_id in [1, 2]:
                        results[f'native_NODDI_{metric}_hyp_roi{roi_id}_mean'] = np.nan
    
    return results


def process_normalized_diffusion(subject_dir, template_hyp_atlas_data):
    """Process normalized/template-space diffusion data from hypothalamus ROIs."""
    results = {}
    
    if template_hyp_atlas_data is None:
        return results
    
    cat12_dir = subject_dir / 'nifti' / 'cat12_v2560'
    mri_dir = cat12_dir / 'mri'
    
    # Define diffusion metrics
    dti_metrics = ['FA', 'MD', 'L1', 'L2', 'L3']
    noddi_metrics = ['ICVF', 'ISOVF', 'OD']
    
    # Process DTI metrics
    dti_dir = mri_dir / 'DTI'
    if dti_dir.exists():
        for metric in dti_metrics:
            # Look for normalized files with 'wr' prefix
            metric_files = list(dti_dir.glob(f'wr*_{metric}.nii')) + list(dti_dir.glob(f'wr*_ECC_{metric}.nii'))
            
            if metric_files:
                try:
                    data, _, _ = load_nifti_data(str(metric_files[0]))
                    if data is not None:
                        # Extract for hypothalamus ROI 1 and 2
                        for roi_id in [1, 2]:
                            mean_val = extract_roi_mean_nonzero(data, template_hyp_atlas_data, roi_id)
                            results[f'normalized_DTI_{metric}_hyp_roi{roi_id}_mean'] = mean_val
                except Exception as e:
                    for roi_id in [1, 2]:
                        results[f'normalized_DTI_{metric}_hyp_roi{roi_id}_mean'] = np.nan
    
    # Process NODDI metrics
    noddi_dir = mri_dir / 'NODDI'
    if noddi_dir.exists():
        for metric in noddi_metrics:
            # Look for normalized files with 'wr' prefix
            if metric == 'ICVF':
                metric_files = list(noddi_dir.glob('wr*NODDI_ICVF.nii')) + list(noddi_dir.glob('wr*NODDI_ficvf.nii'))
            elif metric == 'ISOVF':
                metric_files = list(noddi_dir.glob('wr*NODDI_ISOVF.nii')) + list(noddi_dir.glob('wr*NODDI_fiso.nii'))
            elif metric == 'OD':
                metric_files = list(noddi_dir.glob('wr*NODDI_OD.nii')) + list(noddi_dir.glob('wr*NODDI_odi.nii'))
            
            if metric_files:
                try:
                    data, _, _ = load_nifti_data(str(metric_files[0]))
                    if data is not None:
                        # Extract for hypothalamus ROI 1 and 2
                        for roi_id in [1, 2]:
                            mean_val = extract_roi_mean_nonzero(data, template_hyp_atlas_data, roi_id)
                            results[f'normalized_NODDI_{metric}_hyp_roi{roi_id}_mean'] = mean_val
                except Exception as e:
                    for roi_id in [1, 2]:
                        results[f'normalized_NODDI_{metric}_hyp_roi{roi_id}_mean'] = np.nan
    
    return results


def process_subject(subject_dir, template_hyp_atlas_data):
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
        except Exception as e:
            pass
        
        # Process native-space diffusion
        try:
            native_diffusion_results = process_native_diffusion(subject_dir)
            results.update(native_diffusion_results)
        except Exception as e:
            pass
        
        # Process normalized-space diffusion
        try:
            normalized_diffusion_results = process_normalized_diffusion(subject_dir, template_hyp_atlas_data)
            results.update(normalized_diffusion_results)
        except Exception as e:
            pass
        
        results['processing_status'] = 'completed'
        
    except Exception as e:
        results['processing_status'] = f'failed: {str(e)}'
    
    return results


def main():
    parser = argparse.ArgumentParser(
        description='Extract hypothalamus volumes and diffusion metrics from CAT12 data'
    )
    parser.add_argument('input_path', help='Path to subject directory or parent directory with multiple subjects')
    parser.add_argument('-o', '--output', default='cat12_hypothalamus_data.csv', help='Output CSV filename')
    parser.add_argument('--batch', action='store_true', help='Process multiple subjects in batch mode')
    parser.add_argument('--hyp-atlas', default=HYPOTHALAMUS_ATLAS_PATH, 
                       help='Path to hypothalamus atlas template for normalized analysis')
    parser.add_argument('--verbose', action='store_true', help='Print detailed processing information')
    
    args = parser.parse_args()
    
    # Load template hypothalamus atlas for normalized diffusion analysis
    print("Loading template hypothalamus atlas...")
    template_hyp_atlas_data = None
    
    if Path(args.hyp_atlas).exists():
        template_hyp_atlas_data, _, _ = load_nifti_data(args.hyp_atlas)
        if template_hyp_atlas_data is not None:
            unique_rois = np.unique(template_hyp_atlas_data[template_hyp_atlas_data > 0])
            print(f"✓ Loaded hypothalamus atlas: {args.hyp_atlas}")
            print(f"  Unique ROIs: {unique_rois}")
        else:
            print(f"✗ Failed to load hypothalamus atlas")
    else:
        print(f"✗ Hypothalamus atlas not found: {args.hyp_atlas}")
        print("WARNING: Normalized diffusion analysis will be skipped.\n")
    
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
                
                results = process_subject(subject_dir, template_hyp_atlas_data)
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
        
        results = process_subject(input_path, template_hyp_atlas_data)
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
        
        # Show column groups
        print(f"\nOutput columns:")
        print(f"  - subject_id, processing_status")
        print(f"  - brain_mask_volume_mm3")
        print(f"  - hypothalamus_roi1_volume_mm3, hypothalamus_roi2_volume_mm3")
        print(f"  - gm_volume_mm3, wm_volume_mm3, csf_volume_mm3")
        print(f"  - native_DTI_* and native_NODDI_* metrics (hyp_roi1 and hyp_roi2)")
        print(f"  - normalized_DTI_* and normalized_NODDI_* metrics (hyp_roi1 and hyp_roi2)")
    else:
        print("No subjects found for processing")


if __name__ == "__main__":
    main()