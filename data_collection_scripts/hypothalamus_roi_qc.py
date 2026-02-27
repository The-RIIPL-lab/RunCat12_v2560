#!/usr/bin/env python3
"""
QC overlay script: hypothalamus ROIs on native-space T1w images.

Produces per-subject PNG mosaics:
  Top row    : 4 coronal slices evenly spaced across the full ROI extent
  Bottom row : sagittal slice through the midpoint of ROI 1 (left panel)
               sagittal slice through the midpoint of ROI 2 (right panel)

ROI boundaries are drawn as thin coloured contour outlines, not opaque overlays.

Inputs (inside each subject's CAT12 output folder):
  mri/m*.nii             bias-corrected native-space T1
  whypothalamusAtlas.nii native-space hypothalamus atlas (ROI 1 and 2)

Usage:
  Single subject : python hypothalamus_roi_qc.py /path/to/SUBJECT_ID
  Batch          : python hypothalamus_roi_qc.py --batch /path/to/project/
"""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use('Agg')           # headless — safe for HPC nodes without a display
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import nibabel as nib
import numpy as np
from matplotlib.lines import Line2D


ROI_COLORS = {1: '#FF6B6B', 2: '#4ECDC4'}   # coral-red / teal
ROI_LABELS = {1: 'ROI 1',    2: 'ROI 2'}
CONTOUR_LW = 1.5


# ---------------------------------------------------------------------------
# Path helpers (mirrors extract_roi_values.py)
# ---------------------------------------------------------------------------

def find_cat12_dir(subject_dir):
    """Locate the CAT12 output folder, handling different build numbers."""
    nifti_dir  = Path(subject_dir) / 'nifti'
    candidates = sorted(nifti_dir.glob('cat12_v*'))
    return candidates[0] if candidates else nifti_dir / 'cat12_v2560'


def find_native_t1(mri_dir):
    """Return the bias-corrected native-space T1 (m*.nii, excluding mwp* files)."""
    candidates = [f for f in Path(mri_dir).glob('m*.nii')
                  if not f.name.startswith('mw')]
    return candidates[0] if candidates else None


# ---------------------------------------------------------------------------
# Image utilities
# ---------------------------------------------------------------------------

def load_canonical(filepath):
    """Load a NIfTI volume reoriented to RAS+ so axes are ordered (R, A, S)."""
    return nib.as_closest_canonical(nib.load(str(filepath))).get_fdata()


def get_slice(data, axis, idx):
    """
    Extract a 2D display slice from a RAS-oriented volume.

    Returns an array ready for imshow (default origin='upper') with:
      row 0  = superior (top of image)
      col 0  = patient left  for coronal slices  (axis 1)
      col 0  = posterior     for sagittal slices (axis 0)
    """
    if axis == 0:
        sl = data[idx, :, :]    # sagittal  → shape (A, S)
    elif axis == 1:
        sl = data[:, idx, :]    # coronal   → shape (R, S)
    else:
        sl = data[:, :, idx]    # axial     → shape (R, A)
    return np.flipud(sl.T)      # → (S, *) with S_max at row 0


def roi_extent(atlas, roi_ids, axis):
    """First and last slice indices along *axis* containing any of *roi_ids*."""
    present = np.isin(atlas, roi_ids)
    other   = tuple(a for a in range(3) if a != axis)
    indices = np.where(present.any(axis=other))[0]
    if len(indices) == 0:
        return None, None
    return int(indices[0]), int(indices[-1])


def roi_midpoint(atlas, roi_id, axis):
    """Index along *axis* at the centre of *roi_id*'s spatial extent."""
    other   = tuple(a for a in range(3) if a != axis)
    indices = np.where((atlas == roi_id).any(axis=other))[0]
    if len(indices) == 0:
        return None
    return int(indices[len(indices) // 2])


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def draw_panel(ax, t1_sl, atl_sl, vmin, vmax):
    """
    Render *t1_sl* in grey and overlay each ROI as a thin contour outline.

    Both *t1_sl* and *atl_sl* must already be in display coordinates (i.e.
    the output of get_slice), so that imshow and contour share the same
    coordinate system.
    """
    ax.imshow(t1_sl, cmap='gray', vmin=vmin, vmax=vmax,
              interpolation='bilinear', aspect='equal')

    for roi_id, color in ROI_COLORS.items():
        mask = (atl_sl == roi_id).astype(np.float32)
        if mask.max() > 0:
            ax.contour(mask, levels=[0.5], colors=[color],
                       linewidths=[CONTOUR_LW], zorder=5)

    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)


# ---------------------------------------------------------------------------
# Main per-subject function
# ---------------------------------------------------------------------------

def generate_qc_mosaic(subject_dir, output_path):
    """
    Build and save the QC mosaic for one subject.
    Returns True on success, False if a required file is missing or the ROI
    is not found.
    """
    subject_dir = Path(subject_dir)
    cat12_dir   = find_cat12_dir(subject_dir)
    mri_dir     = cat12_dir / 'mri'

    atlas_path = cat12_dir / 'whypothalamusAtlas.nii'
    t1_path    = find_native_t1(cat12_dir)

    if not atlas_path.exists():
        print(f"  SKIP {subject_dir.name}: atlas not found ({atlas_path})")
        return False
    if t1_path is None:
        print(f"  SKIP {subject_dir.name}: native T1 not found in {cat12_dir}")
        return False

    t1_data  = load_canonical(t1_path)
    atl_data = load_canonical(atlas_path)

    if t1_data.shape != atl_data.shape:
        print(f"  SKIP {subject_dir.name}: shape mismatch "
              f"T1={t1_data.shape} atlas={atl_data.shape}")
        return False

    # Global intensity window from non-zero brain voxels
    brain = t1_data[t1_data > 0]
    vmin  = float(np.percentile(brain, 1.0))
    vmax  = float(np.percentile(brain, 99.5))

    # ---- Coronal slice selection (axis 1 = A direction in RAS) -------------
    c_start, c_end = roi_extent(atl_data, [1, 2], axis=1)
    if c_start is None:
        print(f"  SKIP {subject_dir.name}: no ROI voxels found in atlas")
        return False

    # 4 evenly-spaced slices from first to last coronal extent
    cor_indices = [int(c_start + round((c_end - c_start) * i / 3))
                   for i in range(4)]

    # ---- Sagittal midpoints (axis 0 = R direction in RAS) ------------------
    sag_idx = {roi_id: roi_midpoint(atl_data, roi_id, axis=0)
               for roi_id in (1, 2)}

    # ---- Figure layout ------------------------------------------------------
    fig = plt.figure(figsize=(16, 8), facecolor='black')
    gs  = gridspec.GridSpec(2, 4, figure=fig,
                            hspace=0.06, wspace=0.04,
                            left=0.02, right=0.98, top=0.92, bottom=0.08)

    # Top row: 4 coronal slices
    for col, c_idx in enumerate(cor_indices):
        ax = fig.add_subplot(gs[0, col])
        ax.set_facecolor('black')
        draw_panel(ax,
                   get_slice(t1_data,  axis=1, idx=c_idx),
                   get_slice(atl_data, axis=1, idx=c_idx),
                   vmin, vmax)
        ax.set_title(f'Coronal  {c_idx}', color='white', fontsize=8, pad=3)

    # Bottom row: sagittal slices, each spanning 2 grid columns
    for grid_col, roi_id in enumerate((1, 2)):
        ax    = fig.add_subplot(gs[1, grid_col * 2 : grid_col * 2 + 2])
        ax.set_facecolor('black')
        idx   = sag_idx[roi_id]
        label = f'{ROI_LABELS[roi_id]} — sagittal midpoint'

        if idx is not None:
            draw_panel(ax,
                       get_slice(t1_data,  axis=0, idx=idx),
                       get_slice(atl_data, axis=0, idx=idx),
                       vmin, vmax)
            ax.set_title(f'{label}  (slice {idx})',
                         color='white', fontsize=8, pad=3)
        else:
            ax.text(0.5, 0.5, f'{label}\n(ROI not found)',
                    ha='center', va='center', color='white',
                    transform=ax.transAxes, fontsize=9)
            ax.set_xticks([])
            ax.set_yticks([])

    # Legend
    fig.legend(
        handles=[
            Line2D([0], [0], color=ROI_COLORS[r], lw=CONTOUR_LW,
                   label=ROI_LABELS[r])
            for r in (1, 2)
        ],
        loc='lower center', ncol=2, fontsize=9,
        frameon=False, labelcolor='white',
        bbox_to_anchor=(0.5, 0.01),
    )

    fig.suptitle(f'Hypothalamus ROI QC — {subject_dir.name}',
                 color='white', fontsize=11, y=0.97)

    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(str(output_path), dpi=200, bbox_inches='tight', facecolor='black')
    plt.close(fig)
    print(f"  Saved: {output_path}")
    return True


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Generate hypothalamus ROI QC mosaic images'
    )
    parser.add_argument(
        'input_path',
        help='Subject directory, or parent directory when --batch is used',
    )
    parser.add_argument(
        '-o', '--output-dir', default=None,
        help='Directory to write QC images into '
             '(default: <cat12_dir>/QC_registration/ inside each subject)',
    )
    parser.add_argument(
        '--batch', action='store_true',
        help='Process all subdirectories of input_path',
    )
    args = parser.parse_args()

    input_path = Path(args.input_path)

    if args.batch:
        subject_dirs = sorted(d for d in input_path.iterdir() if d.is_dir())
        print(f"Batch mode: {len(subject_dirs)} folders found")
        ok = skip = errors = 0

        for sdir in subject_dirs:
            cat12_dir = find_cat12_dir(sdir)
            if not cat12_dir.exists():
                skip += 1
                continue
            out_dir     = Path(args.output_dir) if args.output_dir \
                          else cat12_dir / 'QC_registration'
            output_path = out_dir / f'hypothalamus_roi_qc_{sdir.name}.png'
            try:
                if generate_qc_mosaic(sdir, output_path):
                    ok += 1
                else:
                    skip += 1
            except Exception as e:
                print(f"  ERROR {sdir.name}: {e}")
                errors += 1

        print(f"\nDone — completed: {ok}, skipped: {skip}, errors: {errors}")

    else:
        cat12_dir   = find_cat12_dir(input_path)
        out_dir     = Path(args.output_dir) if args.output_dir \
                      else cat12_dir / 'QC_registration'
        output_path = out_dir / f'hypothalamus_roi_qc_{input_path.name}.png'
        generate_qc_mosaic(input_path, output_path)


if __name__ == '__main__':
    main()
