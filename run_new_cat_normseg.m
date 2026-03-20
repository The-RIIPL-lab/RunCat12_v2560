function run_new_cat_normseg(base_dir, output_root, varargin)
    if nargin < 2
        output_root = '';
    end
    % NATIVE-SPACE LABELMAP EXPORT: This script creates native-space labelmaps
    % in Step 5 using the process_atlases_to_native() function.
    % Output files will have 'w' prefix (e.g., wneuromorphometrics.nii)
    %
    % Optional name-value arguments (defaults shown):
    %   'hypothalamus'  true   — include hypothalamusAtlas in warp + atlas list
    %   'surface'       false  — run CAT12 surface mapping
    %   'qc_images'     false  — generate SPM/CAT12 registration QC images
    %
    % Examples:
    %   run_new_cat_normseg('/path/to/subj')
    %   run_new_cat_normseg('/path/to/subj', '', 'surface', true)
    %   run_new_cat_normseg('/path/to/subj', '', 'hypothalamus', false, 'qc_images', true)

    p = inputParser();
    addParameter(p, 'hypothalamus', true);
    addParameter(p, 'surface',      false);
    addParameter(p, 'qc_images',    false);
    parse(p, varargin{:});
    opts = p.Results;

    % Determine repo root from this file's location (works regardless of where repo is cloned)
    repo_root = fileparts(mfilename('fullpath'));
    spm_dir   = fullfile(repo_root, 'libs', 'spm12', 'spm12');
    addpath(spm_dir);
    addpath(repo_root);

    % Atlas list — hypothalamusAtlas included only when requested
    label_map_path  = fullfile(spm_dir, 'toolbox', 'cat12', 'templates_MNI152NLin2009cAsym', '');
    label_maps_files = {'neuromorphometrics.nii','aal3.nii','cobra.nii','hammers.nii','lpba40.nii','ibsr.nii','JHU.nii','brainmask_T1.nii'};
    if opts.hypothalamus
        label_maps_files{end+1} = 'hypothalamusAtlas.nii';
    end

    % Template path for QC
    template_path = fullfile(spm_dir, 'toolbox', 'cat12', 'templates_MNI152NLin2009cAsym', 'Template_0_GS.nii');

    % Initialize SPM
    spm('defaults', 'FMRI');
    spm_jobman('initcfg');

    % Detect CAT12 build number for output folder naming
    try
        [~, cat_build] = cat_version();
        cat_folder_name = sprintf('cat12_v%s', cat_build);
    catch
        cat_folder_name = 'cat12_v2560';  % fallback if detection fails
    end

    % Determine output root: separate dir or same as input
    if isempty(output_root)
        newdir = fullfile(base_dir, 'nifti', cat_folder_name);
    else
        % Extract session name robustly (handles trailing slashes)
        base_dir_clean = regexprep(base_dir, '[/\\]+$', '');
        [~, session_name] = fileparts(base_dir_clean);
        newdir = fullfile(output_root, session_name, 'nifti', cat_folder_name);
    end
    if ~exist(newdir, 'dir')
        mkdir(newdir);
    end

    % Create QC directory only when QC images are requested
    if opts.qc_images
        qc_dir = fullfile(newdir, 'QC_registration');
        if ~exist(qc_dir, 'dir')
            mkdir(qc_dir);
        end
    else
        qc_dir = '';
    end

% Find all T1w images across all projects
t1wfiles = [];

% MARVEL
marvel_files = dir(fullfile([base_dir, '/nifti/'], 'M*-tfl3d*ns.nii'));
if ~isempty(marvel_files)
	    warning('MARVEL T1w image found');
	        t1wfiles = [t1wfiles; marvel_files];
end

% ADRC
adrc_files = dir(fullfile([base_dir, '/nifti/'], '3*-tfl3d116ns.nii'));
if ~isempty(adrc_files)
	    warning('ADRC T1w image found');
	        t1wfiles = [t1wfiles; adrc_files];
end

% SWITCH
switch_files = dir(fullfile([base_dir, '/nifti/'], 'S*-tfl3d*ns.nii'));
if ~isempty(switch_files)
	    warning('SWITCH T1w image found');
	        t1wfiles = [t1wfiles; switch_files];
end

% Error if no files found at all
if isempty(t1wfiles)
	    error('No T1w images found')
end

    % Copy the T1w image to the new directory
    copyfile(fullfile(t1wfiles(1).folder, t1wfiles(1).name), newdir, 'f');

    % Check if processing has already been done
    mri_dir = fullfile(newdir, 'mri');
    y_files = dir(fullfile(mri_dir, 'y_*-tfl3d116ns.nii'));
    iy_files = dir(fullfile(mri_dir, 'iy_*-tfl3d116ns.nii'));

    disp('========================================');
    disp(sprintf('Subject dir: %s', base_dir));
    disp(sprintf('Output folder: %s', cat_folder_name));
    disp(sprintf('Output dir: %s', newdir));
    disp(sprintf('MRI dir: %s', mri_dir));
    disp(sprintf('Input T1: %s', fullfile(newdir, t1wfiles(end).name)));
    if opts.qc_images
        disp(sprintf('QC dir: %s', qc_dir));
    end
    disp(sprintf('Options: hypothalamus=%d  surface=%d  qc_images=%d', ...
        opts.hypothalamus, opts.surface, opts.qc_images));
    disp('========================================');

    % Get DTI, NODDI, and ASL files
    [asl_struct, asl_success] = getASLFiles(base_dir);
    [dti_struct, dti_success] = getDTIFiles(base_dir);
    [noddi_struct, noddi_success] = getNODDIFiles(base_dir);

    % Summary of available data
    disp('Available data:');
    disp(sprintf('  DTI: %s', mat2str(dti_success)));
    disp(sprintf('  NODDI: %s', mat2str(noddi_success)));
    disp(sprintf('  ASL: %s', mat2str(asl_success)));
    disp('========================================');

    % Step 1: CAT12 Segmentation and Normalization (if not already done)
    if isempty(y_files)
        disp('=== STEP 1: Running CAT12 segmentation and normalization ===');
        run_cat12_segmentation(newdir, t1wfiles, spm_dir, opts.surface);

        % Re-check for deformation fields
        y_files = dir(fullfile(mri_dir, 'y_*-tfl3d*.nii'));
        iy_files = dir(fullfile(mri_dir, 'iy_*-tfl3d*ns.nii'));

        if isempty(y_files) || isempty(iy_files)
            warn('CAT12 segmentation failed - deformation fields not created');
        end

        disp('CAT12 segmentation completed successfully');
    else
        disp('=== STEP 1: CAT12 processing already completed, skipping ===');
    end

    % Step 1b: Collect surface-based AAL3 stats (only when surface mapping was run)
    if opts.surface
        disp('=== STEP 1b: Collecting surface AAL3 ROI stats ===');
        collect_surface_aal3_stats(newdir, t1wfiles);
    end

    % Step 2: Process DTI
    if dti_success
        disp('=== STEP 2: Processing DTI ===');
        success = process_dti_data(dti_struct, mri_dir, t1wfiles, newdir, y_files, qc_dir, template_path, opts.qc_images);
        if ~success
            warning('DTI processing failed - continuing with other modalities');
        else
            disp('DTI processing completed successfully');
        end
    else
        disp('=== STEP 2: No DTI data available, skipping ===');
    end

    % Step 3: Process NODDI
    if noddi_success
        disp('=== STEP 3: Processing NODDI ===');
        success = process_noddi_data(noddi_struct, mri_dir, t1wfiles, newdir, y_files, qc_dir, template_path, opts.qc_images);
        if ~success
            warning('NODDI processing failed - continuing with other modalities');
        else
            disp('NODDI processing completed successfully');
        end
    else
        disp('=== STEP 3: No NODDI data available, skipping ===');
    end

    % Step 4: Process ASL
    if asl_success
        disp('=== STEP 4: Processing ASL ===');
        success = process_asl_data(asl_struct, mri_dir, t1wfiles, newdir, y_files, qc_dir, template_path, opts.qc_images);
        if ~success
            warning('ASL processing failed - continuing with atlas processing');
        else
            disp('ASL processing completed successfully');
        end
    else
        disp('=== STEP 4: No ASL data available, skipping ===');
    end

    % Step 5: Process atlases to native space
    % THIS STEP CREATES NATIVE-SPACE LABELMAPS with 'w' prefix
    disp('=== STEP 5: Processing atlases to native space ===');
    process_atlases_to_native(label_map_path, label_maps_files, iy_files, newdir);
    disp('Atlas processing completed successfully');

    % Print summary
    print_processing_summary(newdir, mri_dir, dti_success, noddi_success, asl_success, qc_dir, label_maps_files, opts);

    disp('========================================');
    disp('=== ALL PROCESSING COMPLETE ===');
    disp('========================================');

end

%% CAT12 Segmentation Function
function run_cat12_segmentation(newdir, t1wfiles, spm_dir, run_surface)
    disp('  Building CAT12 batch...');

    clear matlabbatch;

    matlabbatch{1}.spm.tools.cat.estwrite.data = {fullfile(newdir, t1wfiles(end).name)};
    matlabbatch{1}.spm.tools.cat.estwrite.nproc = 8;
    matlabbatch{1}.spm.tools.cat.estwrite.useprior = '';
    matlabbatch{1}.spm.tools.cat.estwrite.opts.tpm = {fullfile(spm_dir, 'tpm', 'TPM.nii')};
    matlabbatch{1}.spm.tools.cat.estwrite.opts.affreg = 'mni';
    matlabbatch{1}.spm.tools.cat.estwrite.opts.biasacc = 0.5;
    matlabbatch{1}.spm.tools.cat.estwrite.extopts.restypes.optimal = [1 0.3];
    matlabbatch{1}.spm.tools.cat.estwrite.extopts.setCOM = 1;
    matlabbatch{1}.spm.tools.cat.estwrite.extopts.APP = 1070;
    matlabbatch{1}.spm.tools.cat.estwrite.extopts.affmod = 0;
    matlabbatch{1}.spm.tools.cat.estwrite.extopts.LASstr = 0.5;
    matlabbatch{1}.spm.tools.cat.estwrite.extopts.LASmyostr = 0;
    matlabbatch{1}.spm.tools.cat.estwrite.extopts.gcutstr = 2;
    matlabbatch{1}.spm.tools.cat.estwrite.extopts.WMHC = 2;
    matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.shooting.shootingtpm = ...
        {fullfile(spm_dir, 'toolbox', 'cat12', 'templates_MNI152NLin2009cAsym', 'Template_0_GS.nii')};
    matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.shooting.regstr = 0.5;
    matlabbatch{1}.spm.tools.cat.estwrite.extopts.vox = 1;
    matlabbatch{1}.spm.tools.cat.estwrite.extopts.bb = 12;
    matlabbatch{1}.spm.tools.cat.estwrite.extopts.SRP = 22;
    matlabbatch{1}.spm.tools.cat.estwrite.extopts.ignoreErrors = 1;
    matlabbatch{1}.spm.tools.cat.estwrite.output.BIDS.BIDSno = 1;

    % Surface mapping — controlled by caller
    if run_surface
        matlabbatch{1}.spm.tools.cat.estwrite.output.surface = 1;
        matlabbatch{1}.spm.tools.cat.estwrite.output.surf_measures = 1;
        disp('  Surface mapping: ENABLED');
    else
        matlabbatch{1}.spm.tools.cat.estwrite.output.surface = 0;
        matlabbatch{1}.spm.tools.cat.estwrite.output.surf_measures = 0;
        disp('  Surface mapping: DISABLED');
    end

    % Atlas settings
    matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.neuromorphometrics = 1;
    matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.lpba40 = 1;
    matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.cobra = 0;
    matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.hammers = 0;
    matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.thalamus = 0;
    matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.thalamic_nuclei = 0;
    matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.suit = 0;
    matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.ibsr = 0;
    matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.aal3 = 1;

    % Tissue segmentation outputs
    matlabbatch{1}.spm.tools.cat.estwrite.output.GM.native = 1;
    matlabbatch{1}.spm.tools.cat.estwrite.output.GM.warped = 1;
    matlabbatch{1}.spm.tools.cat.estwrite.output.GM.mod = 1;
    matlabbatch{1}.spm.tools.cat.estwrite.output.GM.dartel = 0;

    matlabbatch{1}.spm.tools.cat.estwrite.output.WM.native = 1;
    matlabbatch{1}.spm.tools.cat.estwrite.output.WM.warped = 1;
    matlabbatch{1}.spm.tools.cat.estwrite.output.WM.mod = 1;
    matlabbatch{1}.spm.tools.cat.estwrite.output.WM.dartel = 0;

    matlabbatch{1}.spm.tools.cat.estwrite.output.CSF.native = 1;
    matlabbatch{1}.spm.tools.cat.estwrite.output.CSF.warped = 1;
    matlabbatch{1}.spm.tools.cat.estwrite.output.CSF.mod = 1;
    matlabbatch{1}.spm.tools.cat.estwrite.output.CSF.dartel = 0;

    % Other outputs (disabled)
    matlabbatch{1}.spm.tools.cat.estwrite.output.ct.native = 0;
    matlabbatch{1}.spm.tools.cat.estwrite.output.ct.warped = 0;
    matlabbatch{1}.spm.tools.cat.estwrite.output.ct.dartel = 0;
    matlabbatch{1}.spm.tools.cat.estwrite.output.pp.native = 0;
    matlabbatch{1}.spm.tools.cat.estwrite.output.pp.warped = 0;
    matlabbatch{1}.spm.tools.cat.estwrite.output.pp.dartel = 0;
    matlabbatch{1}.spm.tools.cat.estwrite.output.WMH.native = 0;
    matlabbatch{1}.spm.tools.cat.estwrite.output.WMH.warped = 0;
    matlabbatch{1}.spm.tools.cat.estwrite.output.WMH.mod = 0;
    matlabbatch{1}.spm.tools.cat.estwrite.output.WMH.dartel = 0;
    matlabbatch{1}.spm.tools.cat.estwrite.output.SL.native = 0;
    matlabbatch{1}.spm.tools.cat.estwrite.output.SL.warped = 0;
    matlabbatch{1}.spm.tools.cat.estwrite.output.SL.mod = 0;
    matlabbatch{1}.spm.tools.cat.estwrite.output.SL.dartel = 0;
    matlabbatch{1}.spm.tools.cat.estwrite.output.TPMC.native = 0;
    matlabbatch{1}.spm.tools.cat.estwrite.output.TPMC.warped = 0;
    matlabbatch{1}.spm.tools.cat.estwrite.output.TPMC.mod = 0;
    matlabbatch{1}.spm.tools.cat.estwrite.output.TPMC.dartel = 0;
    matlabbatch{1}.spm.tools.cat.estwrite.output.atlas.native = 0;
    matlabbatch{1}.spm.tools.cat.estwrite.output.label.native = 0;
    matlabbatch{1}.spm.tools.cat.estwrite.output.label.warped = 0;
    matlabbatch{1}.spm.tools.cat.estwrite.output.label.dartel = 0;
    matlabbatch{1}.spm.tools.cat.estwrite.output.labelnative = 1;
    matlabbatch{1}.spm.tools.cat.estwrite.output.bias.warped = 1;
    matlabbatch{1}.spm.tools.cat.estwrite.output.las.native = 0;
    matlabbatch{1}.spm.tools.cat.estwrite.output.las.warped = 0;
    matlabbatch{1}.spm.tools.cat.estwrite.output.las.dartel = 0;
    matlabbatch{1}.spm.tools.cat.estwrite.output.jacobianwarped = 0;
    matlabbatch{1}.spm.tools.cat.estwrite.output.warps = [1 1];
    matlabbatch{1}.spm.tools.cat.estwrite.output.rmat = 0;

    disp('  Running CAT12 segmentation (this may take 15-30 minutes)...');
    spm_jobman('run', matlabbatch);
end

%% Surface AAL3 Stats Collection
function collect_surface_aal3_stats(newdir, t1wfiles)
    % Parse the CAT12 surface ROI XML (label/catROIs_*.xml) and write a CSV
    % containing per-ROI thickness, gyrification, and sulcal depth for AAL3.
    % Output: <newdir>/surface_aal3_stats.csv

    label_dir = fullfile(newdir, 'label');
    [~, t1_name] = fileparts(t1wfiles(end).name);

    xml_files = dir(fullfile(label_dir, 'catROIs_*.xml'));
    if isempty(xml_files)
        warning('collect_surface_aal3_stats: no catROIs_*.xml found in %s', label_dir);
        return;
    end

    xml_file = fullfile(xml_files(1).folder, xml_files(1).name);
    disp(sprintf('  Reading surface ROI stats from: %s', xml_file));

    try
        S = cat_io_xml(xml_file);
    catch ME
        warning('collect_surface_aal3_stats: failed to parse XML: %s', ME.message);
        return;
    end

    % Navigate to the AAL3 ROI block.  CAT12 organises results as either
    %   S.ROI.aal3.lh.<metric>  or  S.ROI.lh.aal3.<metric>
    % depending on version.  We try both layouts.
    aal3_lh = [];
    aal3_rh = [];

    if isfield(S, 'ROI')
        R = S.ROI;
        if isfield(R, 'aal3') && isfield(R.aal3, 'lh')
            % Layout A: S.ROI.aal3.lh / S.ROI.aal3.rh
            aal3_lh = R.aal3.lh;
            aal3_rh = R.aal3.rh;
        elseif isfield(R, 'lh') && isfield(R.lh, 'aal3')
            % Layout B: S.ROI.lh.aal3 / S.ROI.rh.aal3
            aal3_lh = R.lh.aal3;
            aal3_rh = R.rh.aal3;
        end
    end

    if isempty(aal3_lh)
        warning('collect_surface_aal3_stats: AAL3 block not found in XML. Available fields: %s', ...
            strjoin(fieldnames(S.ROI), ', '));
        return;
    end

    % Build one flat row: subject_id + lh/rh columns for each metric × ROI
    row = struct();
    row.subject_id = t1_name;

    for hemi_label = {'lh', 'rh'}
        hemi = hemi_label{1};
        if strcmp(hemi, 'lh')
            block = aal3_lh;
        else
            block = aal3_rh;
        end

        metrics = fieldnames(block);
        for m = 1:numel(metrics)
            metric = metrics{m};
            values = block.(metric);
            if isnumeric(values)
                for roi = 1:numel(values)
                    col = sprintf('aal3_%s_roi%d_%s', hemi, roi, metric);
                    row.(col) = values(roi);
                end
            end
        end
    end

    csv_path = fullfile(newdir, 'surface_aal3_stats.csv');
    try
        T = struct2table(row);
        writetable(T, csv_path);
        disp(sprintf('  Surface AAL3 stats written to: %s', csv_path));
        disp(sprintf('  Columns written: %d', width(T)));
    catch ME
        warning('collect_surface_aal3_stats: failed to write CSV: %s', ME.message);
    end
end

%% DTI Processing Function
function success = process_dti_data(dti_struct, mri_dir, t1wfiles, newdir, y_files, qc_dir, template_path, qc_images)
    success = false;

    try
        % Get T1 native space image
        t1_native = fullfile(mri_dir, ['m' t1wfiles(end).name]);
        if ~exist(t1_native, 'file')
            t1_native = fullfile(newdir, t1wfiles(end).name);
        end

        if ~exist(t1_native, 'file')
            error('T1 native space image not found');
        end

        % Create output directory
        dti_norm_dir = fullfile(mri_dir, 'DTI');
        if ~exist(dti_norm_dir, 'dir')
            mkdir(dti_norm_dir);
        end

        % Prepare files for coregistration
        dti_files = struct2cell(dti_struct);
        ref_image = dti_struct.S0;
        source_images = dti_files(~strcmp(dti_files, ref_image));

        % Verify all source files exist
        if ~exist(ref_image, 'file')
            error('DTI reference file (S0) not found: %s', ref_image);
        end
        for i = 1:length(source_images)
            if ~exist(source_images{i}, 'file')
                error('DTI source file not found: %s', source_images{i});
            end
        end

        % SUBSTEP 2A: Coregister DTI to native T1 space
        disp('  Substep 2a: Coregistering DTI to native T1 space...');
        clear matlabbatch;

        matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {t1_native};
        matlabbatch{1}.spm.spatial.coreg.estwrite.source = {ref_image};
        matlabbatch{1}.spm.spatial.coreg.estwrite.other = source_images;
        matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
        matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
        matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
        matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 4;
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';

        spm_jobman('run', matlabbatch);

        % Validate coregistered files exist
        [pth, nam, ext] = fileparts(ref_image);
        coreg_ref = fullfile(pth, ['r' nam ext]);
        if ~exist(coreg_ref, 'file')
            error('Coregistered DTI reference file not created: %s', coreg_ref);
        end
        disp('  DTI coregistration successful');

        % Move native-space (r*) files to DTI_native/
        disp('  Moving native-space DTI (r*) files...');
        dti_native_dir = fullfile(mri_dir, 'DTI_native');
        if ~exist(dti_native_dir, 'dir')
            mkdir(dti_native_dir);
        end

        coregistered_dti_files = cell(size(dti_files));
        for i = 1:length(dti_files)
            [pth, nam, ext] = fileparts(dti_files{i});
            original_r_file = fullfile(pth, ['r' nam ext]);
            [~, fname, fext] = fileparts(original_r_file);
            new_r_file_path = fullfile(dti_native_dir, [fname fext]);

            if ~exist(original_r_file, 'file')
                error('Coregistered file not found: %s', original_r_file);
            end

            [move_success, msg] = movefile(original_r_file, new_r_file_path, 'f');
            if ~move_success
                error('Failed to move file %s to %s: %s', original_r_file, new_r_file_path, msg);
            end

            coregistered_dti_files{i} = new_r_file_path;

            if qc_images
                [~, fname_raw] = fileparts(dti_files{i});
                qc_name = sprintf('QC_DTI_%s_vs_T1_Native.png', fname_raw);
                generate_reg_qc(t1_native, new_r_file_path, fullfile(qc_dir, qc_name));
            end
        end
        disp(sprintf('  Successfully moved %d files to %s', length(dti_files), dti_native_dir));

        % SUBSTEP 2B: Normalize coregistered DTI files
        disp('  Substep 2b: Normalizing coregistered DTI files...');

        if isempty(y_files)
            error('Forward deformation field not found');
        end

        clear matlabbatch;
        matlabbatch{1}.spm.util.defs.comp{1}.def = {fullfile(y_files(1).folder, y_files(1).name)};
        matlabbatch{1}.spm.util.defs.out{1}.pull.fnames = coregistered_dti_files;
        matlabbatch{1}.spm.util.defs.out{1}.pull.savedir.saveusr = {dti_norm_dir};
        matlabbatch{1}.spm.util.defs.out{1}.pull.interp = 4;
        matlabbatch{1}.spm.util.defs.out{1}.pull.mask = 1;
        matlabbatch{1}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
        matlabbatch{1}.spm.util.defs.out{1}.pull.prefix = 'w';

        spm_jobman('run', matlabbatch);

        % Validate normalized files
        [~, nam, ext] = fileparts(ref_image);
        norm_ref = fullfile(dti_norm_dir, ['wr' nam ext]);
        if ~exist(norm_ref, 'file')
            error('Normalized DTI reference file not created: %s', norm_ref);
        end

        disp('  DTI normalization successful');

        if qc_images
            disp('    Generating QC images for normalized DTI files...');
            for i = 1:length(coregistered_dti_files)
                [~, fname, fext] = fileparts(coregistered_dti_files{i});
                norm_file = fullfile(dti_norm_dir, ['w' fname fext]);
                [~, fname_raw] = fileparts(dti_files{i});
                qc_name = sprintf('QC_DTI_%s_vs_Template.png', ['w' fname_raw]);
                generate_reg_qc(template_path, norm_file, fullfile(qc_dir, qc_name));
            end
        end

        success = true;

    catch ME
        warning('DTI processing error: %s', ME.message);
        disp('  Stack trace:');
        for i = 1:length(ME.stack)
            disp(sprintf('    %s (line %d)', ME.stack(i).name, ME.stack(i).line));
        end
        success = false;
    end
end

%% NODDI Processing Function
function success = process_noddi_data(noddi_struct, mri_dir, t1wfiles, newdir, y_files, qc_dir, template_path, qc_images)
    success = false;

    try
        % Get T1 native space image
        t1_native = fullfile(mri_dir, ['m' t1wfiles(end).name]);
        if ~exist(t1_native, 'file')
            t1_native = fullfile(newdir, t1wfiles(end).name);
        end

        if ~exist(t1_native, 'file')
            error('T1 native space image not found');
        end

        % Create output directory
        noddi_norm_dir = fullfile(mri_dir, 'NODDI');
        if ~exist(noddi_norm_dir, 'dir')
            mkdir(noddi_norm_dir);
        end

        % Prepare files for coregistration
        noddi_files = struct2cell(noddi_struct);
        ref_image = noddi_files{1};
        source_images = noddi_files(2:end);

        % Verify all files exist
        if ~exist(ref_image, 'file')
            error('NODDI reference file not found: %s', ref_image);
        end
        for i = 1:length(source_images)
            if ~exist(source_images{i}, 'file')
                error('NODDI source file not found: %s', source_images{i});
            end
        end

        % SUBSTEP 3A: Coregister NODDI to native T1 space
        disp('  Substep 3a: Coregistering NODDI to native T1 space...');
        clear matlabbatch;

        matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {t1_native};
        matlabbatch{1}.spm.spatial.coreg.estwrite.source = {ref_image};
        matlabbatch{1}.spm.spatial.coreg.estwrite.other = source_images;
        matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
        matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
        matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
        matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 4;
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
        matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';

        spm_jobman('run', matlabbatch);

        % Validate coregistered files exist
        [pth, nam, ext] = fileparts(ref_image);
        coreg_ref = fullfile(pth, ['r' nam ext]);
        if ~exist(coreg_ref, 'file')
            error('Coregistered NODDI reference file not created: %s', coreg_ref);
        end
        disp('  NODDI coregistration successful');

        % Move native-space (r*) files to NODDI_native/
        disp('  Moving native-space NODDI (r*) files...');
        noddi_native_dir = fullfile(mri_dir, 'NODDI_native');
        if ~exist(noddi_native_dir, 'dir')
            mkdir(noddi_native_dir);
        end

        coregistered_noddi_files = cell(size(noddi_files));
        for i = 1:length(noddi_files)
            [pth, nam, ext] = fileparts(noddi_files{i});
            original_r_file = fullfile(pth, ['r' nam ext]);
            [~, fname, fext] = fileparts(original_r_file);
            new_r_file_path = fullfile(noddi_native_dir, [fname fext]);

            if ~exist(original_r_file, 'file')
                error('Coregistered file not found: %s', original_r_file);
            end

            [move_success, msg] = movefile(original_r_file, new_r_file_path, 'f');
            if ~move_success
                error('Failed to move file %s to %s: %s', original_r_file, new_r_file_path, msg);
            end

            coregistered_noddi_files{i} = new_r_file_path;

            if qc_images
                [~, fname_raw] = fileparts(noddi_files{i});
                qc_name = sprintf('QC_NODDI_%s_vs_T1_Native.png', fname_raw);
                generate_reg_qc(t1_native, new_r_file_path, fullfile(qc_dir, qc_name));
            end
        end
        disp(sprintf('  Successfully moved %d files to %s', length(noddi_files), noddi_native_dir));

        % SUBSTEP 3B: Normalize coregistered NODDI files
        disp('  Substep 3b: Normalizing coregistered NODDI files...');

        if isempty(y_files)
            error('Forward deformation field not found');
        end

        clear matlabbatch;
        matlabbatch{1}.spm.util.defs.comp{1}.def = {fullfile(y_files(1).folder, y_files(1).name)};
        matlabbatch{1}.spm.util.defs.out{1}.pull.fnames = coregistered_noddi_files;
        matlabbatch{1}.spm.util.defs.out{1}.pull.savedir.saveusr = {noddi_norm_dir};
        matlabbatch{1}.spm.util.defs.out{1}.pull.interp = 4;
        matlabbatch{1}.spm.util.defs.out{1}.pull.mask = 1;
        matlabbatch{1}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
        matlabbatch{1}.spm.util.defs.out{1}.pull.prefix = 'w';

        spm_jobman('run', matlabbatch);

        % Validate normalized files
        [~, nam, ext] = fileparts(ref_image);
        norm_ref = fullfile(noddi_norm_dir, ['wr' nam ext]);
        if ~exist(norm_ref, 'file')
            error('Normalized NODDI reference file not created: %s', norm_ref);
        end

        disp('  NODDI normalization successful');

        if qc_images
            disp('    Generating QC images for normalized NODDI files...');
            for i = 1:length(coregistered_noddi_files)
                [~, fname, fext] = fileparts(coregistered_noddi_files{i});
                norm_file = fullfile(noddi_norm_dir, ['w' fname fext]);
                [~, fname_raw] = fileparts(noddi_files{i});
                qc_name = sprintf('QC_NODDI_%s_vs_Template.png', ['w' fname_raw]);
                generate_reg_qc(template_path, norm_file, fullfile(qc_dir, qc_name));
            end
        end

        success = true;

    catch ME
        warning('NODDI processing error: %s', ME.message);
        disp('  Stack trace:');
        for i = 1:length(ME.stack)
            disp(sprintf('    %s (line %d)', ME.stack(i).name, ME.stack(i).line));
        end
        success = false;
    end
end

%% ASL Processing Function
function success = process_asl_data(asl_struct, mri_dir, t1wfiles, newdir, y_files, qc_dir, template_path, qc_images)
    success = false;

    try
        % Get T1 native space image
        t1_native = fullfile(mri_dir, ['m' t1wfiles(end).name]);
        if ~exist(t1_native, 'file')
            t1_native = fullfile(newdir, t1wfiles(end).name);
        end

        if ~exist(t1_native, 'file')
            error('T1 native space image not found');
        end

        % Create output directory
        asl_norm_dir = fullfile(mri_dir, 'ASL');
        if ~exist(asl_norm_dir, 'dir')
            mkdir(asl_norm_dir);
        end

        % Handle multiple ASL datasets
        all_coregistered_asl_files = {};

        for asl_idx = 1:length(asl_struct)
            disp(sprintf('  Processing ASL dataset %d of %d...', asl_idx, length(asl_struct)));

            asl_files_cell = struct2cell(asl_struct(asl_idx));
            asl_files_current = asl_files_cell(~cellfun('isempty', asl_files_cell));

            if isempty(asl_files_current)
                continue;
            end

            ref_image = asl_struct(asl_idx).M0_masked;
            source_images = asl_files_current(~strcmp(asl_files_current, ref_image));

            % Verify all files exist
            if ~exist(ref_image, 'file')
                error('ASL reference file (M0) not found: %s', ref_image);
            end
            for i = 1:length(source_images)
                if ~exist(source_images{i}, 'file')
                    error('ASL source file not found: %s', source_images{i});
                end
            end

            % SUBSTEP 4A: Coregister ASL to native T1 space
            disp(sprintf('  Substep 4a-%d: Coregistering ASL dataset %d...', asl_idx, asl_idx));
            clear matlabbatch;

            matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {t1_native};
            matlabbatch{1}.spm.spatial.coreg.estwrite.source = {ref_image};
            matlabbatch{1}.spm.spatial.coreg.estwrite.other = source_images;
            matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
            matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
            matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
            matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
            matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 4;
            matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
            matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
            matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';

            spm_jobman('run', matlabbatch);

            % Validate coregistered files exist
            [pth, nam, ext] = fileparts(ref_image);
            coreg_ref = fullfile(pth, ['r' nam ext]);
            if ~exist(coreg_ref, 'file')
                error('Coregistered ASL reference file not created: %s', coreg_ref);
            end
            disp(sprintf('  ASL dataset %d coregistration successful', asl_idx));

            % Move native-space (r*) files to ASL_native_run<N>/
            disp(sprintf('  Moving native-space ASL (r*) files for dataset %d...', asl_idx));
            asl_native_dir = fullfile(mri_dir, sprintf('ASL_native_run%d', asl_idx));
            if ~exist(asl_native_dir, 'dir')
                mkdir(asl_native_dir);
            end

            coregistered_asl_files = cell(size(asl_files_current));
            for i = 1:length(asl_files_current)
                [pth, nam, ext] = fileparts(asl_files_current{i});
                original_r_file = fullfile(pth, ['r' nam ext]);
                [~, fname, fext] = fileparts(original_r_file);
                new_r_file_path = fullfile(asl_native_dir, [fname fext]);

                if ~exist(original_r_file, 'file')
                    error('Coregistered file not found: %s', original_r_file);
                end

                [move_success, msg] = movefile(original_r_file, new_r_file_path, 'f');
                if ~move_success
                    error('Failed to move file %s to %s: %s', original_r_file, new_r_file_path, msg);
                end

                coregistered_asl_files{i} = new_r_file_path;

                if qc_images
                    [~, fname_raw] = fileparts(asl_files_current{i});
                    qc_name = sprintf('QC_ASL_Run%d_%s_vs_T1_Native.png', asl_idx, fname_raw);
                    generate_reg_qc(t1_native, new_r_file_path, fullfile(qc_dir, qc_name));
                end
            end
            disp(sprintf('  Successfully moved %d files to %s', length(asl_files_current), asl_native_dir));

            % SUBSTEP 4B: Normalize coregistered ASL files
            disp(sprintf('  Substep 4b-%d: Normalizing coregistered ASL dataset %d...', asl_idx, asl_idx));

            if isempty(y_files)
                error('Forward deformation field not found');
            end

            clear matlabbatch;
            matlabbatch{1}.spm.util.defs.comp{1}.def = {fullfile(y_files(1).folder, y_files(1).name)};
            matlabbatch{1}.spm.util.defs.out{1}.pull.fnames = coregistered_asl_files;
            matlabbatch{1}.spm.util.defs.out{1}.pull.savedir.saveusr = {asl_norm_dir};
            matlabbatch{1}.spm.util.defs.out{1}.pull.interp = 4;
            matlabbatch{1}.spm.util.defs.out{1}.pull.mask = 1;
            matlabbatch{1}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
            matlabbatch{1}.spm.util.defs.out{1}.pull.prefix = 'w';

            spm_jobman('run', matlabbatch);

            % Validate normalized files
            [~, nam, ext] = fileparts(ref_image);
            norm_ref = fullfile(asl_norm_dir, ['wr' nam ext]);
            if ~exist(norm_ref, 'file')
                error('Normalized ASL reference file not created: %s', norm_ref);
            end

            disp(sprintf('  ASL dataset %d normalization successful', asl_idx));

            if qc_images
                disp(sprintf('    Generating QC images for normalized ASL dataset %d...', asl_idx));
                for i = 1:length(coregistered_asl_files)
                    [~, fname, fext] = fileparts(coregistered_asl_files{i});
                    norm_file = fullfile(asl_norm_dir, ['w' fname fext]);
                    [~, fname_raw] = fileparts(asl_files_current{i});
                    qc_name = sprintf('QC_ASL_Run%d_%s_vs_Template.png', asl_idx, ['w' fname_raw]);
                    generate_reg_qc(template_path, norm_file, fullfile(qc_dir, qc_name));
                end
            end
        end

        disp('  All ASL datasets processed successfully');
        success = true;

    catch ME
        warning('ASL processing error: %s', ME.message);
        disp('  Stack trace:');
        for i = 1:length(ME.stack)
            disp(sprintf('    %s (line %d)', ME.stack(i).name, ME.stack(i).line));
        end
        success = false;
    end
end

%% Atlas Processing Function - Creates Native-Space Labelmaps
function process_atlases_to_native(label_map_path, label_maps_files, iy_files, newdir)
    % This function transforms template-space atlases to native T1 space
    % Output: Native-space labelmap files with 'w' prefix in newdir

    disp('  Processing label maps to native space...');

    if isempty(iy_files)
        error('Inverse deformation field not found');
    end

    for atlas_idx = 1:length(label_maps_files)
        current_atlas = label_maps_files{atlas_idx};
        atlas_file = fullfile(label_map_path, current_atlas);

        % Handle .nii.gz files - SPM may need them uncompressed
        if endsWith(current_atlas, '.nii.gz')
            % Check if uncompressed version exists
            atlas_file_nii = atlas_file(1:end-3); % Remove .gz
            if exist(atlas_file_nii, 'file')
                disp(sprintf('  Using uncompressed version: %s', atlas_file_nii));
                atlas_file = atlas_file_nii;
            elseif exist(atlas_file, 'file')
                % SPM12 can usually handle .nii.gz, but with ,1 volume specification
                atlas_file = [atlas_file ',1'];
                disp(sprintf('  Using gzipped file: %s', current_atlas));
            else
                warning('Atlas file not found: %s', atlas_file);
                continue;
            end
        elseif ~exist(atlas_file, 'file')
            warning('Atlas file not found: %s', atlas_file);
            continue;
        end

        disp(sprintf('  Processing atlas %d/%d: %s', atlas_idx, length(label_maps_files), current_atlas));

        clear matlabbatch;
        matlabbatch{1}.spm.util.defs.comp{1}.def = {fullfile(iy_files(1).folder, iy_files(1).name)};
        matlabbatch{1}.spm.util.defs.out{1}.pull.fnames = {atlas_file};
        matlabbatch{1}.spm.util.defs.out{1}.pull.savedir.saveusr = {newdir};
        matlabbatch{1}.spm.util.defs.out{1}.pull.interp = 0;  % Nearest neighbor for label maps
        matlabbatch{1}.spm.util.defs.out{1}.pull.mask = 1;
        matlabbatch{1}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
        matlabbatch{1}.spm.util.defs.out{1}.pull.prefix = 'w';

        try
            spm_jobman('run', matlabbatch);

            % Verify output was created
            [~, name, ~] = fileparts(strrep(current_atlas, '.gz', ''));
            output_file = fullfile(newdir, ['w' name '.nii']);
            if exist(output_file, 'file')
                disp(sprintf('  Successfully created: %s', output_file));
            else
                warning('Output file not created: %s', output_file);
            end
        catch ME
            warning('Failed to process atlas %s: %s', current_atlas, ME.message);
        end
    end

    disp('  Atlas processing complete');
end

%% Summary Function
function print_processing_summary(newdir, mri_dir, dti_success, noddi_success, asl_success, qc_dir, label_maps_files, opts)
    disp('');
    disp('========================================');
    disp('=== PROCESSING SUMMARY ===');
    disp('========================================');
    disp(sprintf('Output directory: %s', newdir));
    disp(sprintf('MRI outputs: %s', mri_dir));
    disp(sprintf('Options: hypothalamus=%d  surface=%d  qc_images=%d', ...
        opts.hypothalamus, opts.surface, opts.qc_images));
    if opts.qc_images
        disp(sprintf('QC outputs: %s', qc_dir));
    else
        disp('QC outputs: disabled');
    end
    disp('');

    if dti_success
        disp('DTI Processing:');
        disp(sprintf('  Native coregistered (r* prefix): %s', fullfile(mri_dir, 'DTI_native')));
        disp(sprintf('  Normalized (wr* prefix): %s', fullfile(mri_dir, 'DTI')));
    else
        disp('DTI Processing: FAILED or NOT AVAILABLE');
    end
    disp('');

    if noddi_success
        disp('NODDI Processing:');
        disp(sprintf('  Native coregistered (r* prefix): %s', fullfile(mri_dir, 'NODDI_native')));
        disp(sprintf('  Normalized (wr* prefix): %s', fullfile(mri_dir, 'NODDI')));
    else
        disp('NODDI Processing: FAILED or NOT AVAILABLE');
    end
    disp('');

    if asl_success
        disp('ASL Processing:');
        disp(sprintf('  Native coregistered (r* prefix): %s', fullfile(mri_dir, 'ASL_native_run*')));
        disp(sprintf('  Normalized (wr* prefix): %s', fullfile(mri_dir, 'ASL')));
    else
        disp('ASL Processing: FAILED or NOT AVAILABLE');
    end
    disp('');

    disp('Native Space Atlases (LABELMAPS):');
    disp(sprintf('  Location: %s', newdir));
    disp('  Prefix: w*');
    disp('  Files created:');
    for i = 1:length(label_maps_files)
        [~, name, ~] = fileparts(label_maps_files{i});
        disp(sprintf('    - w%s.nii', name));
    end
    disp('');

    if opts.surface
        disp('Surface Processing:');
        disp(sprintf('  Surface files: %s', fullfile(newdir, 'surf')));
        disp(sprintf('  AAL3 stats CSV: %s', fullfile(newdir, 'surface_aal3_stats.csv')));
        disp('');
    end

    disp('Tissue Segmentations:');
    disp(sprintf('  Location: %s', mri_dir));
    disp('  Warped GM: mwp1*');
    disp('  Warped WM: mwp2*');
    disp('');

    disp('Usage Notes:');
    disp('  For native space analysis: Use r* files (now in *_native dirs) with w* atlas files');
    disp('  For template space analysis: Use wr* files with original template atlases');
    if opts.qc_images
        disp('  Check QC_registration/ folder for registration quality control images');
    end
    disp('========================================');
end

%% QC Generation Function
function generate_reg_qc(ref_img, source_img, output_filename)
    disp(sprintf('    Generating QC image: %s', output_filename));
    try
        % Make sure Graphics window is available
        spm_figure('GetWin', 'Graphics');
        spm_figure('Clear', 'Graphics');

        % Check registration - this displays the images in orthogonal views
        spm_check_registration(char({ref_img, source_img}));

        % Print to file
        print(gcf, output_filename, '-dpng', '-r150');
    catch ME
        warning('Failed to generate QC image: %s', ME.message);
    end
end
