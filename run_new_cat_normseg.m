function run_new_cat_normseg(base_dir)
    % NATIVE-SPACE LABELMAP EXPORT: This script creates native-space labelmaps
    % in Step 5 using the process_atlases_to_native() function.
    % Output files will have 'w' prefix (e.g., wneuromorphometrics.nii)

    % Add SPM12 and CAT12 to the MATLAB path
    addpath('/isilon/datalake/riipl/original/DEMONco/Hellcat-12.9/libs/spm12/spm12/');
    addpath('/isilon/datalake/riipl/original/DEMONco/Hellcat-12.9/');

    % Updated atlas list including hypothalamic atlas
    label_map_path='/isilon/datalake/riipl/original/DEMONco/Hellcat-12.9/libs/spm12/spm12/toolbox/cat12/templates_MNI152NLin2009cAsym/';
    label_maps_files = {'neuromorphometrics.nii','aal3.nii','cobra.nii','hammers.nii','lpba40.nii','ibsr.nii','JHU.nii','brainmask_T1.nii','hypothalamusAtlas.nii'};

    % Initialize SPM
    spm('defaults', 'FMRI');
    spm_jobman('initcfg');

    % Create the new directory
    newdir = fullfile(base_dir, 'nifti', 'cat12_v2560');
    if ~exist(newdir, 'dir')
        mkdir(newdir);
    end

    % Find the T1w image
    t1wfiles = dir(fullfile([base_dir, '/nifti/'], '3*-tfl3d116ns.nii*'));
    if isempty(t1wfiles)
        error('No T1w image found');
    end

    % Copy the T1w image to the new directory
    copyfile(fullfile(t1wfiles(1).folder, t1wfiles(1).name), newdir, 'f');

    % Check if processing has already been done
    mri_dir = fullfile(newdir, 'mri');
    y_files = dir(fullfile(mri_dir, 'y_*-tfl3d116ns.nii'));
    iy_files = dir(fullfile(mri_dir, 'iy_*-tfl3d116ns.nii'));

    disp('========================================');
    disp(sprintf('Subject dir: %s', base_dir));
    disp(sprintf('MRI dir: %s', mri_dir));
    disp(sprintf('Input T1: %s', fullfile(newdir, t1wfiles(end).name)));
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
        run_cat12_segmentation(newdir, t1wfiles);
        
        % Re-check for deformation fields
        y_files = dir(fullfile(mri_dir, 'y_*-tfl3d116ns.nii'));
        iy_files = dir(fullfile(mri_dir, 'iy_*-tfl3d116ns.nii'));
        
        if isempty(y_files) || isempty(iy_files)
            error('CAT12 segmentation failed - deformation fields not created');
        end
        disp('CAT12 segmentation completed successfully');
    else
        disp('=== STEP 1: CAT12 processing already completed, skipping ===');
    end

    % Step 2: Process DTI
    if dti_success
        disp('=== STEP 2: Processing DTI ===');
        success = process_dti_data(dti_struct, mri_dir, t1wfiles, newdir, y_files);
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
        success = process_noddi_data(noddi_struct, mri_dir, t1wfiles, newdir, y_files);
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
        success = process_asl_data(asl_struct, mri_dir, t1wfiles, newdir, y_files);
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
    print_processing_summary(newdir, mri_dir, dti_success, noddi_success, asl_success);
    
    disp('========================================');
    disp('=== ALL PROCESSING COMPLETE ===');
    disp('========================================');

end

%% CAT12 Segmentation Function
function run_cat12_segmentation(newdir, t1wfiles)
    disp('  Building CAT12 batch...');
    
    clear matlabbatch;
    
    matlabbatch{1}.spm.tools.cat.estwrite.data = {fullfile(newdir, t1wfiles(end).name)};
    matlabbatch{1}.spm.tools.cat.estwrite.nproc = 8;
    matlabbatch{1}.spm.tools.cat.estwrite.useprior = '';
    matlabbatch{1}.spm.tools.cat.estwrite.opts.tpm = {'/isilon/datalake/riipl/original/DEMONco/Hellcat-12.9/libs/spm12/spm12/tpm/TPM.nii'};
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
    matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.shooting.shootingtpm = {'/isilon/datalake/riipl/original/DEMONco/Hellcat-12.9/libs/spm12/spm12/toolbox/cat12/templates_MNI152NLin2009cAsym/Template_0_GS.nii'};
    matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.shooting.regstr = 0.5;
    matlabbatch{1}.spm.tools.cat.estwrite.extopts.vox = 1;
    matlabbatch{1}.spm.tools.cat.estwrite.extopts.bb = 12;
    matlabbatch{1}.spm.tools.cat.estwrite.extopts.SRP = 22;
    matlabbatch{1}.spm.tools.cat.estwrite.extopts.ignoreErrors = 1;
    matlabbatch{1}.spm.tools.cat.estwrite.output.BIDS.BIDSno = 1;
    
    % Disable surface mapping
    matlabbatch{1}.spm.tools.cat.estwrite.output.surface = 0;
    matlabbatch{1}.spm.tools.cat.estwrite.output.surf_measures = 0;
    
    % Atlas settings
    matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.neuromorphometrics = 1;
    matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.lpba40 = 1;
    matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.cobra = 0;
    matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.hammers = 0;
    matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.thalamus = 0;
    matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.thalamic_nuclei = 0;   
    matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.suit = 0;
    matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.ibsr = 0;
    
    % Tissue segmentation outputs
    matlabbatch{1}.spm.tools.cat.estwrite.output.GM.native = 1;
    matlabbatch{1}.spm.tools.cat.estwrite.output.GM.warped = 1;
    matlabbatch{1}.spm.tools.cat.estwrite.output.GM.mod = 0;
    matlabbatch{1}.spm.tools.cat.estwrite.output.GM.dartel = 0;
    matlabbatch{1}.spm.tools.cat.estwrite.output.WM.native = 1;
    matlabbatch{1}.spm.tools.cat.estwrite.output.WM.warped = 1;
    matlabbatch{1}.spm.tools.cat.estwrite.output.WM.mod = 0;
    matlabbatch{1}.spm.tools.cat.estwrite.output.WM.dartel = 0;
    matlabbatch{1}.spm.tools.cat.estwrite.output.CSF.native = 1;
    matlabbatch{1}.spm.tools.cat.estwrite.output.CSF.warped = 0;
    matlabbatch{1}.spm.tools.cat.estwrite.output.CSF.mod = 0;
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
    matlabbatch{1}.spm.tools.cat.estwrite.output.label.native = 1;
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

%% DTI Processing Function
function success = process_dti_data(dti_struct, mri_dir, t1wfiles, newdir, y_files)
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
        
        % SUBSTEP 2B: Normalize coregistered DTI files
        disp('  Substep 2b: Normalizing coregistered DTI files...');
        
        if isempty(y_files)
            error('Forward deformation field not found');
        end
        
        % Create list of coregistered files and verify they exist
        coregistered_dti_files = cell(size(dti_files));
        for i = 1:length(dti_files)
            [pth, nam, ext] = fileparts(dti_files{i});
            coregistered_dti_files{i} = fullfile(pth, ['r' nam ext]);
            
            if ~exist(coregistered_dti_files{i}, 'file')
                error('Coregistered file not found: %s', coregistered_dti_files{i});
            end
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
function success = process_noddi_data(noddi_struct, mri_dir, t1wfiles, newdir, y_files)
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
        
        % SUBSTEP 3B: Normalize coregistered NODDI files
        disp('  Substep 3b: Normalizing coregistered NODDI files...');
        
        if isempty(y_files)
            error('Forward deformation field not found');
        end
        
        % Create list of coregistered files and verify they exist
        coregistered_noddi_files = cell(size(noddi_files));
        for i = 1:length(noddi_files)
            [pth, nam, ext] = fileparts(noddi_files{i});
            coregistered_noddi_files{i} = fullfile(pth, ['r' nam ext]);
            
            if ~exist(coregistered_noddi_files{i}, 'file')
                error('Coregistered file not found: %s', coregistered_noddi_files{i});
            end
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
function success = process_asl_data(asl_struct, mri_dir, t1wfiles, newdir, y_files)
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
            
            % SUBSTEP 4B: Normalize coregistered ASL files
            disp(sprintf('  Substep 4b-%d: Normalizing coregistered ASL dataset %d...', asl_idx, asl_idx));
            
            if isempty(y_files)
                error('Forward deformation field not found');
            end
            
            % Create list of coregistered files and verify they exist
            coregistered_asl_files = cell(size(asl_files_current));
            for i = 1:length(asl_files_current)
                [pth, nam, ext] = fileparts(asl_files_current{i});
                coregistered_asl_files{i} = fullfile(pth, ['r' nam ext]);
                
                if ~exist(coregistered_asl_files{i}, 'file')
                    error('Coregistered file not found: %s', coregistered_asl_files{i});
                end
            end
            
            % Add to master list
            all_coregistered_asl_files = [all_coregistered_asl_files; coregistered_asl_files];
            
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
function print_processing_summary(newdir, mri_dir, dti_success, noddi_success, asl_success)
    disp('');
    disp('========================================');
    disp('=== PROCESSING SUMMARY ===');
    disp('========================================');
    disp(sprintf('Output directory: %s', newdir));
    disp(sprintf('MRI outputs: %s', mri_dir));
    disp('');
    
    if dti_success
        disp('DTI Processing:');
        disp(sprintf('  Native coregistered (r* prefix): Original DTI directory'));
        disp(sprintf('  Normalized (wr* prefix): %s', fullfile(mri_dir, 'DTI')));
    else
        disp('DTI Processing: FAILED or NOT AVAILABLE');
    end
    disp('');
    
    if noddi_success
        disp('NODDI Processing:');
        disp(sprintf('  Native coregistered (r* prefix): Original NODDI directory'));
        disp(sprintf('  Normalized (wr* prefix): %s', fullfile(mri_dir, 'NODDI')));
    else
        disp('NODDI Processing: FAILED or NOT AVAILABLE');
    end
    disp('');
    
    if asl_success
        disp('ASL Processing:');
        disp(sprintf('  Native coregistered (r* prefix): Original ASL directories'));
        disp(sprintf('  Normalized (wr* prefix): %s', fullfile(mri_dir, 'ASL')));
    else
        disp('ASL Processing: FAILED or NOT AVAILABLE');
    end
    disp('');
    
    disp('Native Space Atlases (LABELMAPS):');
    disp(sprintf('  Location: %s', newdir));
    disp('  Prefix: w*');
    disp('  Files created:');
    disp('    - wneuromorphometrics.nii');
    disp('    - waal3.nii');
    disp('    - wcobra.nii');
    disp('    - whammers.nii');
    disp('    - wlpba40.nii');
    disp('    - wibsr.nii');
    disp('    - wJHU.nii');
    disp('    - wbrainmask_T1.nii');
    disp('    - whypothalamusAtlas.nii');
    disp('');
    
    disp('Tissue Segmentations:');
    disp(sprintf('  Location: %s', mri_dir));
    disp('  Warped GM: mwp1*');
    disp('  Warped WM: mwp2*');
    disp('');
    
    disp('Usage Notes:');
    disp('  For native space analysis: Use r* files with w* atlas files');
    disp('  For template space analysis: Use wr* files with original template atlases');
    disp('========================================');
end