function run_new_cat_normseg(base_dir)

    % Add SPM12 and CAT12 to the MATLAB path
    addpath('/isilon/datalake/riipl/original/DEMONco/Hellcat-12.9/libs/spm12/spm12/');
    addpath('/isilon/datalake/riipl/original/DEMONco/Hellcat-12.9/');

    % Updated atlas list including hypothalamic atlas
    label_map_path='/isilon/datalake/riipl/original/DEMONco/Hellcat-12.9/libs/spm12/spm12/toolbox/cat12/templates_MNI152NLin2009cAsym/';
    label_maps_files = {'neuromorphometrics.nii','aal3.nii','cobra.nii','hammers.nii','lpba40.nii','ibsr.nii','JHU.nii','brainmask_T1.nii','hypothalamusAtlas.nii.gz'};

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

    disp(sprintf("Subject dir %s", base_dir));
    disp(sprintf("MRI dir %s", mri_dir));
    disp(sprintf("Input T1 %s", fullfile(newdir, t1wfiles(end).name)));

    % Get DTI, NODDI, and ASL files
    [asl_struct, asl_success] = getASLFiles(base_dir);
    [dti_struct, dti_success] = getDTIFiles(base_dir);
    [noddi_struct, noddi_success] = getNODDIFiles(base_dir);

    % Initialize batch counter
    batch_idx = 1;
    
    % Step 1: CAT12 Segmentation and Normalization (if not already done)
    if isempty(y_files)
        disp('Running CAT12 segmentation and normalization...');
        
        spm('defaults', 'FMRI');
        spm_jobman('initcfg');

        matlabbatch{batch_idx}.spm.tools.cat.estwrite.data = {fullfile(newdir, t1wfiles(end).name)};
        matlabbatch{batch_idx}.spm.tools.cat.estwrite.nproc = 8;
        matlabbatch{batch_idx}.spm.tools.cat.estwrite.useprior = '';
        matlabbatch{batch_idx}.spm.tools.cat.estwrite.opts.tpm = {'/isilon/datalake/riipl/original/DEMONco/Hellcat-12.9/libs/spm12/spm12/tpm/TPM.nii'};
        matlabbatch{batch_idx}.spm.tools.cat.estwrite.opts.affreg = 'mni';
        matlabbatch{batch_idx}.spm.tools.cat.estwrite.opts.biasacc = 0.5;
        matlabbatch{batch_idx}.spm.tools.cat.estwrite.extopts.restypes.optimal = [1 0.3];
        matlabbatch{batch_idx}.spm.tools.cat.estwrite.extopts.setCOM = 1;
        matlabbatch{batch_idx}.spm.tools.cat.estwrite.extopts.APP = 1070;
        matlabbatch{batch_idx}.spm.tools.cat.estwrite.extopts.affmod = 0;
        matlabbatch{batch_idx}.spm.tools.cat.estwrite.extopts.LASstr = 0.5;
        matlabbatch{batch_idx}.spm.tools.cat.estwrite.extopts.LASmyostr = 0;
        matlabbatch{batch_idx}.spm.tools.cat.estwrite.extopts.gcutstr = 2;
        matlabbatch{batch_idx}.spm.tools.cat.estwrite.extopts.WMHC = 2;
        matlabbatch{batch_idx}.spm.tools.cat.estwrite.extopts.registration.shooting.shootingtpm = {'/isilon/datalake/riipl/original/DEMONco/Hellcat-12.9/libs/spm12/spm12/toolbox/cat12/templates_MNI152NLin2009cAsym/Template_0_GS.nii'};
        matlabbatch{batch_idx}.spm.tools.cat.estwrite.extopts.registration.shooting.regstr = 0.5;
        matlabbatch{batch_idx}.spm.tools.cat.estwrite.extopts.vox = 1.5;
        matlabbatch{batch_idx}.spm.tools.cat.estwrite.extopts.bb = 12;
        matlabbatch{batch_idx}.spm.tools.cat.estwrite.extopts.SRP = 22;
        matlabbatch{batch_idx}.spm.tools.cat.estwrite.extopts.ignoreErrors = 1;
        matlabbatch{batch_idx}.spm.tools.cat.estwrite.output.BIDS.BIDSno = 1;
        
        % Disable surface mapping as requested
        matlabbatch{batch_idx}.spm.tools.cat.estwrite.output.surface = 0;
        matlabbatch{batch_idx}.spm.tools.cat.estwrite.output.surf_measures = 0;
        
        % Atlas settings - keep existing ones
        matlabbatch{batch_idx}.spm.tools.cat.estwrite.output.ROImenu.atlases.neuromorphometrics = 1;
        matlabbatch{batch_idx}.spm.tools.cat.estwrite.output.ROImenu.atlases.lpba40 = 1;
        matlabbatch{batch_idx}.spm.tools.cat.estwrite.output.ROImenu.atlases.cobra = 0;
        matlabbatch{batch_idx}.spm.tools.cat.estwrite.output.ROImenu.atlases.hammers = 0;
        matlabbatch{batch_idx}.spm.tools.cat.estwrite.output.ROImenu.atlases.thalamus = 0;
        matlabbatch{batch_idx}.spm.tools.cat.estwrite.output.ROImenu.atlases.thalamic_nuclei = 0;   
        matlabbatch{batch_idx}.spm.tools.cat.estwrite.output.ROImenu.atlases.suit = 0;
        matlabbatch{batch_idx}.spm.tools.cat.estwrite.output.ROImenu.atlases.ibsr = 0;
        
        % Tissue segmentation outputs
        matlabbatch{batch_idx}.spm.tools.cat.estwrite.output.GM.native = 1;
        matlabbatch{batch_idx}.spm.tools.cat.estwrite.output.GM.warped = 1;  % Enable warped for analysis
        matlabbatch{batch_idx}.spm.tools.cat.estwrite.output.GM.mod = 0;
        matlabbatch{batch_idx}.spm.tools.cat.estwrite.output.GM.dartel = 0;
        matlabbatch{batch_idx}.spm.tools.cat.estwrite.output.WM.native = 1;
        matlabbatch{batch_idx}.spm.tools.cat.estwrite.output.WM.warped = 1;  % Enable warped for analysis
        matlabbatch{batch_idx}.spm.tools.cat.estwrite.output.WM.mod = 0;
        matlabbatch{batch_idx}.spm.tools.cat.estwrite.output.WM.dartel = 0;
        matlabbatch{batch_idx}.spm.tools.cat.estwrite.output.CSF.native = 1;
        matlabbatch{batch_idx}.spm.tools.cat.estwrite.output.CSF.warped = 0;
        matlabbatch{batch_idx}.spm.tools.cat.estwrite.output.CSF.mod = 0;
        matlabbatch{batch_idx}.spm.tools.cat.estwrite.output.CSF.dartel = 0;
        
        % Other outputs (disabled)
        matlabbatch{batch_idx}.spm.tools.cat.estwrite.output.ct.native = 0;
        matlabbatch{batch_idx}.spm.tools.cat.estwrite.output.ct.warped = 0;
        matlabbatch{batch_idx}.spm.tools.cat.estwrite.output.ct.dartel = 0;
        matlabbatch{batch_idx}.spm.tools.cat.estwrite.output.pp.native = 0;
        matlabbatch{batch_idx}.spm.tools.cat.estwrite.output.pp.warped = 0;
        matlabbatch{batch_idx}.spm.tools.cat.estwrite.output.pp.dartel = 0;
        matlabbatch{batch_idx}.spm.tools.cat.estwrite.output.WMH.native = 0;
        matlabbatch{batch_idx}.spm.tools.cat.estwrite.output.WMH.warped = 0;
        matlabbatch{batch_idx}.spm.tools.cat.estwrite.output.WMH.mod = 0;
        matlabbatch{batch_idx}.spm.tools.cat.estwrite.output.WMH.dartel = 0;
        matlabbatch{batch_idx}.spm.tools.cat.estwrite.output.SL.native = 0;
        matlabbatch{batch_idx}.spm.tools.cat.estwrite.output.SL.warped = 0;
        matlabbatch{batch_idx}.spm.tools.cat.estwrite.output.SL.mod = 0;
        matlabbatch{batch_idx}.spm.tools.cat.estwrite.output.SL.dartel = 0;
        matlabbatch{batch_idx}.spm.tools.cat.estwrite.output.TPMC.native = 0;
        matlabbatch{batch_idx}.spm.tools.cat.estwrite.output.TPMC.warped = 0;
        matlabbatch{batch_idx}.spm.tools.cat.estwrite.output.TPMC.mod = 0;
        matlabbatch{batch_idx}.spm.tools.cat.estwrite.output.TPMC.dartel = 0;
        matlabbatch{batch_idx}.spm.tools.cat.estwrite.output.atlas.native = 0;
        matlabbatch{batch_idx}.spm.tools.cat.estwrite.output.label.native = 1;
        matlabbatch{batch_idx}.spm.tools.cat.estwrite.output.label.warped = 0;
        matlabbatch{batch_idx}.spm.tools.cat.estwrite.output.label.dartel = 0;
        matlabbatch{batch_idx}.spm.tools.cat.estwrite.output.labelnative = 1;
        matlabbatch{batch_idx}.spm.tools.cat.estwrite.output.bias.warped = 1;
        matlabbatch{batch_idx}.spm.tools.cat.estwrite.output.las.native = 0;
        matlabbatch{batch_idx}.spm.tools.cat.estwrite.output.las.warped = 0;
        matlabbatch{batch_idx}.spm.tools.cat.estwrite.output.las.dartel = 0;
        matlabbatch{batch_idx}.spm.tools.cat.estwrite.output.jacobianwarped = 0;
        matlabbatch{batch_idx}.spm.tools.cat.estwrite.output.warps = [1 1];  % Forward and inverse deformation fields
        matlabbatch{batch_idx}.spm.tools.cat.estwrite.output.rmat = 0;
        
        batch_idx = batch_idx + 1;
    else
        disp('CAT12 processing already completed, skipping segmentation step.');
        spm('defaults', 'FMRI');
        spm_jobman('initcfg');
    end

    % Step 2: Coregister DTI to native T1 space and then normalize
    if dti_success && isfield(dti_struct, 'S0')
        disp('Processing DTI data...');
        
        % Create DTI normalized directory
        dti_norm_dir = fullfile(mri_dir, 'DTI');
        if ~exist(dti_norm_dir, 'dir'), mkdir(dti_norm_dir); end
        
        % Get T1 native space image for coregistration target
        t1_native = fullfile(mri_dir, ['m' t1wfiles(end).name]);
        if ~exist(t1_native, 'file')
            % If bias-corrected doesn't exist, use original
            t1_native = fullfile(newdir, t1wfiles(end).name);
        end
        
        % Coregister DTI to native T1 space and WRITE coregistered files
        dti_files = struct2cell(dti_struct);
        ref_image = dti_struct.S0;  % Use S0 as reference for coregistration
        source_images = dti_files(~strcmp(dti_files, ref_image));  % All other DTI files
        
        matlabbatch{batch_idx}.spm.spatial.coreg.estwrite.ref = {t1_native};
        matlabbatch{batch_idx}.spm.spatial.coreg.estwrite.source = {ref_image};
        matlabbatch{batch_idx}.spm.spatial.coreg.estwrite.other = source_images;
        matlabbatch{batch_idx}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
        matlabbatch{batch_idx}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
        matlabbatch{batch_idx}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
        matlabbatch{batch_idx}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
        matlabbatch{batch_idx}.spm.spatial.coreg.estwrite.roptions.interp = 4;
        matlabbatch{batch_idx}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
        matlabbatch{batch_idx}.spm.spatial.coreg.estwrite.roptions.mask = 0;
        matlabbatch{batch_idx}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';
        
        batch_idx = batch_idx + 1;
        
        % Copy coregistered DTI files to DTI output directory for easy access
        % (the r* files are created in original directories by estwrite, but let's also copy them)
        dti_output_dir = fullfile(newdir, 'DTI_coregistered');
        if ~exist(dti_output_dir, 'dir'), mkdir(dti_output_dir); end
        
        % Copy the coregistered files that were just created
        for i = 1:length(dti_files)
            [pth, nam, ext] = fileparts(dti_files{i});
            coregistered_file = fullfile(pth, ['r' nam ext]);
            if exist(coregistered_file, 'file')
                copyfile(coregistered_file, dti_output_dir);
            end
        end
        
        % Apply forward deformation to coregistered DTI files for normalization
        if ~isempty(y_files)
            matlabbatch{batch_idx}.spm.util.defs.comp{1}.def = {fullfile(y_files(1).folder, y_files(1).name)};
        else
            matlabbatch{batch_idx}.spm.util.defs.comp{1}.def(1) = cfg_dep('CAT12: Segmentation: Deformation Field', ...
                substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
                substruct('()',{1}, '.','def', '()',{':'}));
        end
        
        % Create list of coregistered DTI files (with 'r' prefix) for normalization
        coregistered_dti_files = cell(size(dti_files));
        for i = 1:length(dti_files)
            [pth, nam, ext] = fileparts(dti_files{i});
            coregistered_dti_files{i} = fullfile(pth, ['r' nam ext]);
        end
        
        matlabbatch{batch_idx}.spm.util.defs.out{1}.pull.fnames = coregistered_dti_files;
        matlabbatch{batch_idx}.spm.util.defs.out{1}.pull.savedir.saveusr = {dti_norm_dir};
        matlabbatch{batch_idx}.spm.util.defs.out{1}.pull.interp = 4;
        matlabbatch{batch_idx}.spm.util.defs.out{1}.pull.mask = 1;
        matlabbatch{batch_idx}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
        matlabbatch{batch_idx}.spm.util.defs.out{1}.pull.prefix = 'w';
        
        batch_idx = batch_idx + 1;
    end

    % Step 3: Coregister NODDI to native T1 space and then normalize
    if noddi_success && isfield(noddi_struct, 'File1')
        disp('Processing NODDI data...');
        
        % Create NODDI normalized directory
        noddi_norm_dir = fullfile(mri_dir, 'NODDI');
        if ~exist(noddi_norm_dir, 'dir'), mkdir(noddi_norm_dir); end
        
        % Get T1 native space image for coregistration target
        t1_native = fullfile(mri_dir, ['m' t1wfiles(end).name]);
        if ~exist(t1_native, 'file')
            t1_native = fullfile(newdir, t1wfiles(end).name);
        end
        
        % Coregister NODDI to native T1 space and WRITE coregistered files
        noddi_files = struct2cell(noddi_struct);
        ref_image = noddi_files{1};  % Use first NODDI file as reference
        source_images = noddi_files(2:end);  % All other NODDI files
        
        matlabbatch{batch_idx}.spm.spatial.coreg.estwrite.ref = {t1_native};
        matlabbatch{batch_idx}.spm.spatial.coreg.estwrite.source = {ref_image};
        matlabbatch{batch_idx}.spm.spatial.coreg.estwrite.other = source_images;
        matlabbatch{batch_idx}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
        matlabbatch{batch_idx}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
        matlabbatch{batch_idx}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
        matlabbatch{batch_idx}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
        matlabbatch{batch_idx}.spm.spatial.coreg.estwrite.roptions.interp = 4;
        matlabbatch{batch_idx}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
        matlabbatch{batch_idx}.spm.spatial.coreg.estwrite.roptions.mask = 0;
        matlabbatch{batch_idx}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';
        
        batch_idx = batch_idx + 1;
        
        % Copy coregistered NODDI files to NODDI output directory for easy access
        noddi_output_dir = fullfile(newdir, 'NODDI_coregistered');
        if ~exist(noddi_output_dir, 'dir'), mkdir(noddi_output_dir); end
        
        % Copy the coregistered files that were just created
        for i = 1:length(noddi_files)
            [pth, nam, ext] = fileparts(noddi_files{i});
            coregistered_file = fullfile(pth, ['r' nam ext]);
            if exist(coregistered_file, 'file')
                copyfile(coregistered_file, noddi_output_dir);
            end
        end
        
        % Apply forward deformation to coregistered NODDI files for normalization
        if ~isempty(y_files)
            matlabbatch{batch_idx}.spm.util.defs.comp{1}.def = {fullfile(y_files(1).folder, y_files(1).name)};
        else
            matlabbatch{batch_idx}.spm.util.defs.comp{1}.def(1) = cfg_dep('CAT12: Segmentation: Deformation Field', ...
                substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
                substruct('()',{1}, '.','def', '()',{':'}));
        end
        
        % Create list of coregistered NODDI files (with 'r' prefix) for normalization
        coregistered_noddi_files = cell(size(noddi_files));
        for i = 1:length(noddi_files)
            [pth, nam, ext] = fileparts(noddi_files{i});
            coregistered_noddi_files{i} = fullfile(pth, ['r' nam ext]);
        end
        
        matlabbatch{batch_idx}.spm.util.defs.out{1}.pull.fnames = coregistered_noddi_files;
        matlabbatch{batch_idx}.spm.util.defs.out{1}.pull.savedir.saveusr = {noddi_norm_dir};
        matlabbatch{batch_idx}.spm.util.defs.out{1}.pull.interp = 4;
        matlabbatch{batch_idx}.spm.util.defs.out{1}.pull.mask = 1;
        matlabbatch{batch_idx}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
        matlabbatch{batch_idx}.spm.util.defs.out{1}.pull.prefix = 'w';
        
        batch_idx = batch_idx + 1;
    end

    % Step 4: Coregister ASL to native T1 space and then normalize
    if asl_success && isfield(asl_struct, 'M0_masked')
        disp('Processing ASL data...');
        
        % Create ASL directories
        asl_norm_dir = fullfile(mri_dir, 'ASL');
        if ~exist(asl_norm_dir, 'dir'), mkdir(asl_norm_dir); end
        
        % Get T1 native space image for coregistration target
        t1_native = fullfile(mri_dir, ['m' t1wfiles(end).name]);
        if ~exist(t1_native, 'file')
            t1_native = fullfile(newdir, t1wfiles(end).name);
        end
        
        % Handle multiple ASL datasets
        all_coregistered_asl_files = {};
        
        for asl_idx = 1:length(asl_struct)
            asl_files_cell = struct2cell(asl_struct(asl_idx));
            asl_files_current = asl_files_cell(~cellfun('isempty', asl_files_cell));
            
            if ~isempty(asl_files_current)
                ref_image = asl_struct(asl_idx).M0_masked;  % Use M0 as reference
                source_images = asl_files_current(~strcmp(asl_files_current, ref_image));
                
                % Coregister ASL to native T1 space and WRITE coregistered files
                matlabbatch{batch_idx}.spm.spatial.coreg.estwrite.ref = {t1_native};
                matlabbatch{batch_idx}.spm.spatial.coreg.estwrite.source = {ref_image};
                matlabbatch{batch_idx}.spm.spatial.coreg.estwrite.other = source_images;
                matlabbatch{batch_idx}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
                matlabbatch{batch_idx}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
                matlabbatch{batch_idx}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
                matlabbatch{batch_idx}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
                matlabbatch{batch_idx}.spm.spatial.coreg.estwrite.roptions.interp = 4;
                matlabbatch{batch_idx}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
                matlabbatch{batch_idx}.spm.spatial.coreg.estwrite.roptions.mask = 0;
                matlabbatch{batch_idx}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';
                
                batch_idx = batch_idx + 1;
                
                % Copy coregistered ASL files to ASL output directory for easy access
                asl_output_dir = fullfile(newdir, 'ASL_coregistered');
                if ~exist(asl_output_dir, 'dir'), mkdir(asl_output_dir); end
                
                % Copy the coregistered files that were just created
                for i = 1:length(asl_files_current)
                    [pth, nam, ext] = fileparts(asl_files_current{i});
                    coregistered_file = fullfile(pth, ['r' nam ext]);
                    if exist(coregistered_file, 'file')
                        copyfile(coregistered_file, asl_output_dir);
                    end
                end
                
                % Create list of coregistered ASL files (with 'r' prefix) for normalization
                coregistered_asl_files = cell(size(asl_files_current));
                for i = 1:length(asl_files_current)
                    [pth, nam, ext] = fileparts(asl_files_current{i});
                    coregistered_asl_files{i} = fullfile(pth, ['r' nam ext]);
                end
                
                % Add to master list
                all_coregistered_asl_files = [all_coregistered_asl_files; coregistered_asl_files];
                
                % Apply forward deformation to coregistered ASL files for normalization
                if ~isempty(y_files)
                    matlabbatch{batch_idx}.spm.util.defs.comp{1}.def = {fullfile(y_files(1).folder, y_files(1).name)};
                else
                    matlabbatch{batch_idx}.spm.util.defs.comp{1}.def(1) = cfg_dep('CAT12: Segmentation: Deformation Field', ...
                        substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
                        substruct('()',{1}, '.','def', '()',{':'}));
                end
                
                matlabbatch{batch_idx}.spm.util.defs.out{1}.pull.fnames = coregistered_asl_files;
                matlabbatch{batch_idx}.spm.util.defs.out{1}.pull.savedir.saveusr = {asl_norm_dir};
                matlabbatch{batch_idx}.spm.util.defs.out{1}.pull.interp = 4;
                matlabbatch{batch_idx}.spm.util.defs.out{1}.pull.mask = 1;
                matlabbatch{batch_idx}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
                matlabbatch{batch_idx}.spm.util.defs.out{1}.pull.prefix = 'w';
                
                batch_idx = batch_idx + 1;
            end
        end
    end

    % Step 5: Apply inverse deformation to bring all label maps to native space
    disp('Processing label maps to native space...');
    
    for atlas_idx = 1:length(label_maps_files)
        current_atlas = label_maps_files{atlas_idx};
        
        % Special handling for hypothalamic atlas (might be .nii.gz)
        atlas_file = fullfile(label_map_path, current_atlas);
        
        if ~isempty(iy_files)
            matlabbatch{batch_idx}.spm.util.defs.comp{1}.def = {fullfile(iy_files(1).folder, iy_files(1).name)};
        else
            matlabbatch{batch_idx}.spm.util.defs.comp{1}.def(1) = cfg_dep('CAT12: Segmentation: Inverse Deformation Field', ...
                substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
                substruct('()',{1}, '.','invdef', '()',{':'}));
        end
        
        matlabbatch{batch_idx}.spm.util.defs.out{1}.pull.fnames = {atlas_file};
        matlabbatch{batch_idx}.spm.util.defs.out{1}.pull.savedir.saveusr = {newdir};
        matlabbatch{batch_idx}.spm.util.defs.out{1}.pull.interp = 0;  % Nearest neighbor for label maps
        matlabbatch{batch_idx}.spm.util.defs.out{1}.pull.mask = 1;
        matlabbatch{batch_idx}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
        matlabbatch{batch_idx}.spm.util.defs.out{1}.pull.prefix = 'w';
        
        batch_idx = batch_idx + 1;
    end

    % Step 6: Apply forward deformation to native tissue segmentations to create normalized versions
    disp('Creating normalized tissue segmentations (wrp*)...');
    
    % Find native tissue segmentations
    gm_native_files = dir(fullfile(mri_dir, 'p1*.nii'));
    wm_native_files = dir(fullfile(mri_dir, 'p2*.nii'));
    csf_native_files = dir(fullfile(mri_dir, 'p3*.nii'));
    
    if ~isempty(gm_native_files) || ~isempty(wm_native_files) || ~isempty(csf_native_files)
        % Collect all native segmentation files
        native_seg_files = {};
        if ~isempty(gm_native_files)
            native_seg_files{end+1} = fullfile(gm_native_files(1).folder, gm_native_files(1).name);
        end
        if ~isempty(wm_native_files)
            native_seg_files{end+1} = fullfile(wm_native_files(1).folder, wm_native_files(1).name);
        end
        if ~isempty(csf_native_files)
            native_seg_files{end+1} = fullfile(csf_native_files(1).folder, csf_native_files(1).name);
        end
        
        % Apply forward deformation to native segmentations
        if ~isempty(y_files)
            matlabbatch{batch_idx}.spm.util.defs.comp{1}.def = {fullfile(y_files(1).folder, y_files(1).name)};
        else
            matlabbatch{batch_idx}.spm.util.defs.comp{1}.def(1) = cfg_dep('CAT12: Segmentation: Deformation Field', ...
                substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), ...
                substruct('()',{1}, '.','def', '()',{':'}));
        end
        
        matlabbatch{batch_idx}.spm.util.defs.out{1}.pull.fnames = native_seg_files;
        matlabbatch{batch_idx}.spm.util.defs.out{1}.pull.savedir.saveusr = {mri_dir};
        matlabbatch{batch_idx}.spm.util.defs.out{1}.pull.interp = 4;  % Trilinear for tissue segmentations
        matlabbatch{batch_idx}.spm.util.defs.out{1}.pull.mask = 1;
        matlabbatch{batch_idx}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
        matlabbatch{batch_idx}.spm.util.defs.out{1}.pull.prefix = 'wr';  % Warped native segmentations
        
        batch_idx = batch_idx + 1;
    end

    % Run all batch jobs
    disp('Running SPM batch jobs...');
    spm_jobman('run', matlabbatch);
    clear matlabbatch;
    
    disp('Processing complete!');
    
    % Print summary of outputs
    disp('=== Processing Summary ===');
    disp(['Output directory: ' newdir]);
    disp(['MRI outputs: ' mri_dir]);
    
    if dti_success
        disp(['DTI coregistered (r* prefix): Original DTI directories + ' fullfile(newdir, 'DTI_coregistered')]);
        disp(['DTI normalized (wr* prefix): ' fullfile(mri_dir, 'DTI')]);
    end
    if noddi_success
        disp(['NODDI coregistered (r* prefix): Original NODDI directory + ' fullfile(newdir, 'NODDI_coregistered')]);
        disp(['NODDI normalized (wr* prefix): ' fullfile(mri_dir, 'NODDI')]);
    end
    if asl_success
        disp(['ASL coregistered (r* prefix): Original ASL directories + ' fullfile(newdir, 'ASL_coregistered')]);
        disp(['ASL normalized (wr* prefix): ' fullfile(mri_dir, 'ASL')]);
    end
    
    disp('Native space atlases (w* prefix) available in output directory');
    disp('Native tissue segmentations (p1*, p2*, p3*) available in mri directory');
    disp('Normalized tissue segmentations (wrp1*, wrp2*, wrp3*) available in mri directory');
    disp('CAT12 template segmentations (mwp1*, mwp2*) available in mri directory');
    disp('');
    disp('For native space analysis: Use r* files with w* atlas files');
    disp('For template space analysis: Use wr* files with original template atlases');
    disp('For tissue reference regions: Use p* (native) or wrp* (normalized) segmentations');

end