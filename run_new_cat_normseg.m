function run_newcat_normseg(base_dir)

    % Add SPM12 and CAT12 to the MATLAB path
    addpath('/isilon/datalake/riipl/scratch/ADRC/Hellcat-12.9/libs/spm12/spm12/');
    addpath('/isilon/datalake/riipl/scratch/ADRC/Hellcat-12.9/');

    % atlas aavl3
    label_map_path='/isilon/datalake/riipl/scratch/ADRC/Hellcat-12.9/libs/spm12/spm12/toolbox/cat12/templates_MNI152NLin2009cAsym/';
    label_maps_files = {'neuromorphometrics.nii','aal3.nii','cobra.nii','hammers.nii','lpba40.nii','ibsr.nii','JHU.nii','brainmask_T1.nii'};

    % Define the base directory of test subject
    %base_dir = '/isilon/datalake/riipl/scratch/ADRC/Hellcat-12.9/testdata/3MCJO002401_125300350_20231002/';

    % Create the new directory
    newdir = fullfile(base_dir, 'nifti', 'cat12_v2560');
    if ~exist(newdir, 'dir')
        mkdir(newdir);
    end

    % Find the T1w image
    t1wfiles = dir(fullfile([base_dir, '/nifti/'], '3*-tfl3d116ns.nii'));
    if isempty(t1wfiles)
        error('No T1w image found');
    end

    % Copy the T1w image to the new directory
    copyfile(fullfile(t1wfiles(1).folder, t1wfiles(1).name), newdir, 'f');

    % Check if a file matching the pattern "iy_*-tfl3d116ns.nii" is present in the subject's nifti/cat12_2550/mri folder
    mri_dir = fullfile(newdir, 'mri');
    y_files = dir(fullfile(mri_dir, 'y_*-tfl3d116ns.nii'));
    iy_files = dir(fullfile(mri_dir, 'iy_*-tfl3d116ns.nii'));

    disp(sprintf("Subject dir %s", base_dir));
    disp(sprintf("MRI dir %s", mri_dir));
    disp(sprintf("Input T1 %s", fullfile(newdir, t1wfiles(1).name)));

    [asl_struct, ~] = getASLFiles(base_dir);
    [dti_struct, ~] = getDTIFiles(base_dir);
    [noddi_struct, ~] = getNODDIFiles(base_dir);

    if ~isempty(y_files)
        disp('File matching the pattern "y_*-tfl3d116ns.nii" is present in the mri folder.');
    else
        disp('No file matching the pattern "y_*-tfl3d116ns.nii" found in the mri folder.');
        % Perform the segmentation and normalization using SPM12 and CAT12
        % This is a placeholder - replace with your actual code
        spm('defaults', 'FMRI');
        spm_jobman('initcfg');

        matlabbatch{1}.spm.tools.cat.estwrite.data = {fullfile(newdir, t1wfiles(1).name)};
        matlabbatch{1}.spm.tools.cat.estwrite.data_wmh = {''};
        matlabbatch{1}.spm.tools.cat.estwrite.nproc = 8;
        matlabbatch{1}.spm.tools.cat.estwrite.useprior = '';
        matlabbatch{1}.spm.tools.cat.estwrite.opts.tpm = {'/isilon/datalake/riipl/scratch/ADRC/Hellcat-12.9/libs/spm12/spm12/tpm/TPM.nii'};
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
        matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.shooting.shootingtpm = {'/isilon/datalake/riipl/scratch/ADRC/Hellcat-12.9/libs/spm12/spm12/toolbox/cat12/templates_MNI152NLin2009cAsym/Template_0_GS.nii'};
        matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.shooting.regstr = 0.5;
        matlabbatch{1}.spm.tools.cat.estwrite.extopts.vox = 1.5;
        matlabbatch{1}.spm.tools.cat.estwrite.extopts.bb = 12;
        matlabbatch{1}.spm.tools.cat.estwrite.extopts.SRP = 22;
        matlabbatch{1}.spm.tools.cat.estwrite.extopts.ignoreErrors = 1;
        matlabbatch{1}.spm.tools.cat.estwrite.output.BIDS.BIDSno = 1;
        matlabbatch{1}.spm.tools.cat.estwrite.output.surface = 0;
        matlabbatch{1}.spm.tools.cat.estwrite.output.surf_measures = 0;
        matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.neuromorphometrics = 0;
        matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.lpba40 = 0;
        matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.cobra = 0;
        matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.hammers = 0;
        matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.thalamus = 0;
        matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.thalamic_nuclei = 0;
        matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.suit = 0;
        matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.ibsr = 0;
        matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.atlases.ownatlas = {''};
        matlabbatch{1}.spm.tools.cat.estwrite.output.GM.native = 1;
        matlabbatch{1}.spm.tools.cat.estwrite.output.GM.mod = 0;
        matlabbatch{1}.spm.tools.cat.estwrite.output.GM.dartel = 0;
        matlabbatch{1}.spm.tools.cat.estwrite.output.WM.native = 1;
        matlabbatch{1}.spm.tools.cat.estwrite.output.WM.mod = 0;
        matlabbatch{1}.spm.tools.cat.estwrite.output.WM.dartel = 0;
        matlabbatch{1}.spm.tools.cat.estwrite.output.CSF.native = 0;
        matlabbatch{1}.spm.tools.cat.estwrite.output.CSF.warped = 0;
        matlabbatch{1}.spm.tools.cat.estwrite.output.CSF.mod = 0;
        matlabbatch{1}.spm.tools.cat.estwrite.output.CSF.dartel = 0;
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
        matlabbatch{1}.spm.tools.cat.estwrite.output.atlas.native = 1;
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
    end

    if ~isempty(iy_files)
        spm('defaults', 'FMRI');
        spm_jobman('initcfg');
        i=1;
    else
        i=2;
    end

    if ~isempty(iy_files)
        matlabbatch{i}.spm.util.defs.comp{1}.def = {fullfile(iy_files(1).folder, iy_files(1).name)};
    else
        matlabbatch{i}.spm.util.defs.comp{1}.def(1) = cfg_dep('CAT12: Segmentation: Inverse Deformation Field', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','invdef', '()',{':'}));
    end
    matlabbatch{i}.spm.util.defs.out{1}.pull.fnames = fullfile(label_map_path, label_maps_files(1));
    matlabbatch{i}.spm.util.defs.out{1}.pull.savedir.saveusr = {newdir};
    matlabbatch{i}.spm.util.defs.out{1}.pull.interp = 0;
    matlabbatch{i}.spm.util.defs.out{1}.pull.mask = 1;
    matlabbatch{i}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
    matlabbatch{i}.spm.util.defs.out{1}.pull.prefix = 'w';
    i=i+1;

    if ~isempty(iy_files)
        matlabbatch{i}.spm.util.defs.comp{1}.def = {fullfile(iy_files(1).folder, iy_files(1).name)};
    else
        matlabbatch{i}.spm.util.defs.comp{1}.def(1) = cfg_dep('CAT12: Segmentation: Inverse Deformation Field', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','invdef', '()',{':'}));
    end
    matlabbatch{i}.spm.util.defs.out{1}.pull.fnames = fullfile(label_map_path, label_maps_files(2));
    matlabbatch{i}.spm.util.defs.out{1}.pull.savedir.saveusr = {newdir};
    matlabbatch{i}.spm.util.defs.out{1}.pull.interp = 0;
    matlabbatch{i}.spm.util.defs.out{1}.pull.mask = 1;
    matlabbatch{i}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
    matlabbatch{i}.spm.util.defs.out{1}.pull.prefix = 'w';
    i=i+1;

    if ~isempty(iy_files)
        matlabbatch{i}.spm.util.defs.comp{1}.def = {fullfile(iy_files(1).folder, iy_files(1).name)};
    else
        matlabbatch{i}.spm.util.defs.comp{1}.def(1) = cfg_dep('CAT12: Segmentation: Inverse Deformation Field', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','invdef', '()',{':'}));
    end
    matlabbatch{i}.spm.util.defs.out{1}.pull.fnames = fullfile(label_map_path, label_maps_files(3));
    matlabbatch{i}.spm.util.defs.out{1}.pull.savedir.saveusr = {newdir};
    matlabbatch{i}.spm.util.defs.out{1}.pull.interp = 0;
    matlabbatch{i}.spm.util.defs.out{1}.pull.mask = 1;
    matlabbatch{i}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
    matlabbatch{i}.spm.util.defs.out{1}.pull.prefix = 'w';
    i=i+1;

    if ~isempty(iy_files)
        matlabbatch{i}.spm.util.defs.comp{1}.def = {fullfile(iy_files(1).folder, iy_files(1).name)};
    else
        matlabbatch{i}.spm.util.defs.comp{1}.def(1) = cfg_dep('CAT12: Segmentation: Inverse Deformation Field', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','invdef', '()',{':'}));
    end
    matlabbatch{i}.spm.util.defs.out{1}.pull.fnames = fullfile(label_map_path, label_maps_files(4));
    matlabbatch{i}.spm.util.defs.out{1}.pull.savedir.saveusr = {newdir};
    matlabbatch{i}.spm.util.defs.out{1}.pull.interp = 0;
    matlabbatch{i}.spm.util.defs.out{1}.pull.mask = 1;
    matlabbatch{i}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
    matlabbatch{i}.spm.util.defs.out{1}.pull.prefix = 'w';
    i=i+1;

    if ~isempty(iy_files)
        matlabbatch{i}.spm.util.defs.comp{1}.def = {fullfile(iy_files(1).folder, iy_files(1).name)};
    else
        matlabbatch{i}.spm.util.defs.comp{1}.def(1) = cfg_dep('CAT12: Segmentation: Inverse Deformation Field', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','invdef', '()',{':'}));
    end
    matlabbatch{i}.spm.util.defs.out{1}.pull.fnames = fullfile(label_map_path, label_maps_files(5));
    matlabbatch{i}.spm.util.defs.out{1}.pull.savedir.saveusr = {newdir};
    matlabbatch{i}.spm.util.defs.out{1}.pull.interp = 0;
    matlabbatch{i}.spm.util.defs.out{1}.pull.mask = 1;
    matlabbatch{i}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
    matlabbatch{i}.spm.util.defs.out{1}.pull.prefix = 'w';
    i=i+1;

    if ~isempty(iy_files)
        matlabbatch{i}.spm.util.defs.comp{1}.def = {fullfile(iy_files(1).folder, iy_files(1).name)};
    else
        matlabbatch{i}.spm.util.defs.comp{1}.def(1) = cfg_dep('CAT12: Segmentation: Inverse Deformation Field', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','invdef', '()',{':'}));
    end
    matlabbatch{i}.spm.util.defs.out{1}.pull.fnames = fullfile(label_map_path, label_maps_files(6));
    matlabbatch{i}.spm.util.defs.out{1}.pull.savedir.saveusr = {newdir};
    matlabbatch{i}.spm.util.defs.out{1}.pull.interp = 0;
    matlabbatch{i}.spm.util.defs.out{1}.pull.mask = 1;
    matlabbatch{i}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
    matlabbatch{i}.spm.util.defs.out{1}.pull.prefix = 'w';
    i=i+1;

    if ~isempty(iy_files)
        matlabbatch{i}.spm.util.defs.comp{1}.def = {fullfile(iy_files(1).folder, iy_files(1).name)};
    else
        matlabbatch{i}.spm.util.defs.comp{1}.def(1) = cfg_dep('CAT12: Segmentation: Inverse Deformation Field', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','invdef', '()',{':'}));
    end
    matlabbatch{i}.spm.util.defs.out{1}.pull.fnames = fullfile(label_map_path, label_maps_files(7));
    matlabbatch{i}.spm.util.defs.out{1}.pull.savedir.saveusr = {newdir};
    matlabbatch{i}.spm.util.defs.out{1}.pull.interp = 0;
    matlabbatch{i}.spm.util.defs.out{1}.pull.mask = 1;
    matlabbatch{i}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
    matlabbatch{i}.spm.util.defs.out{1}.pull.prefix = 'w';
    i=i+1;

    if ~isempty(iy_files)
        matlabbatch{i}.spm.util.defs.comp{1}.def = {fullfile(iy_files(1).folder, iy_files(1).name)};
    else
        matlabbatch{i}.spm.util.defs.comp{1}.def(1) = cfg_dep('CAT12: Segmentation: Inverse Deformation Field', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','invdef', '()',{':'}));
    end
    matlabbatch{i}.spm.util.defs.out{1}.pull.fnames = fullfile(label_map_path, label_maps_files(8));
    matlabbatch{i}.spm.util.defs.out{1}.pull.savedir.saveusr = {newdir};
    matlabbatch{i}.spm.util.defs.out{1}.pull.interp = 0;
    matlabbatch{i}.spm.util.defs.out{1}.pull.mask = 1;
    matlabbatch{i}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
    matlabbatch{i}.spm.util.defs.out{1}.pull.prefix = 'w';
    
    if isfield(dti_struct, 'S0')
        i = i + 1;

        dti_dir = fullfile(mri_dir, 'DTI');
        if ~exist(dti_dir, 'dir')
            mkdir(dti_dir);
        end
        
        dti_files = struct2cell(dti_struct);
        dti_files = dti_files(:);

        if ~isempty(y_files)
            matlabbatch{i}.spm.util.defs.comp{1}.def = {fullfile(y_files(1).folder, y_files(1).name)};
        else
            matlabbatch{i}.spm.util.defs.comp{1}.def(1) = cfg_dep('CAT12: Segmentation: Deformation Field', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','invdef', '()',{':'}));
        end
        matlabbatch{i}.spm.util.defs.out{1}.pull.fnames = dti_files;
        matlabbatch{i}.spm.util.defs.out{1}.pull.savedir.saveusr = {dti_dir};
        matlabbatch{i}.spm.util.defs.out{1}.pull.interp = 4;
        matlabbatch{i}.spm.util.defs.out{1}.pull.mask = 1;
        matlabbatch{i}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
        matlabbatch{i}.spm.util.defs.out{1}.pull.prefix = 'w';
    end

    if isfield(asl_struct, 'M0_masked')
        i = i + 1;

        asl_dir = fullfile(mri_dir, 'ASL');
        if ~exist(asl_dir, 'dir')
            mkdir(asl_dir);
        end

        asl_files = struct2cell(asl_struct);
        asl_files = asl_files(:);

        if ~isempty(y_files)
            matlabbatch{i}.spm.util.defs.comp{1}.def = {fullfile(y_files(1).folder, y_files(1).name)};
        else
            matlabbatch{i}.spm.util.defs.comp{1}.def(1) = cfg_dep('CAT12: Segmentation: Deformation Field', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','invdef', '()',{':'}));
        end
        matlabbatch{i}.spm.util.defs.out{1}.pull.fnames = asl_files;
        matlabbatch{i}.spm.util.defs.out{1}.pull.savedir.saveusr = {asl_dir};
        matlabbatch{i}.spm.util.defs.out{1}.pull.interp = 4;
        matlabbatch{i}.spm.util.defs.out{1}.pull.mask = 1;
        matlabbatch{i}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
        matlabbatch{i}.spm.util.defs.out{1}.pull.prefix = 'w';
    end

    if isfield(noddi_struct, 'File1')
        i = i + 1;

        noddi_dir = fullfile(mri_dir, 'NODDI');
        if ~exist(noddi_dir, 'dir')
            mkdir(noddi_dir);
        end

        noddi_files = struct2cell(noddi_struct);
        noddi_files = noddi_files(:);

        if ~isempty(y_files)
            matlabbatch{i}.spm.util.defs.comp{1}.def = {fullfile(y_files(1).folder, y_files(1).name)};
        else
            matlabbatch{i}.spm.util.defs.comp{1}.def(1) = cfg_dep('CAT12: Segmentation: Deformation Field', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','invdef', '()',{':'}));
        end
        matlabbatch{i}.spm.util.defs.out{1}.pull.fnames = noddi_files;
        matlabbatch{i}.spm.util.defs.out{1}.pull.savedir.saveusr = {noddi_dir};
        matlabbatch{i}.spm.util.defs.out{1}.pull.interp = 4;
        matlabbatch{i}.spm.util.defs.out{1}.pull.mask = 1;
        matlabbatch{i}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
        matlabbatch{i}.spm.util.defs.out{1}.pull.prefix = 'w';
    end

    spm_jobman('run', matlabbatch);
    clear matlabbatch % clear matlabbatch

end