function [dtiFiles, success] = getDTIFiles(inputDataFolder)
    % Initialize output
    dtiFiles = struct();
    success = true;
    
    % Construct the DTI folder path
    dtiFolderPath = fullfile(inputDataFolder, 'DTI');
    
    % Check if the DTI folder exists
    if ~isfolder(dtiFolderPath)
        error('DTI folder not found');
        success = false;
        return;
    end
    
    % Find the S0 file
    s0Files = dir(fullfile(dtiFolderPath, '3*_ECC_S0.nii'));
    if isempty(s0Files)
        error('S0 file not found');
        success = false;
        return;
    end
    
    % Add the S0 file to the struct
    dtiFiles.S0 = fullfile(s0Files(1).folder, s0Files(1).name);
    
    % Find the other DTI files
    otherFiles = dir(fullfile(dtiFolderPath, '3*_ECC_*.nii'));
    
    % Filter out the S0 file
    otherFiles = otherFiles(~strcmp({otherFiles.name}, s0Files(1).name));
    
    % Add the other DTI files to the struct
    for i = 1:length(otherFiles)
        dtiFiles.(['File' num2str(i)]) = fullfile(otherFiles(i).folder, otherFiles(i).name);
    end
end