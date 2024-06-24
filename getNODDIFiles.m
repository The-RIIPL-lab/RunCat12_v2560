function [noddiFiles, success] = getNODDIFiles(inputDataFolder)
    % Initialize output
    noddiFiles = struct();
    success = true;
    
    % Construct the NODDI folder path
    noddiFolderPath = fullfile(inputDataFolder, 'NODDI');
    
    % Check if the NODDI folder exists
    if ~isfolder(noddiFolderPath)
        disp('NODDI folder not found');
        noddiFiles=[];
        success = false;
        return;
    end
    
    % Find the NODDI files
    noddiFilesList = dir(fullfile(noddiFolderPath, '3*-NODDI_*.nii'));
    
    if isempty(noddiFilesList)
        disp('NODDI files not found');
        noddiFiles=[];
        success = false;
        return;
    end
    
    % Add the NODDI files to the struct
    for i = 1:length(noddiFilesList)
        noddiFiles.(['File' num2str(i)]) = fullfile(noddiFilesList(i).folder, noddiFilesList(i).name);
    end
end