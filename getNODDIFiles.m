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
    
    % Find all NODDI files (broad search)
    allNoddiFiles = dir(fullfile(noddiFolderPath, '*-NODDI_*.nii'));
    
    % Filter to keep only files starting with '3' or uppercase letters
    noddiFilesList = [];
    for k = 1:length(allNoddiFiles)
        firstChar = allNoddiFiles(k).name(1);
        if firstChar == '3' || (firstChar >= 'A' && firstChar <= 'Z')
            noddiFilesList = [noddiFilesList; allNoddiFiles(k)];
        end
    end
    
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