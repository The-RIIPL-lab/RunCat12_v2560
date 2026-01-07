function [dtiFiles, success] = getDTIFiles(inputDataFolder)
    % Initialize output
    dtiFiles = struct();
    success = true;
    
    % Construct the DTI folder path
    dtiFolderPath = fullfile(inputDataFolder, 'DTI');
    
    % Check if the DTI folder exists
    if ~isfolder(dtiFolderPath)
        disp('DTI folder not found');
        
        % Check for a folder that matches the pattern DTI_S####
        dtiPatternPath = fullfile(inputDataFolder, 'DTI_S*');
        dtiPatternFiles = dir(dtiPatternPath);
        
        if isempty(dtiPatternFiles)
            disp('No DTI_S#### folder found');
            dtiFiles=[];
            success = false;
            return;
        else
            dtiFolderPath = fullfile(inputDataFolder, dtiPatternFiles(1).name);
            disp(['Using ', dtiFolderPath, ' as DTI folder']);
        end
    end
    
    % Find all S0 files (broad search)
    allS0Files = dir(fullfile(dtiFolderPath, '*_ECC_S0.nii'));
    
    % Filter to keep only files starting with '3' or uppercase letters
    s0Files = [];
    for k = 1:length(allS0Files)
        firstChar = allS0Files(k).name(1);
        if firstChar == '3' || (firstChar >= 'A' && firstChar <= 'Z')
            s0Files = [s0Files; allS0Files(k)];
        end
    end
    
    if isempty(s0Files)
        disp('S0 file not found');
        dtiFiles=[];
        success = false;
        return;
    end
    
    % Add the S0 file to the struct
    dtiFiles.S0 = fullfile(s0Files(1).folder, s0Files(1).name);
    
    % Find all DTI files (broad search)
    allDTIFiles = dir(fullfile(dtiFolderPath, '*_ECC_*.nii'));
    
    % Filter to keep only files starting with '3' or uppercase letters
    % and exclude the S0 file
    otherFiles = [];
    for k = 1:length(allDTIFiles)
        firstChar = allDTIFiles(k).name(1);
        if firstChar == '3' || (firstChar >= 'A' && firstChar <= 'Z')
            % Exclude the S0 file
            if ~strcmp(allDTIFiles(k).name, s0Files(1).name)
                otherFiles = [otherFiles; allDTIFiles(k)];
            end
        end
    end
    
    % Add the other DTI files to the struct
    for i = 1:length(otherFiles)
        dtiFiles.(['File' num2str(i)]) = fullfile(otherFiles(i).folder, otherFiles(i).name);
    end
end