function [aslFiles, success] = getASLFiles(inputDataFolder)
    % Initialize output
    aslFiles = struct();
    success = true;
    
    % Find all ASL folders
    aslFolderPaths = dir(fullfile(inputDataFolder, 'ASL_S*'));
    aslFolderPaths = aslFolderPaths([aslFolderPaths.isdir]); % filter out non-directories
    
    % Check if any ASL folder exists
    if isempty(aslFolderPaths)
        disp('No ASL folder found');
        aslFiles = [];
        success = false;
        return;
    end
    
    % Iterate over each ASL folder
    for j = 1:length(aslFolderPaths)
        aslFolderPath = fullfile(aslFolderPaths(j).folder, aslFolderPaths(j).name);
        
        % Find the M0_masked file
        m0Files = dir(fullfile(aslFolderPath, 'M*_M0_masked.nii'));
        if isempty(m0Files)
            disp(['M0_masked file not found in ' aslFolderPath]);
            aslFiles = [];
            success = false;
            return;
        end
        
        % Add the M0_masked file to the struct
        aslFiles(j).M0_masked = fullfile(m0Files(1).folder, m0Files(1).name);
        
        % Find all CBF files (broad search)
        allCBFFiles = dir(fullfile(aslFolderPath, '*_CBF*.nii'));
        allCVRFiles = dir(fullfile(aslFolderPath, '*_CVR*.nii'));
        
        % Filter to keep only files starting with '3' or uppercase letters (A-Z)
        % This excludes files starting with lowercase like r, m, w
        otherFiles = [];
        for k = 1:length(allCBFFiles)
            firstChar = allCBFFiles(k).name(1);
            % Check if first character is '3' or uppercase letter
            if firstChar == '3' || (firstChar >= 'A' && firstChar <= 'Z')
                % Also exclude the M0_masked file if it matches
                if ~strcmp(allCBFFiles(k).name, m0Files(1).name)
                    otherFiles = [otherFiles; allCBFFiles(k)];
                end
            end
        end

        for k = 1:length(allCVRFiles)
            firstChar = allCVRFiles(k).name(1);
            % Check if first character is '3' or uppercase letter
            if firstChar == '3' || (firstChar >= 'A' && firstChar <= 'Z')
                % Also exclude the M0_masked file if it matches
                if ~strcmp(allCVRFiles(k).name, m0Files(1).name)
                    otherFiles = [otherFiles; allCVRFiles(k)];
                end
            end
        end
        
        % Add the filtered ASL files to the struct
        for i = 1:length(otherFiles)
            aslFiles(j).(['File' num2str(i)]) = fullfile(otherFiles(i).folder, otherFiles(i).name);
        end
    end
end
