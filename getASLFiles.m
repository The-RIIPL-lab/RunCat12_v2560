function [aslFiles, success] = getASLFiles(inputDataFolder)
    % Initialize output
    aslFiles = struct();
    success = true;
    
    % Find all ASL folders
    aslFolderPaths = dir(fullfile(inputDataFolder, 'ASL_S*'));
    aslFolderPaths = aslFolderPaths([aslFolderPaths.isdir]); % filter out non-directories
    
    % Check if any ASL folder exists
    if isempty(aslFolderPaths)
        error('No ASL folder found');
        success = false;
        return;
    end
    
    % Iterate over each ASL folder
    for j = 1:length(aslFolderPaths)
        aslFolderPath = fullfile(aslFolderPaths(j).folder, aslFolderPaths(j).name);
        
        % Find the M0_masked file
        m0Files = dir(fullfile(aslFolderPath, '3*_M0_masked.nii'));
        if isempty(m0Files)
            error(['M0_masked file not found in ' aslFolderPath]);
            success = false;
            return;
        end
        
        % Add the M0_masked file to the struct
        aslFiles(j).M0_masked = fullfile(m0Files(1).folder, m0Files(1).name);
        
        % Find the other ASL files
        otherFiles = dir(fullfile(aslFolderPath, '3*._CBF*.nii'));
        
        % Filter out the M0_masked file
        otherFiles = otherFiles(~strcmp({otherFiles.name}, m0Files(1).name));
        
        % Add the other ASL files to the struct
        for i = 1:length(otherFiles)
            aslFiles(j).(['File' num2str(i)]) = fullfile(otherFiles(i).folder, otherFiles(i).name);
        end
    end
end