% This script takes .mat spectra taken at the same field value and averages
% them together. It saves the averaged data to a specified folder.


% Define the folder containing the .mat files
folderPath = 'P6_ETO_MIR_m';
averagedFolderPath = fullfile(folderPath, 'Averaged_Data');
% Get a list of all .mat files in the folder
files = dir(fullfile(folderPath, '*.mat'));

% Initialize a structure to hold groups of files based on T values
fileGroups = containers.Map();

% Loop through each file to group them based on the T value in the filename
for k = 1:length(files)
    filename = files(k).name;
    % Extract the part of the filename with the T value
    tIndex = find(filename == 'T', 1, 'last'); % Last occurrence of 'T'
    underscoreIndex = find(filename(1:tIndex) == '_', 1, 'last'); % Last occurrence of '_' before 'T'
    tValue = filename(underscoreIndex+1:tIndex+3); % Including 'T' and two digits after it

    % Add the file to the appropriate group
    if isKey(fileGroups, tValue)
        fileGroups(tValue) = [fileGroups(tValue), {filename}];
    else
        fileGroups(tValue) = {filename};
    end
end

% Process each group of files
groupKeys = keys(fileGroups);
for i = 1:length(groupKeys)
    tValue = groupKeys{i};
    filenames = fileGroups(tValue);

    % Initialize a matrix to accumulate the data
    accumulatedData = [];
    numFiles = length(filenames);

    % Load and accumulate the data from each file in the group
    for j = 1:numFiles
        filePath = fullfile(folderPath, filenames{j});
        loadedData = load(filePath); % This loads all variables into a struct

        % Assuming your data variable is named 'data' inside the .mat file
        % Adapt this if your variable name is different
        if isfield(loadedData, 'data')
            data = loadedData.data;
            if isempty(accumulatedData)
                accumulatedData = data;
            else
                accumulatedData(:, 2) = accumulatedData(:, 2) + data(:, 2);
            end
        else
            warning(['Data variable not found in file: ', filenames{j}]);
        end
    end

    % Average the second column
    accumulatedData(:, 2) = accumulatedData(:, 2) / numFiles;

    % Save the averaged data to a new .mat file
    newFilename = [tValue, '_ave.mat'];
    save(fullfile(averagedFolderPath, newFilename), 'accumulatedData');
end

disp('Averaging complete. Averaged files saved.');