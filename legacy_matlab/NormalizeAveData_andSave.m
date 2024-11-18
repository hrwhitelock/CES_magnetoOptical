%This script takes the averaged .mat data, normalizes it with max field average normalization, 
%then  combines it into a single cell array with wavenumber in first column
%and intensity data in subsequent columns. The first row has the field
%value for each spectrum

% Make sure to specify the folder to save to at the end of the script

folderPath = './data/IR_MAgLab_Feb2024/Raw_mat/P6_ETO_FIR_m'; % Select the folder path of the data
% averagedFolderPath = fullfile(folderPath, 'Averaged_Data'); 
folderTitle = strrep(folderPath, '_', '\_');

%=========================================================================
% Process files
unsortedfiles = dir(fullfile(folderPath, '*.mat'));
numericValues = zeros(length(unsortedfiles), 1);

for k = 1:length(unsortedfiles)
    filename = unsortedfiles(k).name;
    tIndex = find(filename == 'T', 1, 'last');
    numericStr = filename(1:tIndex-1);
    numericValues(k) = str2double(numericStr);
end

[~, sortedIdx] = sort(numericValues);
files = unsortedfiles(sortedIdx);

%==========================================================================
% Create labels
labels = cell(1, length(files) + 1);
labels{1} = 'Wavenumbers';  % Set the first label for wavenumbers
for k = 1:length(files)
    filename = files(k).name;
    numericPart = regexp(filename, '\d*\.?\d+', 'match');  % Extract full numeric values including decimals
    if ~isempty(numericPart)
        labels{k + 1} = numericPart{1};  % Use just the numeric part as the label
    else
        labels{k + 1} = ['Data_' num2str(k)];  % Provide a default label if no numeric part is found
    end
end

%==========================================================================
% Load and combine data
firstFile = load(fullfile(files(1).folder, files(1).name));
wavenumbers = firstFile.data(:, 1); 
combined_array = wavenumbers;

for k = 1:length(files)
    dataFile = load(fullfile(files(k).folder, files(k).name));
    combined_array = [combined_array, dataFile.data(:, 2)];
end

% Normalize and adjust data
minValue = min(min(combined_array(:, 2:end)));
combined_array(:, 2:end) = combined_array(:, 2:end) + abs(minValue) + 0.0000001;
maxValues = max(combined_array(:, 2:end), [], 2);
normalized_combined_array = [combined_array(:, 1), 1 - (combined_array(:, 2:end) ./ maxValues)];

%==========================================================================
%==========================================================================
% Convert the normalized_combined_array to a cell array and insert labels
finalData = num2cell(normalized_combined_array);
finalData = [labels; finalData];  % Prepend the labels row to the top of the data array

% Display the finalData to verify structure
%disp('Final data array with labels:');
%disp(finalData(1:5, :));  % Display the first 5 rows for checking
% Save structure as .mat file
newFilename = [folderPath, '_MaxFieldNormAve.mat'];
save(fullfile(folderPath, newFilename), 'finalData');
disp('Normalization complete. Normalized array saved.');