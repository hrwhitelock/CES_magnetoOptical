% This script plots raw .mat spectra. 

folderPath = 'P1_EGO_FIR_m'; %Specify the folder containing .mat data files

folderTitle = strrep(folderPath, '_', '\_'); %The folder name is used in plot titles

enableXYPlot = true; %xy plot of all raw spectra on top of each other
enable3Dplot = true; %3D plot showing cascade of spectra. Rotate plot to observe
%=========================================================================
% Begin by sorting the files by their field values going from 0 to max
% field

% Get a list of all .mat files in the folder. 
unsortedfiles = dir(fullfile(folderPath, '*.mat'));

% Initialize an array to hold the numeric values for sorting
numericValues = zeros(length(unsortedfiles), 1);

% Loop through each file to extract the numeric value from the filename
for k = 1:length(unsortedfiles)
    filename = unsortedfiles(k).name;
    % Extract the part of the filename between the last underscore and 'T'
    underscoreIndex = find(filename == '_', 1, 'last');  % Last occurrence of '_'
    tIndex = find(filename == 'T', 1, 'last');           % Last occurrence of 'T'
        % Ensure that both 'T' and '_' are found
    if isempty(tIndex) || isempty(underscoreIndex)
        warning(['Skipping file due to unexpected filename format: ', filename]);
        continue;
    end
    numericStr = filename(underscoreIndex+1:tIndex-1);
    numericValues(k) = str2double(numericStr);          % Convert to double for sorting
end

% Remove zeros from numericValues and corresponding files
validIdx = numericValues ~= 0;
numericValues = numericValues(validIdx);
files = unsortedfiles(validIdx);

% Sort the numeric values, and use the index to sort the files array
[~, sortedIdx] = sort(numericValues);
files = unsortedfiles(sortedIdx);

% Display sorted filenames for verification
disp('Sorted file names:');
for k = 1:length(files)
    disp(files(k).name);
end
disp(['Number of files: ', num2str(length(files))]);


%==========================================================================
% Create a list of labels from the sorted file names, i.e. the field values.
% Labels are used for the plot legends 

labels = cell(length(files), 1);
for k = 1:length(files)
    filename = files(k).name;
    underscoreIndex = find(filename == '_', 1, 'last');
    tIndex = find(filename == 'T', 1, 'first');  % First occurrence of 'T'
    if ~isempty(tIndex)
        labels{k} = filename(underscoreIndex+1:tIndex);
    end
end

% Display the labels for verification
disp('Labels:');
disp(labels);

%==========================================================================
% Create an array containing all the sorted data. The first column
% contains the wavenumber values.

% Initialize the combined array with the wavenumbers from the first file
firstFile = load(fullfile(files(1).folder, files(1).name));
firstFileData = firstFile.data;  % Assuming the data variable is called 'data'
wavenumbers = firstFileData(:, 1); % Assuming the wavenumber is the first column
combined_array = wavenumbers;

% Loop through each sorted file to append the second column to the combined array
for k = 1:length(files)
    dataFile = load(fullfile(files(k).folder, files(k).name));
    data = dataFile.data;  % Assuming the data variable is called 'data'
    combined_array = [combined_array, data(:, 2)]; % Append the second column
end

wavenumbers = combined_array(:, 1);

% Display the combined array size
disp('Combined array size:');
disp(size(combined_array));

%%
if enableXYPlot
    
    num_columns = size(combined_array, 2);
    colors = parula(num_columns); % You can choose 'parula', 'jet', 'viridis', etc.
    
    % Create a new figure
    figure;
    hold on;
    
    % Plot each column from column 2 onward
    for k = 2:size(combined_array, 2)
        y = combined_array(:, k);
        plot(wavenumbers, y, 'DisplayName', labels{k-1}, 'Color', colors(k, :));  % k-1 because labels correspond to columns 2 onward
    end
    hold off;

    % Add labels and title
    xlabel('Wavenumber (cm^{-1})', 'Interpreter', 'tex');
    ylabel('Intensity (a.u.)');
    title(['Raw Data', folderTitle]);
    legend show;  % Display the legend
    grid on;      % Add grid
    box on;       % Add box
end


%%
if enable3Dplot;
    intensity = combined_array(:, 2:end);
    % Extract numerical values for the Z-axis from the labels
    z_values = zeros(1, length(labels));
    for k = 1:length(labels)
        % Extract the part of the label that contains the numeric value
        z_str = regexp(labels{k}, '\d*\.?\d+', 'match');
        z_values(k) = str2double(z_str{1});
    end

    % Create a colormap
    colors = parula(size(intensity, 2));

    % Create a new figure
    figure;
    hold on;

    % Plot each column ofintensity data with corresponding Z position
    for k = 1:size(intensity, 2)
        z = z_values(k) * ones(size(wavenumbers)); % Z values are constant for each column
        y = wavenumbers; % Y-axis is the wavenumbers
        x = intensity(:, k); % X-axis is the intensity data
        plot3(x, y, z, 'Color', colors(k, :), 'DisplayName', labels{k});
    end

    % Add labels and title
    xlabel('Intensity (a.u.)');
    ylabel('Wavenumber (cm^{-1})' ,'Interpreter', 'tex');
    zlabel('$\mu_0$H (T)', 'Interpreter', 'latex');
    title(['3D Raw Data ', folderTitle]);
    legend('show');
    grid on;
    box on;

    hold off;
end
