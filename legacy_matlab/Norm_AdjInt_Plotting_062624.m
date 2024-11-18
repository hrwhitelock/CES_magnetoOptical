a% This script works with raw .mat spectra. It applies a +y offset to get rid of small negative 
% values in the data and then applies a max field normalization. 
% It gives options for plotting and further data manipulation

folderPath = 'P1_KeRSe2_FIR_2trial_m'; %Specify the folder containing .mat data files
averagedFolderPath = fullfile(folderPath, 'Averaged_Data'); %Specify folder of averaged .mat data

folderTitle = strrep(folderPath, '_', '\_'); %The folder name is used in plot titles

enableMaxFieldNormAve = false; %xy plot of the MaxField Normalization of the Averaged data
enableStdScale = false; %xy plot of MFNA scaled by standard deviation 
enableZeroFieldNorm = false; %xy plot of MFNA data divided by the averaged zero field data
enableUnchangingFeature = false; %xy plot of data scaled according to variance at each cm-1 from mean of all values at that cm-1
enable3Dplot = true; %Plots MFNA in 3D with field value along z
enable3DAdjusted = false; %Plots MFNA with UnchangingFeature scaling in 3D
enable3DwithEigenSpectrum = false; %Plots MFNA and MFNA w UnchangingFeature along with the full spectrum .fig (must specify .fig file name)
enableSurfPlot = false; % Plots two surface plots of wavenumber vs. field with intensity as color scale, one is MFNA and other is MFNA with UnchangingFeature scaling
enableSurfPlotwEigenSpectrum = false; %Plots MFNA and UnchangingFeature surface plots w full spectrum overlayed (must specify .fig file name)
%=========================================================================
% Begin by sorting the files by their field values going from 0 to max
% field

% Get a list of all .mat files in the folder. 
unsortedfiles = dir(fullfile(averagedFolderPath, '*.mat'));

% Initialize an array to hold the numeric values for sorting
numericValues = zeros(length(unsortedfiles), 1);

% Loop through each file to extract the numeric value from the filename
for k = 1:length(unsortedfiles)
    filename = unsortedfiles(k).name;
    % Extract the part of the filename between the last underscore and 'T'
    tIndex = find(filename == 'T', 1, 'last'); % Index of last occurrence of 'T'
    numericStr = filename(1:tIndex-1); %Extract the value
    numericValues(k) = str2double(numericStr); % Convert to double for sorting
end


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
    tIndex = find(filename == 'T', 1, 'first');  % First occurrence of 'T'
    if ~isempty(tIndex)
        labels{k} = filename(1:tIndex);
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
firstFileData = firstFile.accumulatedData;  % Assuming the data variable is called 'data'
wavenumbers = firstFileData(:, 1); % Assuming the wavenumber is the first column
combined_array = wavenumbers;

% Loop through each sorted file to append the second column to the combined array
for k = 1:length(files)
    dataFile = load(fullfile(files(k).folder, files(k).name));
    data = dataFile.accumulatedData;  % Assuming the data variable is called 'data'
    combined_array = [combined_array, data(:, 2)]; % Append the second column
end

% Display the combined array size
disp('Combined array size:');
disp(size(combined_array));

% combined_array has the wavenumbers in the first column and intensity values in subsequent columns

%==========================================================================
% Add a +y offset to all the data to get rid of any negative intensity
% values. This should be a very small value << 1.  

% Find the minimum value in columns 2 onward
minValue = min(min(combined_array(:, 2:end)));

% Display the minimum value. Verify that the value is small. The offset
% should be negligible to the data overall
disp(['Minimum intensity value ', num2str(minValue)]);

% Shift all values in columns 2 onward by the shift value plus a little
% extra so there are no zero values that will mess up normalization
combined_array(:, 2:end) = combined_array(:, 2:end) + abs(minValue) + 0.0000001;

%==========================================================================
% Apply a max field normalization and invert the result

% Max field normalization takes the max value for each wavenumber across
% all fields and divides all data at each wavenumber by the corresponding 
% max value

% Find the maximum value in each row from columns 2 onward (ignore
% wavenumbers in columns 1)
maxValues = max(combined_array(:, 2:end), [], 2);

% Normalize the values in combined_array from columns 2 onward by dividing by the maximum values
normalized_combined_array = combined_array; % Initialize the normalized array   
normalized_combined_array(:, 2:end) = combined_array(:, 2:end) ./ maxValues;
normalized_combined_array(:, 2:end) = 1 - normalized_combined_array(:, 2:end);

%==========================================================================
% Plot the Averaged Max Field Normalized data
if enableMaxFieldNormAve
    % Extract the x-axis data (wavenumbers)
    x = normalized_combined_array(:, 1);
    
    num_columns = size(normalized_combined_array, 2);
    colors = parula(num_columns); % You can choose 'parula', 'jet', 'viridis', etc.
    
    % Create a new figure
    figure;
    hold on;
    
    % Plot each column from column 2 onward
    for k = 2:size(normalized_combined_array, 2)
        y = normalized_combined_array(:, k);
        plot(x, y, 'DisplayName', labels{k-1}, 'Color', colors(k, :));  % k-1 because labels correspond to columns 2 onward
    end
    hold off;

    % Add labels and title
    xlabel('Wavenumber (cm^{-1})', 'Interpreter', 'tex');
    ylabel('Intensity (a.u.)');
    title(['MaxFieldNorm of Averaged Data w +y offset ', folderTitle]);
    legend show;  % Display the legend
    grid on;      % Add grid
    box on;       % Add box
end

%==========================================================================
%Plot the standard deviation of the intensity along with the normalized data 
% This attempts to add more clarity to the shifting peaks but it doesn't
% work very well

if enableStdScale
    % Extract the wavenumbers (first column)
    wavenumbers = normalized_combined_array(:, 1);

    % Extract the normalized intensity values (columns 2 onward)
    normalized_intensity = normalized_combined_array(:, 2:end);

    % Calculate the standard deviation of intensity at each wavenumber across all spectra
    intensity_std = std(normalized_intensity, 0, 2);
    std2 = intensity_std;
    
    % This scales the data by the standard deviation 
    % Adjust the scaling for best results
    scale_std = intensity_std * 100;
    normalized_intensity = normalized_intensity .* scale_std;
    std2 = intensity_std * 50;
         
    
    num_columns = size(normalized_intensity, 2);
    colors = parula(num_columns); % You can choose 'parula', 'jet', 'viridis', etc.

    % Create a new figure
    figure;
    hold on;

    % Plot each column of scaled data with corresponding labels
    num_columns = size(normalized_intensity, 2);
    for k = 1:num_columns
        plot(wavenumbers, normalized_intensity(:, k), 'DisplayName', labels{k}, 'Color', colors(k, :));
    end

    % Plot the standard deviation
    plot(wavenumbers, std2, 'k', 'LineWidth', 2, 'DisplayName', 'Standard Deviation');

    % Add labels and title
    xlabel('Wavenumber (cm^{-1})', 'Interpreter', 'tex');
    ylabel('Scaled Intensity (a.u.)');
    title(['MaxField Normalization with Standard Deviation Scaling ', folderTitle]);
    legend('show');
    grid on;
    box on;

    hold off;
end

%=========================================================================
% This takes the normalized data and subtracts the zero field (column 2)
% data, attempting to add more clarity to data. 

if enableZeroFieldNorm
    % Extract the x-axis data (wavenumbers)
    x = normalized_combined_array(:, 1);
    
    % Subtract column 2 from columns 3 onward
    subtracted_data = normalized_combined_array(:, 3:end) ./ normalized_combined_array(:, 2);

    % Define the colormap
    num_columns = size(subtracted_data, 2);
    colors = parula(num_columns); % You can choose 'parula', 'jet', 'viridis', etc.

    % Create a new figure
    figure;

    % Plot each column of subtracted data with corresponding labels
    hold on;
    for k = 1:num_columns
        y = subtracted_data(:, k);
        plot(x, y, 'DisplayName', labels{k+1}, 'Color', colors(k, :));  % k+1 because labels correspond to columns 3 onward
    end
    hold off;

    % Add labels and title
    xlabel('Wavenumber (cm^{-1})', 'Interpreter', 'tex');
    ylabel('Intensity (a.u.)');
    title(['MaxFieldNorm of Averaged Data with Y Offset ', folderTitle]);
    legend show;  % Display the legend
    grid on;      % Add grid
    box on;       % Add box
end

%==========================================================================
% Scales data according to variance at each cm-1 from mean of all values at that cm-1

if enableUnchangingFeature
    % Extract the wavenumbers (first column)
    wavenumbers = normalized_combined_array(:, 1);

    % Extract the normalized intensity values (columns 2 onward)
    normalized_intensity = normalized_combined_array(:, 2:end);

    % Initialize the adjusted intensity array
    adjusted_intensity = normalized_intensity;

    % Loop through each row
    for i = 1:size(normalized_intensity, 1)
        % Calculate the mean of the current row
        row_mean = mean(normalized_intensity(i, :));

        % Calculate the percent variance for each value in the row
        percent_variance = abs(normalized_intensity(i, :) - row_mean) / row_mean;

        % Scale the data according to its variance from the mean. Low
        % varying points are set to zero
        for j = 1:size(normalized_intensity, 2)
            if percent_variance(j) <= 0.5 %Adjust value for best results
                adjusted_intensity(i, j) = 0;
            else 
                adjusted_intensity(i,j) = normalized_intensity(i, j) * percent_variance(j); %can adjust scaling
            end
        end
    end

    % Combine adjusted data with the first column (wavenumbers)
    adjusted_combined_array = [wavenumbers, adjusted_intensity];

    % Define the colormap
    num_columns = size(adjusted_combined_array, 2) - 1;  % Exclude the wavenumber column
    colors = parula(num_columns);  % You can choose 'parula', 'jet', 'viridis', etc.

    % Create a new figure
    figure;
    hold on;

    % Plot each column of adjusted intensity data with corresponding labels
    for k = 1:num_columns
        y = adjusted_combined_array(:, k + 1);  % Offset by 1 to account for wavenumbers
        plot(wavenumbers, adjusted_intensity(:, k), 'DisplayName', labels{k}, 'Color', colors(k, :));
    end

    % Add labels and title
    xlabel('Wavenumber (cm^{-1})', 'Interpreter', 'tex');
    ylabel('Intensity (a.u.)');
    title(['Variation from Mean Adjusted Intensity Data ', folderTitle]);
    legend('show');
    grid on;
    box on;

    hold off;
end

%==========================================================================
%3D plot of MFNA data with no additional scaling applied

if enable3Dplot;

    % Extract the wavenumbers (first column) as Y
    wavenumbers = normalized_combined_array(:, 1);

    % Extract the normalized intensity values (columns 2 onward)
    normalized_intensity = normalized_combined_array(:, 2:end);

    % Extract numerical values for the Z-axis from the labels
    z_values = zeros(1, length(labels));
    for k = 1:length(labels)
        % Extract the part of the label that contains the numeric value
        z_str = regexp(labels{k}, '\d*\.?\d+', 'match');
        z_values(k) = str2double(z_str{1});
    end

    % Create a colormap
    colors = parula(size(normalized_intensity, 2));

    % Create a new figure
    figure;
    hold on;

    % Plot each column of normalized intensity data with corresponding Z position
    for k = 1:size(normalized_intensity, 2)
        z = z_values(k) * ones(size(wavenumbers)); % Z values are constant for each column
        y = wavenumbers; % Y-axis is the wavenumbers
        x = normalized_intensity(:, k); % X-axis is the normalized intensity data
        plot3(x, y, z, 'Color', colors(k, :), 'DisplayName', labels{k});
    end

    % Add labels and title
    xlabel('Intensity (a.u.)');
    ylabel('Wavenumber (cm^{-1})' ,'Interpreter', 'tex');
    zlabel('$\mu_0$H (T)', 'Interpreter', 'latex');
    title(['MaxFieldNorm ', folderTitle]);
    legend('show');
    grid on;
    box on;

    hold off;
end

%==========================================================================
% 3D plot of data scaled according to UnchangingFeature percent variance 
if enable3DAdjusted
    % Extract the wavenumbers (first column)
    wavenumbers = normalized_combined_array(:, 1);

    % Extract the normalized intensity values (columns 2 onward)
    normalized_intensity = normalized_combined_array(:, 2:end);

    % Initialize the adjusted intensity array
    adjusted_intensity = normalized_intensity;

    % Loop through each row to adjust the intensity values
    for i = 1:size(normalized_intensity, 1)
        % Calculate the mean of the current row
        row_mean = mean(normalized_intensity(i, :));

        % Calculate the percent variance for each value in the row
        percent_variance = abs(normalized_intensity(i, :) - row_mean) / row_mean;

        % Scale the data according to its variance from the mean. Low
        % varying points are set to zero
        for j = 1:size(normalized_intensity, 2)
            if percent_variance(j) < 0.1 %Adjust value for best results
                adjusted_intensity(i, j) = 0;
            else
                adjusted_intensity(i, j) = normalized_intensity(i, j) * percent_variance(j); % Can adjust scaling
            end
        end
    end

    % Extract numerical values for the Z-axis from the labels
    z_values = zeros(1, length(labels));
    for k = 1:length(labels)
        % Extract the part of the label that contains the numeric value
        z_str = regexp(labels{k}, '\d*\.?\d+', 'match');
        z_values(k) = str2double(z_str{1});
    end

    % Create a colormap
    colors = parula(size(adjusted_intensity, 2));

    % Create a new figure
    figure;
    hold on;

    % Plot each column of adjusted intensity data with corresponding Z position
    for k = 1:size(adjusted_intensity, 2)
        z = z_values(k) * ones(size(wavenumbers)); % Z values are constant for each column
        y = wavenumbers; % Y-axis is the wavenumbers
        x = adjusted_intensity(:, k); % X-axis is the adjusted intensity data
        plot3(x, y, z, 'Color', colors(k, :), 'DisplayName', labels{k});
    end

    % Add labels and title
    xlabel('Adjusted Intensity (a.u.)');
    ylabel('Wavenumber (cm^{-1})', 'Interpreter', 'tex');
    zlabel('$\mu_0$H (T)', 'Interpreter', 'latex');
    title(['MaxFieldNormAve Variation from Mean Adjusted ', folderTitle]);
    legend('show');
    grid on;
    box on;

    hold off;
end

%=========================================================================
if enable3DwithEigenSpectrum
    % Must specify the path to your .fig file (this expects a figure with
    % wavenumber on y axis and H[T] on x axis)
    
    figFilePath = 'Normalized_Data/P1_EGO_FIR/EGO_FullSpectrum_b.fig';
    
    % Open the figure
    fig = openfig(figFilePath, 'invisible'); % Load the figure without displaying it

    % Find all line objects in the figure
    lines = findobj(fig, 'Type', 'line');

    % Initialize cell arrays to store the data
    xData = cell(length(lines), 1);
    yData = cell(length(lines), 1);

    % Loop through each line object and extract the X and Y data
    for k = 1:length(lines)
        xData{k} = get(lines(k), 'XData');
        yData{k} = get(lines(k), 'YData');
    end

    % Close the figure (optional)
    close(fig);

    % Remove the first dataset
    xData(1) = [];
    yData(1) = [];

    % Crop the data. Adjust bounds as necessary
    max_x = 20;
    max_y = 10000;
    cropped_xData = cellfun(@(x) x(x <= max_x), xData, 'UniformOutput', false);
    cropped_yData = cellfun(@(y) y(y <= max_y), yData, 'UniformOutput', false);

    % Combine X and Y data into two-column datasets
    datasets = cell(length(cropped_xData), 1);
    for k = 1:length(cropped_xData)
        % Ensure both x and y data have the same length after cropping
        len = min(length(cropped_xData{k}), length(cropped_yData{k}));
        datasets{k} = [cropped_xData{k}(1:len)', cropped_yData{k}(1:len)']; % Transpose to make sure X and Y are columns
    end

    % Extract the wavenumbers (first column) as Y
    wavenumbers = normalized_combined_array(:, 1);

    % Extract the normalized intensity values (columns 2 onward)
    normalized_intensity = normalized_combined_array(:, 2:end);

    % Extract numerical values for the Z-axis from the labels
    z_values = zeros(1, length(labels));
    for k = 1:length(labels)
        % Extract the part of the label that contains the numeric value
        z_str = regexp(labels{k}, '\d*\.?\d+', 'match');
        z_values(k) = str2double(z_str{1});
    end

    % Create a colormap
    colors = parula(size(normalized_intensity, 2));

    % Create a new figure
    figure;
    hold on;

    % Plot each column of normalized intensity data with corresponding Z position
    for k = 1:size(normalized_intensity, 2)
        z = z_values(k) * ones(size(wavenumbers)); % Z values are constant for each column
        y = wavenumbers; % Y-axis is the wavenumbers
        x = normalized_intensity(:, k); % X-axis is the normalized intensity data
        plot3(x, y, z, 'Color', colors(k, :), 'DisplayName', labels{k});
    end

    % Plot the 2D data on the y-z plane
    for k = 1:length(datasets)
        y = datasets{k}(:, 2); % Y-axis is the cropped yData
        z = datasets{k}(:, 1); % Z-axis is the cropped xData
        x = zeros(size(y)); % X-axis is zero to place it on the y-z plane
        plot3(x, y, z, 'k', 'LineWidth', 2, 'DisplayName', ['Eigen Data ' num2str(k)]);
    end

    % Add labels and title
    xlabel('Intensity (a.u.)');
    ylabel('Wavenumber (cm^{-1})', 'Interpreter', 'tex');
    zlabel('$\mu_0$H (T)', 'Interpreter', 'latex');
    title(['MaxFieldNormAve w Eigen Spectrum', folderTitle]);
    legend('show');
    grid on;
    box on;

    hold off;

    % Initialize the adjusted intensity array
    adjusted_intensity = normalized_intensity;

    % Loop through each row to adjust the intensity values
    for i = 1:size(normalized_intensity, 1)
        % Calculate the mean of the current row
        row_mean = mean(normalized_intensity(i, :));

        % Calculate the percent variance for each value in the row
        percent_variance = abs(normalized_intensity(i, :) - row_mean) / row_mean;

        % Scale the data according to its variance from the mean. Low
        % varying points are set to zero
        for j = 1:size(normalized_intensity, 2)
            if percent_variance(j) < 0.1 %Adjust value for best results
                adjusted_intensity(i, j) = 0;
            else
                adjusted_intensity(i, j) = normalized_intensity(i, j) * percent_variance(j); % Can adjust scaling
            end
        end
    end

    % Extract numerical values for the Z-axis from the labels
    z_values = zeros(1, length(labels));
    for k = 1:length(labels)
        % Extract the part of the label that contains the numeric value
        z_str = regexp(labels{k}, '\d*\.?\d+', 'match');
        z_values(k) = str2double(z_str{1});
    end

    % Create a colormap
    colors = parula(size(adjusted_intensity, 2));

    % Create a new figure
    figure;
    hold on;

    % Plot each column of adjusted intensity data with corresponding Z position
    for k = 1:size(adjusted_intensity, 2)
        z = z_values(k) * ones(size(wavenumbers)); % Z values are constant for each column
        y = wavenumbers; % Y-axis is the wavenumbers
        x = adjusted_intensity(:, k); % X-axis is the adjusted intensity data
        plot3(x, y, z, 'Color', colors(k, :), 'DisplayName', labels{k});
    end

    % Plot the 2D data on the y-z plane
    for k = 1:length(datasets)
        y = datasets{k}(:, 2); % Y-axis is the cropped yData
        z = datasets{k}(:, 1); % Z-axis is the cropped xData
        x = zeros(size(y)); % X-axis is zero to place it on the y-z plane
        plot3(x, y, z, 'k', 'LineWidth', 2, 'DisplayName', ['Eigen Data ' num2str(k)]);
    end

    % Add labels and title
    xlabel('Adjusted Intensity (a.u.)');
    ylabel('Wavenumber (cm^{-1})', 'Interpreter', 'tex');
    zlabel('$\mu_0$H (T)', 'Interpreter', 'latex');
    title(['Adjusted Intensity MaxFielNormAve with Eigen Spectrum', folderTitle]);
    legend('show');
    grid on;
    box on;

    hold off;
end
%=========================================================================
% Plots a surface plot of wavenumber vs. field with intensity as color
% scale
if enableSurfPlot

    % Extract the wavenumbers (first column) as Y
    wavenumbers = normalized_combined_array(:, 1);

    % Extract the normalized intensity values (columns 2 onward)
    normalized_intensity = normalized_combined_array(:, 2:end);

    % Initialize the adjusted intensity array
    adjusted_intensity = normalized_intensity;

    % Loop through each row to adjust the intensity values
    for i = 1:size(normalized_intensity, 1)
        % Calculate the mean of the current row
        row_mean = mean(normalized_intensity(i, :));

        % Calculate the percent variance for each value in the row
        percent_variance = abs(normalized_intensity(i, :) - row_mean) / row_mean;

       % Scale the data according to its variance from the mean. Low
        % varying points are set to zero
        for j = 1:size(normalized_intensity, 2)
            if percent_variance(j) < 0.1 %Adjust value for best results
                adjusted_intensity(i, j) = 0;
            else
                adjusted_intensity(i, j) = normalized_intensity(i, j) * percent_variance(j); % Can adjust scaling
            end
        end
    end

    % Extract numerical values for the magnetic field (T) from the labels
    magnetic_field = zeros(1, length(labels));
    for k = 1:length(labels)
        % Extract the part of the label that contains the numeric value
        field_str = regexp(labels{k}, '\d*\.?\d+', 'match');
        magnetic_field(k) = str2double(field_str{1});
    end

    % Create a 2D grid for wavenumbers and magnetic field
    [Field, Wavenumber] = meshgrid(magnetic_field, wavenumbers);

    % Create a new figure for the surface plot
    figure;
    rot_adj_int = flipud(rot90(adjusted_intensity));
    % Plot the surface plot
    imagesc(wavenumbers, magnetic_field, rot_adj_int);
    set(gca, 'YDir', 'normal'); % Correct the Y-axis direction
    colormap parula; % Apply the colormap
    colorbar; % Display the color bar

    % Add labels and title
    ylabel('$\mu_0$H (T)', 'Interpreter', 'latex');
    xlabel('Wavenumber (cm^{-1})', 'Interpreter', 'tex');
    title(['Adjusted Intensity MaxFieldNormAve ', folderTitle]);

    
    
    % Create a new figure for the surface plot
    figure;
    % Rotates data, optional for how plot is presented
    norm_int_rot = flipud(rot90(normalized_intensity));
    % Plot the surface plot
    imagesc(wavenumbers, magnetic_field, norm_int_rot);
    set(gca, 'YDir', 'normal'); % Correct the Y-axis direction
    colormap parula; % Apply the colormap
    colorbar; % Display the color bar

    % Add labels and title
    ylabel('$\mu_0$H (T)', 'Interpreter', 'latex');
    xlabel('Wavenumber (cm^{-1})', 'Interpreter', 'tex');
    title(['MaxFieldNormAve ', folderTitle]);
 
end
%=========================================================================

if enableSurfPlotwEigenSpectrum
    
    %Specify the figure file name for the full spectrum
    figFilePath = 'Normalized_Data/P1_KeRSe2_FIR_2trial/KErSe2_FullSpec_Peak.fig';

    % Open the figure
    fig = openfig(figFilePath, 'invisible'); % Load the figure without displaying it

    % Find all line objects in the figure
    lines = findobj(fig, 'Type', 'line');

    % Initialize cell arrays to store the data
    xData = cell(length(lines), 1);
    yData = cell(length(lines), 1);

    % Loop through each line object and extract the X and Y data
    for k = 1:length(lines)
        xData{k} = get(lines(k), 'XData');
        yData{k} = get(lines(k), 'YData');
    end

    % Close the figure (optional)
    close(fig);

    % Remove the first dataset
    xData(1) = [];
    yData(1) = [];

    % Crop the data. Adjust as necessary
    max_x = 20;
    max_y = 10000;
    cropped_xData = cellfun(@(x) x(x <= max_x), xData, 'UniformOutput', false);
    cropped_yData = cellfun(@(y) y(y <= max_y), yData, 'UniformOutput', false);

    % Combine X and Y data into two-column datasets
    datasets = cell(length(cropped_xData), 1);
    for k = 1:length(cropped_xData)
        % Ensure both x and y data have the same length after cropping
        len = min(length(cropped_xData{k}), length(cropped_yData{k}));
        datasets{k} = [cropped_xData{k}(1:len)', cropped_yData{k}(1:len)']; % Transpose to make sure X and Y are columns
    end

    % Extract the wavenumbers (first column) as Y
    wavenumbers = normalized_combined_array(:, 1);

    % Extract the normalized intensity values (columns 2 onward)
    normalized_intensity = normalized_combined_array(:, 2:end);

    % Initialize the adjusted intensity array
    adjusted_intensity = normalized_intensity;

    % Loop through each row to adjust the intensity values
    for i = 1:size(normalized_intensity, 1)
        % Calculate the mean of the current row
        row_mean = mean(normalized_intensity(i, :));

        % Calculate the percent variance for each value in the row
        percent_variance = abs(normalized_intensity(i, :) - row_mean) / row_mean;

       % Scale the data according to its variance from the mean. Low
        % varying points are set to zero
        for j = 1:size(normalized_intensity, 2)
            if percent_variance(j) < 0.1 %Adjust value for best results
                adjusted_intensity(i, j) = 0;
            else
                adjusted_intensity(i, j) = normalized_intensity(i, j) * percent_variance(j); % Can adjust scaling
            end
        end
    end

    % Extract numerical values for the magnetic field (T) from the labels
    magnetic_field = zeros(1, length(labels));
    for k = 1:length(labels)
        % Extract the part of the label that contains the numeric value
        field_str = regexp(labels{k}, '\d*\.?\d+', 'match');
        magnetic_field(k) = str2double(field_str{1});
    end

    % Create a 2D grid for wavenumbers and magnetic field
    [Field, Wavenumber] = meshgrid(magnetic_field, wavenumbers);

    % Create a new figure for the surface plot
    figure;

    % Plot the surface plot
    imagesc(magnetic_field, wavenumbers, adjusted_intensity);
    set(gca, 'YDir', 'normal'); % Correct the Y-axis direction
    colormap parula; % Apply the colormap
    colorbar; % Display the color bar

    % Add labels and title
    xlabel('$\mu_0$H (T)', 'Interpreter', 'latex');
    ylabel('Wavenumber (cm^{-1})', 'Interpreter', 'tex');
    title(['Adjusted Intensity MaxFieldNormAve with Eigen Spectrum', folderTitle]);

    % Overlay the cropped 2D data with maximum intensity
    hold on;
    max_intensity = max(adjusted_intensity(:));
    for k = 1:length(datasets)
        x = datasets{k}(:, 1); % Magnetic field values
        y = datasets{k}(:, 2); % Wavenumber values
        plot(x, y, 'r--', 'LineWidth', 1.5, 'DisplayName', ['Eigen Data ' num2str(k)]);
    end
    legend('show');
    hold off; 
    
    % Create a new figure for the surface plot
    figure;

    % Plot the surface plot
    imagesc(magnetic_field, wavenumbers, normalized_intensity);
    set(gca, 'YDir', 'normal'); % Correct the Y-axis direction
    colormap parula; % Apply the colormap
    colorbar; % Display the color bar

    % Add labels and title
    xlabel('$\mu_0$H (T)', 'Interpreter', 'latex');
    ylabel('Wavenumber (cm^{-1})', 'Interpreter', 'tex');
    title(['MaxFieldNormAve with Eigen Spectrum', folderTitle]);

    % Overlay the cropped 2D data with maximum intensity
    hold on;
    max_intensity = max(adjusted_intensity(:));
    for k = 1:length(datasets)
        x = datasets{k}(:, 1); % Magnetic field values
        y = datasets{k}(:, 2); % Wavenumber values
        plot(x, y, 'r--', 'LineWidth', 1.5, 'DisplayName', ['Eigen Data ' num2str(k)]);
    end
    legend('show');
    hold off; 
end
 
