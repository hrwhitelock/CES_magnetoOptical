%This script converts .txt files of wavenumber vs. intensity data into .mat
%files. It will convert any .txt files within the specified folder. 

OpenFolder = 'BG'; %folder of .txt data files

SaveFolder = 'BG_m'; %folder to save .mat files

files = dir(fullfile(OpenFolder, '*.txt'));

for k = 1:length(files)
    % Construct the full file name and load the data
    fileName = fullfile(OpenFolder, files(k).name);
    data = load(fileName);  % Adjust this line if the data needs different preprocessing
    
    % Construct the new file name for the .mat file
    [~, name, ~] = fileparts(files(k).name);
    matFileName = fullfile(SaveFolder, [name '.mat']);
    
    % Save the data to a .mat file
    save(matFileName, 'data');
end

disp('Finished converting files');