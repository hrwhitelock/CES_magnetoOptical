
function dataStruct = createData()
%lets start with loading raw data, we should be able to put everything into
% a multi dimensional array, with mag field as an index

% the goal here is to do bg treatment for each field separately, before
% adding together - norm shouldn't really matter
% todo: sort files by field first!!


%lets make this a function 

openFolder = uigetdir; %folder of .mat data files

% saveFolder = 'data/RAW'; %folder to save .mat files

files = dir(fullfile(openFolder, '*.mat'));
firstFile = load(fullfile(files(1).folder, files(1).name));
dataStruct.wavenumbers = firstFile.data(:, 1); %
dataArr = zeros(length(files), length(dataStruct.wavenumbers));


i = 1;
for file = 1:length(files)
    fname = fullfile(openFolder, files(file).name);
    mat = load(fullfile(openFolder, files(file).name));
    % get field to load into separate field arr
    tIndex = find(fname == 'T', 1, 'last'); % finds where T is
    underscoreIndex = find(fname(1:tIndex) == '_', 1, 'last'); % Last occurrence of '_' before 'T'
    field = fname(underscoreIndex+1:tIndex-1); % Including 'T' and two digits after it
    dataStruct.rawField(i) = str2double(field);
    dataArr(i, :) = mat.data(:,2); 
    i = i +1;
end

%now we sort by field

[dataStruct.rawField, I] = sort(dataStruct.rawField);

dataStruct.rawData = dataArr(I, :);


% quickly norm each spec to 1
for i = 1:length(dataStruct.rawData(:,1))
    maxVal = max(dataStruct.rawData(i,:)); 
    for j = 1:length(dataStruct.rawData(1,:))
        dataStruct.rawData(i,j) = dataStruct.rawData(i,j)/maxVal; 
    end
end

