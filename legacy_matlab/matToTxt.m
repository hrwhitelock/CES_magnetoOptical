% convert back to txt because I don't get matlab & would prefer more basic
% format
openFolder = 'data/IR_MAgLab_Feb2024/Raw_mat/P1_KeRSe2_FIR_m'; %folder of .mat data files

% saveFolder = 'data/RAW'; %folder to save .mat files

files = dir(fullfile(openFolder, '*.mat'));

for file = 1:(length(files))
    fname = fullfile(openFolder, files(file).name);
    mat = load(fullfile(openFolder, files(file).name));
    % Do some stuff
    %Conten = who;
    F_name = append(erase(fname, '.mat'), '.txt');
    writematrix(mat.data, F_name)
end