% load all data into mat files
clear; 

path ="/Users/hopeless/Desktop/LeeLab/data/IR_MagLab_Feb2024/Raw_mat/"; 
folders = ["P1_EGO_FIR_m" "P1_EGO_MIR_m" "P1_KeRSe2_FIR_2trial.m" "P1_KErSe2_FIR_m" "P2_ErGaO_FIR_m" "P2_ErGaO_MIR_m" "P2_NES_FIR_m" "P3_ErVO4_FIR_800_10_m" "P3_ErVO4_FIR_m" "P3_ErVO4_MIR_m" "P3_EVO_FIR_m" "P4_Er166_FIR_m" "P4_Er166_MIR_m" "P5_Dy166_FIR_m" "P5_Dy166_MIR_m" "P6_ETO_FIR_m"];
 
for i = 1:length(folders)
    fullpath = strcat(path, folders(i)); 
    disp(fullpath)
    createData(fullpath)
    % makeAndSave(erase(folders(i),'_m '), fullpath)
end