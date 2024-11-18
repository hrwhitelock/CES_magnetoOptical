%make and save data
% does not work on ETO_FIR, not sure why??


function dataStruct = makeAndSave(fname)

%first call creatData
dataStruct = createData(); 
dataStruct = normalize(dataStruct); 
dataStruct = avgByField(dataStruct);

save(fname, "dataStruct"); 