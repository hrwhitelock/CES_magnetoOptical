%subtract off Zerofield
%only doing this for ZF norm
function dataStruct = ZFSub(dataStruct)

normSpec = dataStruct.avgData0field(1, :); % this chooses the zero field normalization
normArr = zeros(size(dataStruct.avgData0field));
for i = 2:length(dataStruct.avgField)
    for j =1:length(normSpec)
        normArr(i,j) = dataStruct.avgData0field(i, j)/normSpec(j);
    end
end

maxVal = max(normArr, [], 'all');

for i = 1:length(normArr(:,1))
    for j=1:length(normSpec)
        normArr(i,j) = normArr(i,j)/maxVal;
    end
end

dataStruct.data0fieldSub = normArr; 
dataStruct.field0removed = dataStruct.avgField(2:end); 
