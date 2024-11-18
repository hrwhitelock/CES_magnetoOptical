% here's a little function I'm going to use to average the normalized data
% and then subtract. is this good? probably not. will I do it? yes!

function dataStruct = shenanigans(dataStruct)

dataStruct.shitData = dataStruct.avgData0field; 
for k = 1:1
    tempData = dataStruct.shitData(1,:); 
    tempField = dataStruct.avgField;
    counter = 1; %lazy this is redundant
    
    %okay so I'm going to create an average spectrum
    for i = 2:length(tempField)
        for j =  1:length(tempData)
            tempData(j) = tempData(j) + dataStruct.shitData(i,j);
            dataStruct.shitData(i,j) = dataStruct.shitData(i,j)-tempData(j)/2; 
        end
        counter = counter+1; 
        tempData = dataStruct.shitData(i,:);
    end
    for  i = 1:length(tempData)
        tempData(i) = tempData(i)/counter; 
    end
    % now I'm going to take that average spec and subtract it from each spec
    % this is good and makes sense and is not bad
    for i = 1:length(tempField)
        for j =  1:length(tempData)
            dataStruct.shitData(i,j) = dataStruct.shitData(i,j)-tempData(j);
        end
    end
    
    minVal = min(dataStruct.shitData, [], 'all');

end
for i = 1:length(tempField)
    for j =  1:length(tempData)
        dataStruct.shitData(i,j) = dataStruct.shitData(i,j)+abs(minVal);
    end
end
maxVal = max(dataStruct.shitData, [], 'all');

for i = 1:length(tempField)
    for j =  1:length(tempData)
        dataStruct.shitData(i,j) = dataStruct.shitData(i,j)/abs(maxVal);
    end
end