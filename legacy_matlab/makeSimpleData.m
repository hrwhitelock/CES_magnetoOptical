% okay, so now lets figure out how to recreate the data with only the
% values I care about

%assume we have a dataStruct loaded

data = dataStruct.avgData0field; 
simpleData = zeros(size(data)); 
% meanVal = mean(data, 'all'); 
for i=1:length(data(:,1))
    meanValHorizontal = mean(data(i,:)); 
    for j =1:length(data(1,:))
        meanVal = mean(data(:, j)); 
        if abs(data(i,j)-meanVal) > .03 && data(i,j)<meanVal
            simpleData(i,j) =1-data(i,j);
        else
            simpleData(i,j) = 1-meanVal; 
        end
    end

    if abs(data(i,j)-meanVal) > .05
        simpleData(i,j) =1-data(i,j);
    else
        simpleData(i,j) = 1-meanVal;  
    end
end

dataStruct.simpleData = simpleData; 
