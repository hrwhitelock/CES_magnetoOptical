
%turn this into a function
function dataStruct = normalize(dataStruct)

rawData = dataStruct.rawData;
%%

%%
% choose which spec to use as normalization
% I choose 0 field, for no reason other than I can
% first step is to take all the ZF data and average together
tempData =rawData(1,:); 
counter =1; 
field = dataStruct.rawField;
for i = 2:(length(dataStruct.rawField))
    if field(i)<0.01
        counter = counter +1; 
        for j= 1:length(tempData)
            tempData(j) = tempData(j) + rawData(i,j);
        end
    end
end
for j= 1:length(tempData)
    tempData(j) = tempData(j)/counter;
end
normSpec = tempData; % this chooses the zero field normalization
normArr = zeros(size(dataStruct.rawData)); % initializes dummy array
%dummy array will be left with zeros for the ZF values - maybe that's good,
%maybe its not.. 
% disp(counter)
%actually norm the data - for now subtract, I think? 
for i = counter+1:length(normArr(:,1))
    for j =1:length(normSpec)
        temp = dataStruct.rawData(i, j)-normSpec(j); 
        normArr(i,j) = temp;
    end
end


% 
% %so, this step is experimental, turn off n on as you see fit
% % I'm starting by subtracting off the average spectrum - I'm choosing
% % subtraction instead of division because I don't want small number
% % division issues
% % 
avgArr = normArr(1,:); 
tempField = dataStruct.rawField;
counter2 = 1; %lazy this is redundant
%okay so I'm going to create an average spectrum
for i = counter+1:length(normArr(:,1))
    for j =  1:length(avgArr)
        avgArr(j) = avgArr(j) + normArr(i,j); 
    end
    counter2 = counter2+1; 
end
for  i = 1:length(avgArr)
    avgArr(i) = avgArr(i)/counter2; 
end
%now subtract it off the rawData
%I'm using the raw data array here so that this can be easily turned on and
%off

for i = 1:length(normArr(:,1))-counter
    for j =  1:length(avgArr)
        normArr(i,j) = normArr(i,j)-avgArr(j);
    end
end

%now lets normalize to 1 for the vibes
maxVal = max(normArr, [], 'all'); 
disp(maxVal) % this might be zero and might be the problem
minVal = min(normArr, [], 'all'); 
disp(minVal)
for i = 1:length(normArr(:,1))-counter
    for j =1:length(normSpec)
        normArr(i,j) = normArr(i,j)/maxVal;
    end
end

%shift from 0 for fun
minVal = min(normArr, [], 'all'); 
for i = 1:length(normArr(:,1))
    for j =1:length(normSpec)
        normArr(i,j) = normArr(i,j) + abs(minVal);
    end
end

dataStruct.Norm0field = normArr(counter+1:end, :);
dataStruct.field = dataStruct.rawField(counter+1:end); %removes ZF for yee-haw






% okay, now let's use emily's method of picking the minimum value for each
% column and dividing by that - i don't think this is a good idea
% 
% normArr = zeros(length(dataStruct.rawData(:,1)), length(dataStruct.rawData(1,:)));
% for i = 1:length(dataStruct.field)
%     for j =1:length(normSpec)
%         minVal = min(dataStruct.rawData(:,j));
%         normArr(i,j) = dataStruct.rawData(i, j)/minVal;
%     end
% end
% 
% dataStruct.NormMin = normArr;