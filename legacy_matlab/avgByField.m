% make this a function

function dataStruct = avgByField(dataStruct)

% now lets average field vals

%ass data loaded as struct
%structs are goated ngl

counter = 1; 
field = dataStruct.field;
avgData = dataStruct.Norm0field(2,:); %start at 2, our first pass is always sacrificial
tempData = avgData; % redundant, lazy

%first we do the 0 field normalization

dataStruct.avgData0field = [];
dataStruct.avgField = [];
for i = 2:(length(dataStruct.field))
    if field(i)==field(i-1)
        counter = counter +1;
        for j= 1:length(avgData)
            tempData(j) = tempData(j) + dataStruct.Norm0field(i,j);
        end
    else 
        for j= 1:length(avgData)
            tempData(j) = tempData(j)/counter;
        end
        dataStruct.avgData0field = [dataStruct.avgData0field;tempData];
        dataStruct.avgField = [dataStruct.avgField,dataStruct.field(i-1)];
        tempData = dataStruct.Norm0field(i,:);
        counter = 1;
    end
end

% now we do emily's method

% counter = 1; 
% field = dataStruct.field;
% avgData = dataStruct.NormMin(1,:);
% tempData = avgData; % redundant, lazy
% 
% dataStruct.avgDataMin = [];
% dataStruct.avgField = [];
% for i = 2:(length(dataStruct.field))
%     if field(i)==field(i-1)
%         counter = counter +1;
%         for j= 1:length(avgData)
%             tempData(j) = tempData(j) + dataStruct.NormMin(i,j);
%         end
%     else 
%         for j= 1:length(avgData)
%             tempData(j) = tempData(j)/counter;
%         end
%         dataStruct.avgDataMin = [dataStruct.avgDataMin;tempData];
%         dataStruct.avgField = [dataStruct.avgField,dataStruct.field(i-1)];
%         tempData = dataStruct.NormMin(i,:);
%         counter = 1;
%     end
% end