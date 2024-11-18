% find peaks by single spec

function dataStruct = findPeaks(dataStruct)

%assumes you have simplified data to make this a little easier

data = (2-dataStruct.avgData0field(:,:))*1000; 
field = dataStruct.avgField; 
 
dataStruct.pks = NaN(size(data())); 
break1 = find(dataStruct.wavenumbers>10);  
break2 = find(dataStruct.wavenumbers<120); 

pksWave = 0; 
pksField = 0; 

for i = 1:length(field)
    [pks, locs] = findpeaks(data(i,break1(1):break2(end)), 'MinPeakProminence',50, 'MinPeakDistance', .5, 'NPeaks', 15, 'MinPeakHeight', 150);
    findpeaks(data(i,break1(1):break2(end)))
    for j = 1:length(locs)
        dataStruct.pks(i,locs(j)+break1(1)) = pks(j);
        pks(j) = field(i); 
        pksWave = [pksWave, dataStruct.wavenumbers(locs(j))]; 
        pksField = [pksField, pks(j)];
    end
end

pks2d= dataStruct.pks; 

for i = 1:length(pks2d(:,1))
    for j = 1:length(pks2d(1,:))
        if not(isnan(dataStruct.pks(i,j)))
            pks2d(i,j) = field(i); 
        end
    end
end

dataStruct.pks2d = pks2d; 

dataStruct.pksWave = pksWave; 
dataStruct.pksField = pksField; 




