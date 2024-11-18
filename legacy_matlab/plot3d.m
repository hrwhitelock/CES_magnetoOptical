function plot3d(dataStruct)

figure(); 

%colors = parula

X = dataStruct.wavenumbers; 
Y = dataStruct.avgField; 
Z = dataStruct.avgData0field; 

for i = 1:length(Y)
    hold on; 
    y = Y(i)*ones(length(X)); 
    plot3(X, y, Z(i, :))
    disp(i)
end

