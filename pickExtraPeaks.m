function pickExtraPeaks(dataStruct)

% first we do the plotting
figure(); 
hold on; 

figure(); 
hold on; 
[X,Y] = meshgrid(dataStruct.wavenumbers, dataStruct.avgField);
c = contourf(X,Y, 1-dataStruct.avgData0field,200, 'Linestyle', 'none');
view(2)
clim([-.2,.8])
colormap(winter);
title('KErSe2 w/just ZF subtraction')
xlabel('wavenumber (units???)')
ylabel('field')

scatter(dataStruct.wavenumbers, dataStruct.pks2d, 'red')

