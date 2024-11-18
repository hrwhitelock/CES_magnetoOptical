% plot scratch

%do this so you don't eat up all your ram lol
% 
clear; 
clc; 
fname = uigetfile();


dataStruct = load(fname);
dataStruct = dataStruct.dataStruct; 
%okey, so here I'm going to load all the files







%% make surface plot
figure(); 
hold on; 

[X,Y] = meshgrid(dataStruct.wavenumbers, dataStruct.avgField);
colormap(jet);
s = surf(X,Y, (dataStruct.avgData0field), 'Linestyle', 'none');
s.EdgeColor='none';
% clim([1.2 1.8]);
% zlim([-1 1]);
% zscale('log')
% xlim([0 120])
title('surface plot KErSe2 FIR 2 with peaks')
xlabel('wavenumber (units???)')
ylabel('field')
zlabel('intensity (abs or trans?)')
% 
% [X,Y] = meshgrid(dataStruct.wavenumbers, dataStruct.avgField);
% scatter3(X, Y, dataStruct.pks, 'red')

%% plot single curves
figure(); 
hold on; 
for i = 1:length(dataStruct.avgField)
    % if mod(i,2) ==0
        plot(dataStruct.wavenumbers, i*.5+ dataStruct.avgData0field(i, :))
        scatter(dataStruct.wavenumbers, i*.5+dataStruct.pks(i,:))
    % end
end


%% scatter plot for peaks in 2d
figure(); 
hold on; 
for i = 1:length(dataStruct.avgField)
    % if mod(i,2) ==0
        % plot(dataStruct.wavenumbers, i*.5+ dataStruct.pks)
        scatter(dataStruct.wavenumbers, dataStruct.pks(i,:))
    % end
end

%% make surface plot shit data
figure(); 
hold on; 
[X,Y] = meshgrid(dataStruct.wavenumbers, dataStruct.avgField);
colormap(jet)
% cmp = colormap; 
% cmp = flipud(cmp); 
% colormap(cmp)
s = surf(X,Y, dataStruct.simpleData, 'Linestyle', 'none');
s.EdgeColor='none';
clim([0.3 .8]);
% zlim([-1 1]);
% zscale('log')
title('surface plot KErSe2 FIR 2')
xlabel('wavenumber (units???)')
ylabel('field')
zlabel('intensity (abs or trans?)')

%scatter3(dataStruct.wavenumbers(locs_x), dataStruct.avgField(locs_y), pks)

%% make raw surface plot
figure(); 

[X,Y] = meshgrid(dataStruct.wavenumbers, dataStruct.rawField);
colormap(winter);
s = surf(X,Y, dataStruct.rawData, 'Linestyle', 'none');
s.EdgeColor='none';
% clim([-0.5 0.5]);
% zlim([-1 1]);
title('surface plot KErSe2 FIR 2')
xlabel('wavenumber (units???)')
ylabel('field')
zlabel('intensity (abs or trans?)')


%% make waterfall
figure(); 
[X,Y] = meshgrid(dataStruct.wavenumbers, dataStruct.rawField);

waterfall(X,Y, dataStruct.avgData0field)
%clim([-1 1]);
% zlim([0 2]);
title('Waterfall plot KErSe2 FIR 2')
xlabel('wavenumber (units???)')
ylabel('field')
zlabel('intensity (abs or trans?)')

%% make contour

figure(); 
hold on; 
[X,Y] = meshgrid(dataStruct.wavenumbers, dataStruct.avgField);
c = contourf(X,Y, (2-dataStruct.avgData0field), 150, 'Linestyle', 'none');
view(2)
% clim([1250, 1400])
xlim([4 120])
colormap(cool);
title('CsErSe2 zf subtraction')
xlabel('wavenumber (units???)')
ylabel('field')
scatter(dataStruct.wavenumbers, dataStruct.pks2d, 'black')

%% make contour raw data

figure(); 
[X,Y] = meshgrid(dataStruct.wavenumbers, dataStruct.rawField);
c = contourf(X,Y, dataStruct.rawData, 'Linestyle', 'none');
colormap(winter);
title('KErSe2 w/shenanigans normalization')
xlabel('wavenumber (units???)')
ylabel('field')

