% final fig generation script 12/16/24

% let's start with making the paper figures
%% load my spectroscopy data
filename = 'spectroscopy_calculation_hopes_params.h5';

info = h5info(filename, '/');
my_spec_data = struct; 
for i = 1:length(info.Datasets)
    % Get the dataset name
    datasetName = info.Datasets(i).Name;
    % Load the dataset into a variable
    datasetData = h5read(filename, ['/', datasetName]);
    % Optionally, store the data in a struct (to make it easy to access)
    my_spec_data.(datasetName) = datasetData;
end

%% load Allens spectroscopy data
filename = 'spectroscopy_calculation_Allens_params.h5';

info = h5info(filename, '/');
allens_spec_data = struct; 
for i = 1:length(info.Datasets)
    % Get the dataset name
    datasetName = info.Datasets(i).Name;
    % Load the dataset into a variable
    datasetData = h5read(filename, ['/', datasetName]);
    % Optionally, store the data in a struct (to make it easy to access)
    allens_spec_data.(datasetName) = datasetData;
end

%% 
fig = figure;
ax1 = subplot(2,1,1);
title('B-axis IR')
hold on; 
contourf(my_spec_data.IR_B_field, my_spec_data.IR_B_wavenums, my_spec_data.IR_dataB', 100, 'LineStyle', 'none')
for i = 2:16%size(arrC, 2)
    plot(my_spec_data.calc_field, my_spec_data.linesB(:,i),  'r:', 'LineWidth', 0.5);
end
clim([0 1])
ylim([0 100])
xlim([0 17.5])
colormap(ax1, flipud(parula))
ax2 = subplot(2,1,2);
hold on; 
contourf(my_spec_data.calc_field, my_spec_data.calc_wavenums, my_spec_data.simulated_IR_B, 100,'LineStyle', 'none', 'EdgeSmoothing', 'none');
colormap(ax4, flipud(parula))
clim([0 1])

for i = 2:16%size(arrC, 2)
    plot(my_spec_data.calc_field, my_spec_data.linesB(:,i),  'r:', 'LineWidth', 0.5);
end

clim([0 1])
ylim([0 100])
xlim([0 17.5])

%%
ax2 = subplot(2,3,2);hold on; box on; 
contourf(my_spec_data.IR_C_field, my_spec_data.IR_C_wavenums, my_spec_data.IR_dataC', 100, 'LineStyle', 'none')
for i = 2:16%size(arrC, 2)
    plot(my_spec_data.calc_field,my_spec_data.linesC(:,i),  'r:', 'LineWidth', 0.5);
end
title('C-axis IR')
colormap(ax2, flipud(parula))
clim([0 1])
ylim([0 100])
xlim([0 17.5])

ax5 = subplot(2,3,5);
hold on; box on; 
contourf(my_spec_data.calc_field, my_spec_data.calc_wavenums, my_spec_data.simulated_IR_C, 100, 'LineStyle', 'none');
colormap(ax5, flipud(parula))
for i = 2:16%size(arrC, 2)
    plot(my_spec_data.calc_field,my_spec_data.linesC(:,i),  'r:', 'LineWidth', 0.5);
end
clim([0 1])
ylim([0 100])
xlim([0 17.5])




% fig = figure;
ax3 = subplot(2,3,3);box on; hold on; 
contourf(my_spec_data.raman_field, my_spec_data.raman_wavenums, my_spec_data.ramanData', 100, 'LineStyle', 'none')
for i = 2:16%size(arrC, 2)
    plot(my_spec_data.calc_field,my_spec_data.linesC(:,i),  'r:', 'LineWidth', 0.5);
end
ylim([0 100])
xlim([0 17.5])
colormap(ax3, 'cool')
title ('C-axis Raman')
ax6 = subplot(2,3,6);
hold on; box on; 
contourf(my_spec_data.calc_field, my_spec_data.calc_wavenums, my_spec_data.simulated_raman, 100, 'LineStyle', 'none');
colormap(ax6, 'cool')
% let's add some lines
for i = 2:16%size(arrC, 2)
    plot(my_spec_data.calc_field,my_spec_data.linesC(:,i),  'r:', 'LineWidth', 0.5);
end
clim([0 1])
xlim([0 17.5])
ylim([0 100])




ax3.Position(2)=ax6.Position(2)+ax6.Position(4)+.02;
ax1.Position(2)=ax4.Position(2)+ax4.Position(4)+.02;
% 
% ax4.Position(2)=ax6.Position(2)+ax6.Position(4);
ax2.Position(2)=ax5.Position(2)+ax5.Position(4)+.02;
% 
% 
% ax4.Position(1)=ax3.Position(1)+ax3.Position(3);
ax2.Position(1)=ax1.Position(1)+ax1.Position(3);
ax5.Position(1)=ax4.Position(1)+ax4.Position(3);
linkaxes([ax1, ax2, ax3, ax4, ax5, ax6], 'xy')
% ylabel('Energy [cm^{-1}]')
% xlabel('H [t]')
set(ax2,'Yticklabel',[]) 
set(ax1,'Xticklabel',[])
set(ax2,'Xticklabel',[])

set(ax5,'Yticklabel',[]) 
% set(ax5,'Xticklabel',[])


set(ax3,'Xticklabel',[])

%% now make fig 3
% load my magnetic data
filename = 'magnetic_calculation_hopes_params.h5';

info = h5info(filename, '/');
allens_spec_data = struct; 
for i = 1:length(info.Datasets)
    % Get the dataset name
    datasetName = info.Datasets(i).Name;
    % Load the dataset into a variable
    datasetData = h5read(filename, ['/', datasetName]);
    % Optionally, store the data in a struct (to make it easy to access)
    my_mag_data.(datasetName) = datasetData;
end
%% load allens data
filename = 'magnetic_calculation_Allens_params.h5';

info = h5info(filename, '/');
allens_mag_data = struct; 
for i = 1:length(info.Datasets)
    % Get the dataset name
    datasetName = info.Datasets(i).Name;
    % Load the dataset into a variable
    datasetData = h5read(filename, ['/', datasetName]);
    % Optionally, store the data in a struct (to make it easy to access)
    allens_mag_data.(datasetName) = datasetData;
end
%%
% Plot data
figure;
subplot(2,3,1)
hold on; grid on; box on; 
% Plot experimental data
plot(my_mag_data.magnetization20K(:,1), my_mag_data.magnetization20K(:,2)/1.37, 'o', 'DisplayName', '20K MPMS data');
plot(my_mag_data.magnetization6K(:,1), my_mag_data.magnetization6K(:,2)/1.37, 'o', 'DisplayName', '6K MPMS data');
plot(my_mag_data.magnetization2K(:,1), my_mag_data.magnetization2K(:,2)/1.37, 'o', 'DisplayName', '2K MPMS data');
plot(my_mag_data.CESMHdata(:,7)./1e4, my_mag_data.CESMHdata(:,8), 'b.', 'DisplayName', 'From Allens paper');

% Plot MFT data
H = horzcat(linspace(0,1,100), linspace(1.01,10, 100));
plot(H, my_mag_data.tempMagC(:,11), '-', 'DisplayName', 'MFT 2K');
plot(H, my_mag_data.tempMagC(:,12), '-', 'DisplayName', 'MFT 6K');
plot(H, my_mag_data.tempMagC(:,13), '-', 'DisplayName', 'MFT 20K');

% Set axis limits (optional)
xlim([0, 6]);
ylim([0, 7]);

% Add legend and labels
legend('show');
title('C Magnetization');
xlabel('Field (T)');
ylabel('Magnetization (uB/Er)');
hold off;

% %% get ab plane
% Load data from HDF5 file

% Plot the data
subplot(2,3,4);  grid on; box on; 
hold on;

% MPMS data
plot(my_mag_data.magnetizationAB20K(:,1), my_mag_data.magnetizationAB20K(:,2)/1.37, 'o', 'DisplayName', '20K MPMS data');
plot(my_mag_data.magnetizationAB6K(:,1), my_mag_data.magnetizationAB6K(:,2)/1.37, 'o', 'DisplayName', '6K MPMS data');
plot(my_mag_data.magnetizationAB2K(:,1), my_mag_data.magnetizationAB2K(:,2)/1.37, 'o', 'DisplayName', '2K MPMS data');

% Allen's paper data
plot(my_mag_data.CESMHdata(:, 1) / 1e4, my_mag_data.CESMHdata(:, 2), 'b.', 'DisplayName', 'From Allen''s paper');

% MFT and calculated magnetization
plot(H, my_mag_data.tempMagB(:,11), '-', 'DisplayName', 'MFT 2K');
plot(H, my_mag_data.tempMagB(:,12), '-', 'DisplayName', 'MFT 6K');
plot(H, my_mag_data.tempMagB(:,13), '-', 'DisplayName', 'MFT 20K');

xlim([0, 6]);
ylim([0, 7]);

title('AB plane magnetization');
xlabel('Field (T)');
ylabel('Magnetization (\mu_B/Er)');
legend show;
hold off;


% %% dmdh with data

% Plot settings
labels = my_mag_data.dmdhLabels;
n = length(labels);
colors = cool(n); % wanted inferno, can't find it :(

subplot(1,3,2); grid on; box on; 
hold on;
dmdh = my_mag_data.dmdhC(:, [1,3:11]);
for i = 1:n
    % Normalize yArrs and dmdH for plotting
    dmData = my_mag_data.dmdhData{i};
    dmData = dmData / max(dmData);
    dm = dmdh(:,i);
    % dm = dm + abs(min(dm));
    dm = dm / max(dm);

    % Offset for better visualization
    offset = (i - 1) * 0.5;

    % Plot experimental data
    plot(my_mag_data.dmdhField{i}, dmData + offset, 'DisplayName', labels{i}, 'Color', colors(i, :));

    % Plot calculated data
    plot(H, dm + offset, '--', 'DisplayName', [labels{i} , ' (calc)'], 'Color', colors(i, :));
    text(9, offset+0.2, labels{i}, 'FontSize', 9);
end

title('dM/dH ');
ylabel('dM/dH (arb)');
xlabel('Field (T)');
hold off;

% %% make integrated scm1 data

integratedMag = {}; 
% now we integrate yArrs
for i = 1:n
    x = my_mag_data.dmdhField{i}; 
    y = my_mag_data.dmdhData{i};
    y = y/max(y); 
    
    % first sort x
    [x, inds] = sort(x);
    y = y(inds); 
    temporary = cumsum(y); 
    integratedMag{end+1} = temporary; 
end

% Plotting
n = length(labels);
colors = cool(n+1);

subplot(1,3,3); grid on; box on; 
hold on;
for i = [1, 6, 9]
    % Sort and normalize x and y data
    x = my_mag_data.dmdhField{i}; 
    y = my_mag_data.dmdhData{i};
    [x, sortIdx] = sort(x);
    y = y(sortIdx) / max(y);
    intMag = integratedMag{i} / max(integratedMag{i}); % Normalize integratedMag
    
    % Plot data
    plot(x, intMag , 'DisplayName', labels{i}, 'Color', colors(i+1, :));
    plot(H, my_mag_data.tempMagC(:, i) / max(my_mag_data.tempMagC(:,i)), '--', 'Color', colors(i+1, :));
end

% Add labels, title, and legend
title('Integrated dM/dH');
xlabel('Field (T)');
ylabel('Magnetization (arb)');
legend('show');
hold off;

%% same figs for allen
fig = figure;
ax1 = subplot(2,3,1);
title('B-axis IR')
hold on; 
contourf(allens_spec_data.IR_B_field, allens_spec_data.IR_B_wavenums, allens_spec_data.IR_dataB', 100, 'LineStyle', 'none')
for i = 2:16%size(arrC, 2)
    plot(allens_spec_data.calc_field, allens_spec_data.linesB(:,i),  'r:', 'LineWidth', 0.5);
end
clim([0 1])
ylim([0 100])
xlim([0 17.5])
colormap(ax1, flipud(parula))
ax4 = subplot(2,3,4);
hold on; 
contourf(allens_spec_data.calc_field, allens_spec_data.calc_wavenums, allens_spec_data.simulated_IR_B, 100,'LineStyle', 'none');
colormap(ax4, flipud(parula))
clim([0 1])

for i = 2:16%size(arrC, 2)
    plot(allens_spec_data.calc_field, allens_spec_data.linesB(:,i),  'r:', 'LineWidth', 0.5);
end

clim([0 1])
ylim([0 100])
xlim([0 17.5])


ax2 = subplot(2,3,2);hold on; box on; 
contourf(allens_spec_data.IR_C_field, allens_spec_data.IR_C_wavenums, allens_spec_data.IR_dataC', 100, 'LineStyle', 'none')
for i = 2:16%size(arrC, 2)
    plot(allens_spec_data.calc_field,allens_spec_data.linesC(:,i),  'r:', 'LineWidth', 0.5);
end
title('C-axis IR')
colormap(ax2, flipud(parula))
clim([0 1])
ylim([0 100])
xlim([0 17.5])

ax5 = subplot(2,3,5);
hold on; box on; 
contourf(allens_spec_data.calc_field, allens_spec_data.calc_wavenums, allens_spec_data.simulated_IR_C, 100, 'LineStyle', 'none');
colormap(ax5, flipud(parula))
for i = 2:16%size(arrC, 2)
    plot(allens_spec_data.calc_field,allens_spec_data.linesC(:,i),  'r:', 'LineWidth', 0.5);
end
clim([0 1])
ylim([0 100])
xlim([0 17.5])




% fig = figure;
ax3 = subplot(2,3,3);box on; hold on; 
contourf(allens_spec_data.raman_field, allens_spec_data.raman_wavenums, allens_spec_data.ramanData', 100, 'LineStyle', 'none')
for i = 2:16%size(arrC, 2)
    plot(allens_spec_data.calc_field,allens_spec_data.linesC(:,i),  'r:', 'LineWidth', 0.5);
end
ylim([0 100])
xlim([0 17.5])
colormap(ax3, 'cool')
title ('C-axis Raman')
ax6 = subplot(2,3,6);
hold on; box on; 
contourf(allens_spec_data.calc_field, allens_spec_data.calc_wavenums, allens_spec_data.simulated_raman, 100, 'LineStyle', 'none');
colormap(ax6, 'cool')
% let's add some lines
for i = 2:16%size(arrC, 2)
    plot(allens_spec_data.calc_field,allens_spec_data.linesC(:,i),  'r:', 'LineWidth', 0.5);
end
clim([0 1])
xlim([0 17.5])
ylim([0 100])




ax3.Position(2)=ax6.Position(2)+ax6.Position(4)+.02;
ax1.Position(2)=ax4.Position(2)+ax4.Position(4)+.02;
% 
% ax4.Position(2)=ax6.Position(2)+ax6.Position(4);
ax2.Position(2)=ax5.Position(2)+ax5.Position(4)+.02;
% 
% 
% ax4.Position(1)=ax3.Position(1)+ax3.Position(3);
ax2.Position(1)=ax1.Position(1)+ax1.Position(3);
ax5.Position(1)=ax4.Position(1)+ax4.Position(3);
linkaxes([ax1, ax2, ax3, ax4, ax5, ax6], 'xy')
% ylabel('Energy [cm^{-1}]')
% xlabel('H [t]')
set(ax2,'Yticklabel',[]) 
set(ax1,'Xticklabel',[])
set(ax2,'Xticklabel',[])

set(ax5,'Yticklabel',[]) 
% set(ax5,'Xticklabel',[])


set(ax3,'Xticklabel',[])

%% make my lines mft vs no mft

figure; box on; grid on; hold on; 
for i = 2:16%size(arrC, 2)
    if i == 2
        plot(my_spec_data.calc_field, my_spec_data.linesC(:,i),  'r', 'LineWidth', 0.5, 'DisplayName','MFT');
        plot(my_spec_data.calc_field, my_spec_data.linesC_nomft(:,i),  'b', 'LineWidth', 0.5, 'DisplayName','no MFT');
    else
        plot(my_spec_data.calc_field, my_spec_data.linesC(:,i),  'r', 'LineWidth', 0.5);
        plot(my_spec_data.calc_field, my_spec_data.linesC_nomft(:,i),  'b', 'LineWidth', 0.5);
    end
end
xlabel('H [T]')
ylabel('Energy [cm^{-1}')
legend()

%% make allen lines mft vs no mft

figure; box on; grid on; hold on; 
for i = 2:16%size(arrC, 2)
    if i == 2
        plot(allens_spec_data.calc_field,allens_spec_data.linesC(:,i),  'r', 'LineWidth', 0.5, 'DisplayName','MFT');
        plot(allens_spec_data.calc_field,allens_spec_data.linesC_nomft(:,i),  'b', 'LineWidth', 0.5, 'DisplayName','no MFT');
    else
        plot(allens_spec_data.calc_field,allens_spec_data.linesC(:,i),  'r', 'LineWidth', 0.5);
        plot(allens_spec_data.calc_field,allens_spec_data.linesC_nomft(:,i),  'b', 'LineWidth', 0.5);
    end
end
xlabel('H [T]')
ylabel('Energy [cm^{-1}')
legend()

%% okay, so now let's make a susceptibility fig
% load my magnetic data
filename = 'susceptibility_calculation_hopes_params.h5';

info = h5info(filename, '/');
for i = 1:length(info.Datasets)
    % Get the dataset name
    datasetName = info.Datasets(i).Name;
    % Load the dataset into a variable
    datasetData = h5read(filename, ['/', datasetName]);
    % Optionally, store the data in a struct (to make it easy to access)
    my_sus_data.(datasetName) = datasetData;
end
%% load allens data
filename = 'susceptibility_calculation_Allens_params.h5';

info = h5info(filename, '/'); 
for i = 1:length(info.Datasets)
    % Get the dataset name
    datasetName = info.Datasets(i).Name;
    % Load the dataset into a variable
    datasetData = h5read(filename, ['/', datasetName]);
    % Optionally, store the data in a struct (to make it easy to access)
    allens_sus_data.(datasetName) = datasetData;
end
%% plot temp dependent sus
figure; hold on; box on; grid on; 
for i = 1: length(my_sus_data.susC(1,:))
    plot(my_sus_data.temps, 1./my_sus_data.susC(:,i), 'DisplayName', [num2str(my_sus_data.fieldVals(i)), 'T'])
end
% add data
for i= 1:length(my_sus_data.data_sus_C)
    plot(my_sus_data.data_temps_C{i}, my_sus_data.data_sus_C{i}*1.37, 'b.', 'DisplayName', my_sus_data.clabels{i}); 

end
title('Susceptibility, c-axis, my params')
legend(); 
xlabel('Temperature [K]'); 
ylabel('1/\chi')
%% ab plane
figure; hold on; box on; grid on; 
for i = 1: length(my_sus_data.susB(1,:))
    plot(my_sus_data.temps, 1./my_sus_data.susB(:,i), 'DisplayName', [num2str(my_sus_data.fieldVals(i)), 'T'])
end
% add data
for i= 1:length(my_sus_data.data_sus_AB)
    plot(my_sus_data.data_temps_AB{i}, my_sus_data.data_sus_AB{i}*1.37, 'b.', 'DisplayName', my_sus_data.blabels{i}); 

end
title('Susceptibility, ab plane, my params')
legend(); 
xlabel('Temperature [K]'); 
ylabel('1/\chi')

%% Jz test
% Load data from HDF5
fileName = 'jjz_test.h5';
H = h5read(fileName, '/H');
Jvals = h5read(fileName, '/Jvals');
tempMagC = h5read(fileName, '/tempMagC');
Mdata2K_x = h5read(fileName, '/Mdata2K_x');
Mdata2K_y = h5read(fileName, '/Mdata2K_y');
temperature = h5readatt(fileName, '/', 'temperature');

% Plot magnetization curves
figure;
hold on;
grid on;

% Plot calculated magnetization curves
numCurves = size(tempMagC, 2);
for i = 1:numCurves
    plot(H, tempMagC(:, i), 'DisplayName', sprintf('JJz = %.3g', Jvals(i)));
end

% Overlay experimental 2K MPMS data
plot(Mdata2K_x, Mdata2K_y / 1.37, 'bo', 'DisplayName', '2K MPMS data');

% Customize plot
xlim([0, 8]);
title(sprintf('Magnetization vs Field at T = %.2f', temperature));
xlabel('Field (T)');
ylabel('Magnetization (\mu_B/Er)');
legend('show');
hold off;

%% make temp dependence for M_c vs H
figure;
hold on;
grid on;

% Plot calculated magnetization curves
numCurves = size(my_mag_data.tempMagC, 2);
for i = 1:numCurves
    plot(my_mag_data.H, my_mag_data.tempMagC(:, i), 'DisplayName', num2str(my_mag_data.temps(i)));
end
legend()
xlabel('H[T]')
ylabel('M \mu_B/Er')
title('My params')

%% same for ab
figure;
hold on;
grid on;

% Plot calculated magnetization curves
numCurves = size(my_mag_data.tempMagB, 2);
for i = 1:numCurves
    plot(my_mag_data.H, my_mag_data.tempMagB(:, i), 'DisplayName', num2str(my_mag_data.temps(i)));
end
legend()
xlabel('H[T]')
ylabel('M \mu_B/Er')
title('My params')