%% make c axis magnetization fig
% Load data from the HDF5 file
info = h5info('magnetization_data_calculation_and_allen_c_axis.h5', '/');
data = struct; 

for i = 1:length(info.Datasets)
    % Get the dataset name
    datasetName = info.Datasets(i).Name;
    
    % Load the dataset into a variable
    datasetData = h5read('magnetization_data_calculation_and_allen_c_axis.h5', ['/', datasetName]);
    
    % Optionally, store the data in a struct (to make it easy to access)
    data.(datasetName) = datasetData;

end

% Extract the data
magF = data.magF;
MFTField = data.MFTField;
temperature = data.temperature;
field = data.field;
myCaxisMagnetization = data.myCaxisMagnetization;
allenCaxisMagnetization = data.allenCaxisMagnetization;
myMFTCaxis = data.myMFTCaxis;
allenMFTCaxis = data.allenMFTCaxis;
CESMHdata = data.CESMHdata;
CESMTdata = data.CESMTdata;

% Plot data

fig = figure(); 
fontsize(fig, 24, "points")
hold on; grid on;

% Plot CESMHdata (experimental data)
plot(CESMHdata(7,:)./ 1e4, CESMHdata(8,:), 'b.', 'DisplayName', 'Data ($H \\parallel c$)');

% Plot MFT and Allen MFT curves 
plot(MFTField, myMFTCaxis, '-', 'DisplayName', 'Raman fit B params');
plot(MFTField, allenMFTCaxis, '-', 'DisplayName', 'Neutron fit B params');

% Plot magnetization curves without MFT correction
plot(magF, myCaxisMagnetization, '--', 'DisplayName', 'Raman fit Bparams, no MFT');
plot(magF, allenCaxisMagnetization, '--', 'DisplayName', 'Neutron fit Bparams, no MFT');

% Optional: Adjust axis limits if needed
xlim([0, 6]);
ylim([0, 10]);

% Add legend and labels
legend('show');
title('C Magnetization');
xlabel('Field (T)');
ylabel('Magnetization (uB/Er)');

%% make mpms magnetization fig
info = h5info('magnetization_data_caluclation_mpms_data_c_axis.h5');

% Load datasets
magF_mpms = h5read('magnetization_data_caluclation_mpms_data_c_axis.h5', '/magF_mpms');
MFTField_mpms = h5read('magnetization_data_caluclation_mpms_data_c_axis.h5', '/MFTField_mpms');
magnetization2K = h5read('magnetization_data_caluclation_mpms_data_c_axis.h5', '/magnetization2K');
magnetization6K = h5read('magnetization_data_caluclation_mpms_data_c_axis.h5', '/magnetization6K');
magnetization20K = h5read('magnetization_data_caluclation_mpms_data_c_axis.h5', '/magnetization20K');
MFT2K = h5read('magnetization_data_caluclation_mpms_data_c_axis.h5', '/MFT2K');
MFT6K = h5read('magnetization_data_caluclation_mpms_data_c_axis.h5', '/MFT6K');
MFT20K = h5read('magnetization_data_caluclation_mpms_data_c_axis.h5', '/MFT20K');
Mdata2K = h5read('magnetization_data_caluclation_mpms_data_c_axis.h5', '/Mdata2K');
Mdata6K = h5read('magnetization_data_caluclation_mpms_data_c_axis.h5', '/Mdata6K');
Mdata20K = h5read('magnetization_data_caluclation_mpms_data_c_axis.h5', '/Mdata20K');
CESMHdata = h5read('magnetization_data_caluclation_mpms_data_c_axis.h5', '/CESMHdata');

% Plot data
figure;
hold on; grid on; box on; 
% Plot experimental data
plot(Mdata20K(:,1), Mdata20K(:,2)/1.35, 'o', 'DisplayName', '20K MPMS data');
plot(Mdata6K(:,1), Mdata6K(:,2)/1.35, 'o', 'DisplayName', '6K MPMS data');
plot(Mdata2K(:,1), Mdata2K(:,2)/1.35, 'o', 'DisplayName', '2K MPMS data');
plot(CESMHdata(:,7)./1e4, CESMHdata(:,8), 'b.', 'DisplayName', 'From Allens paper');

% Plot MFT data
plot(MFTField_mpms, MFT2K, '-', 'DisplayName', 'MFT 2K');
plot(MFTField_mpms, MFT6K, '-', 'DisplayName', 'MFT 6K');
plot(MFTField_mpms, MFT20K, '-', 'DisplayName', 'MFT 20K');

% Plot magnetization curves without MFT correction
plot(magF_mpms, magnetization2K, '--', 'DisplayName', '2K, no MFT');
plot(magF_mpms, magnetization6K, '--', 'DisplayName', '6K, no MFT');
plot(magF_mpms, magnetization20K, '--', 'DisplayName', '20K, no MFT');

% Set axis limits (optional)
xlim([-6, 6]);
ylim([-10, 10]);

% Add legend and labels
legend('show');
title('C Magnetization');
xlabel('Field (T)');
ylabel('Magnetization (uB/Er)');
hold off;

%% make temp dependent m vs h
% Load HDF5 data
temps = h5read('M_vs_H_temperature_dependence.h5', '/temps');
MFTField = h5read('M_vs_H_temperature_dependence.h5', '/MFTField');
tempMag = h5read('M_vs_H_temperature_dependence.h5', '/tempMag');

% Create the plot
figure;
hold on;
cmap = jet(length(temps)); % Colormap for consistency

for i = 1:length(temps)
    plot(MFTField, tempMag(:, i), 'DisplayName', sprintf('%.3f K', temps(i)), 'Color', cmap(i, :));
end

legend show;
xlim([0 9]);
title({'C axis magnetization MFT', ...
       'calculated from Raman fit B params', ...
       'test B60 = <INSERT VALUE>'});
xlabel('Field (T)');
ylabel('Magnetization');
hold off;

%% make dmdH
% Load HDF5 data
temps = h5read('dMdH_temperature_dependence.h5', '/temps');
MFTField = h5read('dMdH_temperature_dependence.h5', '/MFTField');
dmdH = h5read('dMdH_temperature_dependence.h5', '/dmdH');

% Set up the colormap for consistency
n = length(temps);
colors = jet(n);

% Create the plot
figure;
hold on;

for i = 1:n
    plot(MFTField, dmdH(:, i), 'DisplayName', sprintf('%.3f K', temps(i)), 'Color', colors(i, :));
end

legend show;
xlabel('Field (T)');
ylabel('dM/dH');
set(gca, 'YScale', 'log'); % Set y-axis to logarithmic scale
title({'C axis dM/dH', 'calculated from Raman fit B params'});

% xlim([0 7]); % Set x-axis limits
% ylim([-1 12]); % Set y-axis limits

hold off;

%% dmdh with data

% Load HDF5 data
temps = h5read('dMdH_temperature_dependence_withdata.h5', '/temps');
MFTField = h5read('dMdH_temperature_dependence_withdata.h5', '/MFTField');
dmdH = h5read('dMdH_temperature_dependence_withdatas.h5', '/dmdH');
xArrs = h5read('dMdH_temperature_dependence_withdata.h5', '/xArrs');
yArrs = h5read('dMdH_temperature_dependence_withdata.h5', '/yArrs');
labels = h5read('dMdH_temperature_dependence_withdata.h5', '/labels');
% labels = cellstr(char(labels')); % Convert labels to cell array of strings

% Plot settings
n = length(labels);
colors = jet(n); % wanted inferno, can't find it :(

figure;
hold on;

for i = 1:n
    % Normalize yArrs and dmdH for plotting
    y = yArrs{i} / max(yArrs{i});
    dm = log(dmdH(:, i));
    dm = dm + abs(min(dm));
    dm = dm / max(dm);
    
    % Offset for better visualization
    offset = (i - 1) * 0.5;
    
    % Plot experimental data
    plot(xArrs{i}, y + offset, 'DisplayName', labels{i}, 'Color', colors(i, :));
    
    % Plot calculated data
    plot(MFTField, dm + offset, '--', 'DisplayName', [labels{i} , ' (calc)'], 'Color', colors(i, :));
end

title('dM/dH from SCM1 \n calculated dM/dH in dotted line');
ylabel('dM/dH (arb)');
xlabel('Field (T)');
legend show;
hold off;

%% susceptibility

% Load HDF5 data
xArrs = h5read('susceptibility_wide_temp_range.h5', '/xArrs');
yArrs = h5read('susceptibility_wide_temp_range.h5', '/yArrs');
labels = h5read('susceptibility_wide_temp_range.h5', '/labels');
temps = h5read('susceptibility_wide_temp_range.h5', '/temps');
myinv01T = h5read('susceptibility_wide_temp_range.h5', '/myinv01T');
neutroninv01T = h5read('susceptibility_wide_temp_range.h5', '/neutroninv01T');
myinv1T = h5read('susceptibility_wide_temp_range.h5', '/myinv1T');
neutroninv1T = h5read('susceptibility_wide_temp_range.h5', '/neutroninv1T');
myinv3T = h5read('susceptibility_wide_temp_range.h5', '/myinv3T');
neutroninv3T = h5read('susceptibility_wide_temp_range.h5', '/neutroninv3T');
susinvPCF = h5read('susceptibility_wide_temp_range.h5', '/susinvPCF');


% Create the plot
figure;
hold on;

% Plot experimental data
for i = 1:length(labels)
    x = sort(xArrs{i}); % Sort x data
    [~, inds] = sort(xArrs{i}); % Sorting indices
    y = yArrs{i}(inds); % Sort y accordingly
    plot(x, y * 1.35, 'DisplayName', labels{i});
end

% Plot CESMT data (replace with your actual data and constant)
plot(CESMTdata(:, 12), 1 ./ CESMTdata(:, 12) * SCF, 'DisplayName', 'C-axis data from Allen''s paper');

% Plot simulated data for different fields
plot(temps, myinv01T, '--', 'DisplayName', 'Raman B params MFT 0.1T');
plot(temps, neutroninv01T, '-.', 'DisplayName', 'Neutrons B params MFT 0.1T');
plot(temps, myinv1T, '--', 'DisplayName', 'Raman B params MFT 1T');
plot(temps, neutroninv1T, '-.', 'DisplayName', 'Neutrons B params MFT 1T');
plot(temps, myinv3T, '--', 'DisplayName', 'Raman B params MFT 3T');
plot(temps, neutroninv3T, '-.', 'DisplayName', 'Neutrons B params MFT 3T');
plot(temps, susinvPCF, '--', 'DisplayName', 'Raman B params no MFT, 0.1T');

% Customize the plot
title('Calculated MFT Susceptibility');
xlabel('Temperature (K)');
ylabel('1/\chi');
legend show;
xlim([0, 200]);
ylim([0, 10]);
hold off;

%% make integrated scm1 data

% Load HDF5 data
xArrs = h5read('dMdH_temperature_dependence_withdata.h5', '/xArrs');
yArrs = h5read('dMdH_temperature_dependence_withdata.h5', '/yArrs');
labels = h5read('dMdH_temperature_dependence_withdata.h5', '/labels');
magF = h5read('dMdH_temperature_dependence_withdata.h5', '/magF');
tempMag = h5read('dMdH_temperature_dependence_withdata.h5', '/tempMag');

integratedMag = {}; 
% now we integrate yArrs
for i = 1:length(yArrs)
    x = xArrs{i}; 
    y = yArrs{i};
    y = y/max(y); 
    
    % first sort x
    [x, inds] = sort(x);
    y = y(inds); 
    temporary = cumsum(y); 
    integratedMag{end+1} = temporary; 
end

% Plotting
n = length(labels);
colors = jet(n);

figure;
hold on;
for i = 1:n
    % Sort and normalize x and y data
    x = xArrs{i};
    y = yArrs{i};
    [x, sortIdx] = sort(x);
    y = y(sortIdx) / max(y);
    intMag = integratedMag{i} / max(integratedMag{i}); % Normalize integratedMag
    
    % Plot data
    plot(x, intMag + i * 0.5, 'DisplayName', labels{i}, 'Color', colors(i, :));
    plot(magF, tempMag(i) / max(tempMag(i)) + i * 0.5, '--', 'Color', colors(i, :));
end

% Add labels, title, and legend
title('Integrated chi(H) \n Numerically Integrated from SCM1 Data \n Calculated Curve in Dotted');
xlabel('Field (T)');
ylabel('Magnetization (arb)');
legend('show');
hold off;

%% low temp chi(H)
% Load data from HDF5 file
xArrs = h5read('susceptibility_low_temp.h5', '/xArrs');
yArrs = h5read('susceptibility_low_temp.h5', '/yArrs');
labels = h5read('susceptibility_low_temp.h5', '/labels');
fields = h5read('susceptibility_low_temp.h5', '/fields');
tempArr = h5read('susceptibility_low_temp.h5', '/tempArr');
susArr = h5read('susceptibility_low_temp.h5', '/susArr');

% Set up colormap
nLabels = length(labels);
dataColors = cool(nLabels); % Colormap for experimental data
nFields = length(fields);
calcColors = cool(nFields); % Colormap for calculated data

% Create figure
figure;
hold on;

% Plot experimental data
yyaxis left;
for i = 1:nLabels
    x = xArrs{i};
    y = yArrs{i};
    
    % Sort data
    [x, sortIdx] = sort(x);
    y = y(sortIdx);
    
    % Normalize y data
    y = y / y(end);
    
    % Check for repeated labels to use consistent colors
    if i > 1 && strcmp(labels{i}, labels{i-1})
        plotColor = dataColors(i-1, :);
    else
        plotColor = dataColors(i, :);
    end
    
    plot(x, y, 'o', 'DisplayName', labels{i}, 'Color', plotColor);
end
ylabel('Data Chi (arbitrarily scaled)');
xlabel('Field (T)');

% Plot calculated susceptibility
yyaxis right;
for i = 1:nFields
    sus = susArr(:, i);
    plot(tempArr, sus, '--', 'DisplayName', sprintf('%.1fT', fields(i)), 'Color', calcColors(i, :));
end
ylabel('Calculated Chi');

% Add title and legend
title('Field-Dependent Susceptibility Analysis');
legend('show');
hold off;

%% AB plane susceptibility
% Load data from HDF5 file
xArrs = h5read('susceptibility_AB_plane.h5', '/xArrs');
yArrs = h5read('susceptibility_AB_plane.h5', '/yArrs');
labels = h5read('susceptibility_AB_plane.h5', '/labels');
temps = h5read('susceptibility_AB_plane.h5', '/temps');
myinv01T = h5read('susceptibility_AB_plane.h5', '/myinv01T');
myinv1T = h5read('susceptibility_AB_plane.h5', '/myinv1T');
myinv6T = h5read('susceptibility_AB_plane.h5', '/myinv6T');
neutroninv = h5read('susceptibility_AB_plane.h5', '/neutroninv');
susinvPCF = h5read('susceptibility_AB_plane.h5', '/susinvPCF');

% Create figure
figure;
hold on;

% Plot experimental data
for i = 1:length(labels)
    x = xArrs{i};
    y = yArrs{i};
    plot(x, y.*1.35, 'DisplayName', labels{i});
end

% Plot calculated susceptibility
plot(temps, myinv01T, '--', 'DisplayName', 'Raman B params MFT, 0.1T');
plot(temps, myinv1T, '--', 'DisplayName', 'Raman B params MFT, 1T');
plot(temps, myinv6T, '--', 'DisplayName', 'Raman B params MFT, 6T');
plot(temps, neutroninv, '-.', 'DisplayName', 'neutrons B params MFT');
plot(temps, susinvPCF, '--', 'DisplayName', 'Raman B params no MFT');

% Add labels, legend, and title
title('calculated MFT susceptibility at 0.1T (AB plane)');
xlabel('Temperature (K)');
ylabel('1/Chi');
legend('show');
xlim([0, 200]);
ylim([0, 10]);
hold off;