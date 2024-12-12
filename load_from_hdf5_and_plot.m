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
info = h5info('magnetization_data_calculation_mpms_data_c_axis.h5');

% Load datasets
magF_mpms = h5read('magnetization_data_calculation_mpms_data_c_axis.h5', '/magF_mpms');
MFTField_mpms = h5read('magnetization_data_calculation_mpms_data_c_axis.h5', '/MFTField_mpms');
magnetization2K = h5read('magnetization_data_calculation_mpms_data_c_axis.h5', '/magnetization2K');
magnetization6K = h5read('magnetization_data_calculation_mpms_data_c_axis.h5', '/magnetization6K');
magnetization20K = h5read('magnetization_data_calculation_mpms_data_c_axis.h5', '/magnetization20K');
MFT2K = h5read('magnetization_data_calculation_mpms_data_c_axis.h5', '/MFT2K');
MFT6K = h5read('magnetization_data_calculation_mpms_data_c_axis.h5', '/MFT6K');
MFT20K = h5read('magnetization_data_calculation_mpms_data_c_axis.h5', '/MFT20K');
Mdata2K = h5read('magnetization_data_calculation_mpms_data_c_axis.h5', '/Mdata2K');
Mdata6K = h5read('magnetization_data_calculation_mpms_data_c_axis.h5', '/Mdata6K');
Mdata20K = h5read('magnetization_data_calculation_mpms_data_c_axis.h5', '/Mdata20K');
CESMHdata = h5read('magnetization_data_calculation_mpms_data_c_axis.h5', '/CESMHdata');

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
figure; grid on; box on; 
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
% figure; grid on; box on; 
hold on;

for i = 1:n
    plot(MFTField, dmdH(:, i), 'DisplayName', sprintf('a %.3f K', temps(i)), 'Color', colors(i, :));
end

legend show;
xlabel('Field (T)');
ylabel('dM/dH');
set(gca, 'YScale', 'log'); % Set y-axis to logarithmic scale
title({'C axis dM/dH', 'calculated from Raman fit B params'});

% xlim([0 7]); % Set x-axis limits
% ylim([-1 12]); % Set y-axis limits

% hold off;

%% dmdh with data

% Load HDF5 data
temps = h5read('dmdh_with_scm1_data.h5', '/temps');
MFTField = h5read('dmdh_with_scm1_data.h5', '/MFTField');
dmdH = h5read('dmdh_with_scm1_data.h5', '/dmdH');
xArrs = h5read('dmdh_with_scm1_data.h5', '/xArrs');
yArrs = h5read('dmdh_with_scm1_data.h5', '/yArrs');
labels = h5read('dmdh_with_scm1_data.h5', '/labels');
% labels = cellstr(char(labels')); % Convert labels to cell array of strings

% Plot settings
n = length(labels);
colors = jet(n); % wanted inferno, can't find it :(

figure; grid on; box on; 
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

%% now plot dmdh with no offset
fname = '/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/calculated_data/dMdH_temperature_dependence_withdata.h5';
temps = h5read(fname, '/temps');
MFTField = h5read(fname, '/MFTField');
dmdH = h5read(fname, '/dmdH');
xArrs = h5read(fname, '/xArrs');
yArrs = h5read(fname, '/yArrs');
labels = h5read(fname, '/labels');
% labels = cellstr(char(labels')); % Convert labels to cell array of strings

% Plot settings
n = length(labels);
colors = jet(n); % wanted inferno, can't find it :(

figure; grid on; box on; 
hold on;

for i = 1:n
    % Normalize yArrs and dmdH for plotting
    y = yArrs{i}; %/ max(yArrs{i}); %normalized data dmdh
    dm = dmdH(:, i); % calcluated dmdh
    % dm = dm / max(dm);
    
    % Plot experimental data
    plot(xArrs{i}, y , 'DisplayName', ['data ', labels{i}], 'Color', colors(i, :));
    
    % % Plot calculated data
    % plot(MFTField, dm, '--', 'DisplayName', [labels{i} , ' (calc)'], 'Color', colors(i, :));
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
myinv0T = h5read('susceptibility_wide_temp_range.h5', '/myinv0T');
neutroninv0T = h5read('susceptibility_wide_temp_range.h5', '/neutroninv0T');
myinv01T = h5read('susceptibility_wide_temp_range.h5', '/myinv01T');
neutroninv01T = h5read('susceptibility_wide_temp_range.h5', '/neutroninv01T');
myinv1T = h5read('susceptibility_wide_temp_range.h5', '/myinv1T');
neutroninv1T = h5read('susceptibility_wide_temp_range.h5', '/neutroninv1T');
myinv3T = h5read('susceptibility_wide_temp_range.h5', '/myinv3T');
neutroninv3T = h5read('susceptibility_wide_temp_range.h5', '/neutroninv3T');
susinvPCF = h5read('susceptibility_wide_temp_range.h5', '/susinvPCF');
CESMTdata = h5read('susceptibility_wide_temp_range.h5', '/CESMTdata');

Na = 6.02214076e23 ;
SCF = 1/(1.07828221e24/Na);
% Create the plot
figure; grid on; box on; 
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
plot(temps, myinv0T, '--', 'DisplayName', 'Raman B params MFT 0T');
plot(temps, neutroninv0T, '-.', 'DisplayName', 'Neutrons B params MFT 0T');
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

figure; grid on; box on; 
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
    plot(magF, tempMag(:, i) / max(tempMag(:,i)) + i * 0.5, '--', 'Color', colors(i, :));
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
figure; grid on; box on; 
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
mysus0T = h5read('susceptibility_AB_plane.h5', '/mysus0T');
mysus01T = h5read('susceptibility_AB_plane.h5', '/mysus01T');
mysus1T = h5read('susceptibility_AB_plane.h5', '/mysus1T');
mysus6T = h5read('susceptibility_AB_plane.h5', '/mysus6T');
neutronsus = h5read('susceptibility_AB_plane.h5', '/neutronsus01T');
susPCF = h5read('susceptibility_AB_plane.h5', '/sussusPCF');

% Create figure
figure; grid on; box on; 
hold on;

% Plot experimental data
for i = 1:length(labels)
    x = xArrs{i};
    y = yArrs{i};
    plot(x, y.*1.35, 'DisplayName', [labels{i}, 'MPMS data']);
end

% Plot calculated susceptibility
plot(temps, 1./mysus0T, '--', 'DisplayName', 'Raman B params MFT, 0T');
plot(temps, 1./mysus01T, '--', 'DisplayName', 'Raman B params MFT, 0.1T');
plot(temps, 1./mysus1T, '--', 'DisplayName', 'Raman B params MFT, 1T');
plot(temps, 1./mysus6T, '--', 'DisplayName', 'Raman B params MFT, 6T');
plot(temps, 1./neutronsus, '-.', 'DisplayName', 'neutrons B params MFT');
plot(temps, -1./susPCF, '--', 'DisplayName', 'Raman B params no MFT');

% Add labels, legend, and title
title('calculated MFT susceptibility at 0.1T (AB plane)');
xlabel('Temperature (K)');
ylabel('1/Chi');
legend('show');
xlim([0, 200]);
ylim([0, 10]);
hold off;

%% get mpms ab plane magnetization
% Load data from HDF5 file
magF_mpms = h5read('magnetization_data_calculation_mpms_data_AB_plane.h5', '/magF_mpms');
MFTField_mpms = h5read('magnetization_data_calculation_mpms_data_AB_plane.h5', '/MFTField_mpms');
magnetization2K = h5read('magnetization_data_calculation_mpms_data_AB_plane.h5', '/magnetization2K');
magnetization6K = h5read('magnetization_data_calculation_mpms_data_AB_plane.h5', '/magnetization6K');
magnetization20K = h5read('magnetization_data_calculation_mpms_data_AB_plane.h5', '/magnetization20K');
MFT2K = h5read('magnetization_data_calculation_mpms_data_AB_plane.h5', '/MFT2K');
MFT6K = h5read('magnetization_data_calculation_mpms_data_AB_plane.h5', '/MFT6K');
MFT20K = h5read('magnetization_data_calculation_mpms_data_AB_plane.h5', '/MFT20K');
Mdata2K = h5read('magnetization_data_calculation_mpms_data_AB_plane.h5', '/Mdata2K');
Mdata6K = h5read('magnetization_data_calculation_mpms_data_AB_plane.h5', '/Mdata6K');
Mdata20K = h5read('magnetization_data_calculation_mpms_data_AB_plane.h5', '/Mdata20K');
CESMHdata = h5read('magnetization_data_calculation_mpms_data_AB_plane.h5', '/CESMHdata');

% Plot the data
figure; grid on; box on; 
hold on;

% MPMS data
plot(Mdata20K(:, 1), Mdata20K(:, 2) / 1.35, 'o', 'DisplayName', '20K MPMS data');
plot(Mdata6K(:, 1), Mdata6K(:, 2) / 1.35, 'o', 'DisplayName', '6K MPMS data');
plot(Mdata2K(:, 1), Mdata2K(:, 2) / 1.35, 'o', 'DisplayName', '2K MPMS data');

% Allen's paper data
plot(CESMHdata(:, 1) / 1e4, CESMHdata(:, 2), 'b.', 'DisplayName', 'From Allen''s paper');

% MFT and calculated magnetization
plot(MFTField_mpms, MFT2K, '-', 'DisplayName', 'MFT 2K');
plot(MFTField_mpms, MFT6K, '-', 'DisplayName', 'MFT 6K');
plot(MFTField_mpms, MFT20K, '-', 'DisplayName', 'MFT 20K');
plot(magF_mpms, magnetization2K, '--', 'DisplayName', '2K, no MFT');
plot(magF_mpms, magnetization6K, '--', 'DisplayName', '6K, no MFT');
plot(magF_mpms, magnetization20K, '--', 'DisplayName', '20K, no MFT');

% Customize the plot
xlim([0, 8]);
title('AB plane magnetization');
xlabel('Field (T)');
ylabel('Magnetization (\mu_B/Er)');
legend show;
hold off;

%% test out the susceptibility issue
% Load data from HDF5 file
h5file = 'test_susceptibility.h5';

% Read datasets

myinv01T = h5read(h5file, '/myinv01T');
neutroninv01T = h5read(h5file, '/neutroninv01T');
mysus01T = h5read(h5file, '/mysus01T');
neutronsus01T = h5read(h5file, '/neutronsus01T');
temps = h5read(h5file, '/temps');
mvsHinvchi = h5read(h5file, '/mvsHinvchi');
mvsHchi = h5read(h5file, '/mvsHchi');
mvsHtemps = h5read(h5file, '/mvsHtemps');
% Recreate the figure
figure; grid on; box on; 
hold on;

% Plot manually extracted data
plot(mvsHtemps, mvsHinvchi, 'o', 'DisplayName', 'inv chi extracted manually from Mvs H');

% Plot the calculated data
plot(temps, mysus01T, '--', 'DisplayName', 'Raman B params MFT 0.1T');
plot(temps, neutronsus01T, '-.', 'DisplayName', 'neutrons B params MFT 0.1T');

% Formatting
title('Calculated MFT susceptibility at 0.1T');
xlabel('Temperature (K)');
ylabel('1/\chi');
legend('show');

hold off;

%% mosre sus diagnosis
h5file = 'sus_diagnosis.h5';

% Read datasets
mag = h5read(h5file, '/mag');
field = h5read(h5file, '/field');
m = h5read(h5file, '/m');

figure; grid on; box on;  
plot(field, m, 'DisplayName',   'mft magnetism')

%% fine spaced M vs H
h5file = 'M_vs_H_temperature_dependence_fine_spacing.h5';

% Read datasets
temps = h5read(h5file, '/temps');
field = h5read(h5file, '/f');
tempMag = h5read(h5file, '/tempMag');

figure; grid on; box on; 
hold on;
cmap = jet(length(temps)); % Colormap for consistency

for i = 1:length(temps)
    plot(field, tempMag(:, i), 'DisplayName', sprintf('%.3f K', temps(i)), 'Color', cmap(i, :));
end

% legend show;
% xlim([0 9]);
title({'C axis magnetization MFT', ...
       'calculated from Raman fit B params', ...
       'test B60 = <INSERT VALUE>'});
xlabel('Field (T)');
ylabel('Magnetization');
hold off;

%% M vs H for ab plane
h5file = 'M_vs_H_temperature_dependence_AB_plane.h5';

% Read datasets
temps = h5read(h5file, '/temps');
field = h5read(h5file, '/MFTField');
tempMag = h5read(h5file, '/tempMag');

figure; grid on; box on; 
hold on;
cmap = jet(length(temps)); % Colormap for consistency

for i = 1:length(temps)
    plot(field, tempMag(:, i), 'DisplayName', sprintf('%.3f K', temps(i)), 'Color', cmap(i, :));
end

% legend show;
% xlim([0 9]);
title({'B axis magnetization MFT', ...
       'Above 1K due to convergence issue'});
xlabel('Field (T)');
ylabel('Magnetization');
hold off;

%% dmdh AB plane
h5file = 'dMdH_temperature_dependence_AB_plane.h5'; 

% Read datasets
temps = h5read(h5file, '/temps');
field = h5read(h5file, '/MFTField');
dMdH = h5read(h5file, '/dMdH');

figure; grid on; box on; 
hold on;
cmap = jet(length(temps)); % Colormap for consistency

for i = 1:length(temps)
    plot(field, dMdH(:, i), 'DisplayName', sprintf('%.3f K', temps(i)), 'Color', cmap(i, :));
end

% legend show;
% xlim([0 9]);
title({'B axis dMdH MFT', ...
       'Above 1K due to convergence issue'});
xlabel('Field (T)');
ylabel('Magnetization');
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot the new mft code
% Load data from HDF5 file
filename = 'hopes_MFT_c_axis_calculation.h5';

% Read datasets
MFTField = h5read(filename, '/MFTField');
magnetization2K = h5read(filename, '/magnetization2K');
magnetization6K = h5read(filename, '/magnetization6K');
magnetization20K = h5read(filename, '/magnetization20K');
allenMFTCaxis = h5read(filename, '/allenMFTCaxis');
CESMHdata = h5read(filename, '/CESMHdata');
tempMagC = h5read(filename, '/tempMagC');
labels = h5read(filename, '/labels');
H = h5read(filename, '/calculationH'); 

% Read attributes
B20 = h5readatt(filename, '/', 'B20');
B40 = h5readatt(filename, '/', 'B40');
B43 = h5readatt(filename, '/', 'B43');
B60 = h5readatt(filename, '/', 'B60');
B63 = h5readatt(filename, '/', 'B63');
B66 = h5readatt(filename, '/', 'B66');
Jz = h5readatt(filename, '/', 'Jz');
JzAllen = h5readatt(filename, '/', 'JzAllen');

% Prepare colors for plotting
nColors = size(tempMagC, 2); % Number of temperature magnetizations
colors = parula(nColors); % MATLAB colormap

% Create figure
figure;
hold on;
grid on;

% % Plot MPMS data
plot(magnetization20K(:, 1), magnetization20K(:, 2) / 1.35, 'o', 'DisplayName', '20K MPMS data');
plot(magnetization6K(:, 1), magnetization6K(:, 2) / 1.35, 'o', 'DisplayName', '6K MPMS data');
plot(magnetization2K(:, 1), magnetization2K(:, 2) / 1.35, 'o', 'DisplayName', '2K MPMS data');

% Plot CES MH data
plot(CESMHdata( :, 7) / 1e4, CESMHdata( :, 8), 'b.', 'DisplayName', 'from Allens paper');

% Plot Allen's 2K MFT calculation
plot(MFTField, allenMFTCaxis, 'b--', 'DisplayName', 'Allens 2K MFT calculation');

% Plot temperature-dependent magnetizations
for i = 1:nColors
    plot(H, tempMagC(:, i), '-', 'Color', colors(i, :), 'DisplayName', strcat('My params my code', labels(i)));
end

% Set plot limits
xlim([0, 7]);
ylim([0, 8]);

% Add legend
legend('Location', 'best');

% Add title and labels
title(sprintf('C magnetization\n B20 = %.7f B40 = %.7f B43 = %.7f\n B60 = %.7f B63 = %.7f B66 = %.7f\n Jz = %.7f JzAllen = %.7f', ...
    B20, B40, B43, B60, B63, B66, Jz, JzAllen));
xlabel('Field (T)');
ylabel('Magnetization (\mu_B/Er)');

% now load my params allens code
filename = 'M_vs_H_temperature_dependence_myParams_allensCode.h5';

% Read datasets
temps = h5read(filename, '/temps');
tempMag = h5read(filename, '/tempMag');
MFTField = h5read(filename, '/MFTField');

nColors = size(tempMag, 2); % Number of temperature magnetizations
colors = parula(nColors); 

for i = 1:nColors
    plot(MFTField, tempMag(:, i), '-', 'Color', colors(i, :), 'DisplayName', ['my B params allens code', num2str(temps(i)), 'K']);
end

% now load my params allens code
filename = 'M_vs_H_temperature_dependence_allens_code_allens_params.h5';

% Read datasets
temps = h5read(filename, '/temps');
tempMag = h5read(filename, '/tempMagAllen');
MFTField = h5read(filename, '/MFTField');

nColors = size(tempMag, 2); % Number of temperature magnetizations
colors = parula(nColors); 

for i = 1:nColors
    plot(MFTField,tempMag(:, i), '-', 'Color', colors(i, :), 'DisplayName', ['Allens B params allens code', num2str(temps(i)), 'K']);
end
%% plot the new mft code
% Load data from HDF5 file
filename = 'hopes_MFT_b_axis_calculation.h5';

% Read datasets
MFTField = h5read(filename, '/MFTField');
magnetization2K = h5read(filename, '/magnetization2K');
magnetization6K = h5read(filename, '/magnetization6K');
magnetization20K = h5read(filename, '/magnetization20K');
allenMFTBaxis = h5read(filename, '/allenMFTBaxis');
CESMHdata = h5read(filename, '/CESMHdata');
tempMagB = h5read(filename, '/tempMagB');
labels = h5read(filename, '/labels');
H = h5read(filename, '/calculationH'); 

% Read attributes
B20 = h5readatt(filename, '/', 'B20');
B40 = h5readatt(filename, '/', 'B40');
B43 = h5readatt(filename, '/', 'B43');
B60 = h5readatt(filename, '/', 'B60');
B63 = h5readatt(filename, '/', 'B63');
B66 = h5readatt(filename, '/', 'B66');
Jperp = h5readatt(filename, '/', 'Jperp');
JperpAllen = h5readatt(filename, '/', 'JperpAllen');

% Prepare colors for plotting
nColors = size(tempMagB, 2); % Number of temperature magnetizations
colors = parula(nColors); % MATLAB colormap

% Create figure
figure;
hold on;
grid on;

% Plot MPMS data
plot(magnetization20K(:, 1), magnetization20K(:, 2) / 1.35, 'o', 'DisplayName', '20K MPMS data');
plot(magnetization6K(:, 1), magnetization6K(:, 2) / 1.35, 'o', 'DisplayName', '6K MPMS data');
plot(magnetization2K(:, 1), magnetization2K(:, 2) / 1.35, 'o', 'DisplayName', '2K MPMS data');

% Plot CES MH data
plot(CESMHdata( :, 1) / 1e4, CESMHdata( :, 2), 'b.', 'DisplayName', 'from Allens paper');

% Plot Allen's 2K MFT calculation
plot(MFTField, allenMFTBaxis, 'b--', 'DisplayName', 'Allens 2K MFT calculation');

% Plot temperature-dependent magnetizations
for i = 1:nColors
    plot(H, tempMagB(:, i), '-', 'Color', colors(i, :), 'DisplayName', strcat('my params my code', labels(i)));
end

% Set plot limits
xlim([0, 7]);
ylim([0, 8]);

% Add legend
legend('Location', 'best');

% Add title and labels
title(sprintf('B magnetization\n B20 = %.7f B40 = %.7f B43 = %.7f\n B60 = %.7f B63 = %.7f B66 = %.7f\n Jperp = %.7f JperpAllen = %.7f', ...
    B20, B40, B43, B60, B63, B66, Jperp, JperpAllen));
xlabel('H [T]');
ylabel('Magnetization [\mu_B/Er]');

% now load my params allens code
filename = 'M_vs_H_temperature_dependence_AB_plane_myParams_Allens_code.h5';

% Read datasets
temps = h5read(filename, '/temps');
tempMag = h5read(filename, '/tempMag');
MFTField = h5read(filename, '/MFTField');

nColors = size(tempMag, 2); % Number of temperature magnetizations
colors = parula(nColors); 

for i = 1:nColors
    plot(MFTField, tempMag(:, i), '-', 'Color', colors(i, :), 'DisplayName', ['my B params allens code', num2str(temps(i)), 'K']);
end

% now load allens params allens code
filename = 'M_vs_H_temperature_dependence_AB_plane_allenParams_Allens_code.h5';

% Read datasets
temps = h5read(filename, '/temps');
tempMag = h5read(filename, '/tempMag');
MFTField = h5read(filename, '/MFTField');

nColors = size(tempMag, 2); % Number of temperature magnetizations
colors = parula(nColors); 

for i = 1:nColors
    plot(MFTField, tempMag(:, i), '-', 'Color', colors(i, :), 'DisplayName', ['Allens B params allens code', num2str(temps(i)), 'K']);
end

%% looking at J because the c axis 2K magnetization is so bad
% Load data from the HDF5 file
filename = 'checking_J_B60_3079.h5';

% Read datasets
testJz = h5read(filename, '/testJz');
mftArr = h5read(filename, '/mftArr');
MFTField = h5read(filename, '/MFTField');
magF_mpms = h5read(filename, '/magF_mpms');
magnetization2K = h5read(filename, '/magnetization2K');
MData2K = h5read(filename, '/MData2K');

% Read attributes
B60 = h5readatt(filename, '/', 'B60');

% Create figure
figure;
hold on;
grid on;

% Plot molecular field theory calculations
for i = 1:length(testJz)
    plot(MFTField, mftArr(:, i), '--', 'DisplayName', ['Jz = ' num2str(testJz(i))]);
end

% Plot 2K magnetization data
plot(magF_mpms, magnetization2K, '-', 'DisplayName', 'No MFT calculation');
plot(MData2K(:, 1), MData2K(:, 2) / 1.35, 'o', 'DisplayName', '2K MPMS data');


% Set plot limits
xlim([0, 7]);
ylim([0, 8]);

% Add legend
legend('Location', 'best');

% Add title and labels
title(sprintf(['C magnetization\n B60 = %.5f\n Jz Test = [' repmat('%.3e ', 1, length(testJz)) ']'], ...
    B60, testJz));
xlabel('Field (T)');
ylabel('Magnetization (\mu_B/Er)');

hold off;

%% testing Jz with multiple B60
% Define parameters
B60Arr = [3.25e-6, 3.2e-6, 3.154e-6, 3.1e-6, 3.08e-6, 3.06e-6, 3.05e-6, 3.03e-6];
numPlots = length(B60Arr);
fnames = {'dmdh_B60_test325e-06.h5','dmdh_B60_test32e-06.h5','dmdh_B60_test3154e-06.h5','dmdh_B60_test31e-06.h5', 'dmdh_B60_test308e-06.h5', 'dmdh_B60_test306e-06.h5', 'dmdh_B60_test305e-06.h5', 'dmdh_B60_test303e-06.h5'};
% Initialize figure
figure;
tiledlayout(ceil(numPlots / 2), 2, 'TileSpacing', 'Compact');

% Loop through each B60 value
for i = 1:numPlots
    % Construct the filename for the HDF5 file
    B60 = B60Arr(i);
    fname = fnames{i}; 
    
    % Read data from HDF5 file
    testJz = h5read(fname, '/testJz');
    H = h5read(fname, '/H');
    dmdhArr = h5read(fname, '/dmdhArr');
    B60_read = h5readatt(fname, '/', 'B60');
    
    % Ensure the B60 value matches
    if abs(B60 - B60_read) > 1e-10
        error('Mismatch in B60 values for file: %s', fname);
    end
    
    % Create subplot
    nexttile;
    hold on;
    grid on;
    
    % Plot dM/dH for each test Jz value
    for j = 1:length(testJz)
        plot(H, dmdhArr(:, j), 'DisplayName', sprintf('Jz = %.3e', testJz(j)));
    end
    
    % Add vertical line at 5.4
    xline(5.4, '--', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1, 'DisplayName', 'H = 5.4 T');
    
    % Add labels, legend, and title
    xlabel('Field (T)');
    ylabel('dM/dH');
    legend('Location', 'best', 'FontSize', 8);
    title(sprintf('dM/dH for B60 = %.6e', B60), 'FontSize', 10);
    
    hold off;
end

% Add overall title
sgtitle('dM/dH vs Field (T) for various B60 values', 'FontSize', 12);



%% testing Jz with multiple B60
% Define parameters
fnames = {'dmdh_B60_test325e-06.h5','dmdh_B60_test32e-06.h5','dmdh_B60_test3154e-06.h5','dmdh_B60_test31e-06.h5', 'dmdh_B60_test308e-06.h5', 'dmdh_B60_test306e-06.h5', 'dmdh_B60_test305e-06.h5', 'dmdh_B60_test303e-06.h5'};
% Initialize figure
figure;
tiledlayout(4, 2, 'TileSpacing', 'Compact');
testJz = {}; 
dmdhArr = {}; 
B60 = []; 
%let's start by looping through the files and getting the data
for i = 1:length(fnames)
    fname = fnames{i}; 

    % read data
    testJz{end+1} = h5read(fname, '/testJz');
    dmdhArr{end+1} = h5read(fname, '/dmdhArr');
    B60(end+1) = h5readatt(fname, '/', 'B60');
    H = h5read(fname, '/H');

end

for i = 1:length(testJz{1,1})
    % Construct the filename for the HDF5 file
    
    % Create subplot
    nexttile;
    hold on;
    grid on;
    
    % Plot dM/dH for each B60 value
    for j = 1:length(B60)
        plot(H, dmdhArr{1,j}(:,i), 'DisplayName', sprintf('B60 = %.3e', B60(j)));
    end
    
    % Add vertical line at 5.4
    xline(5.4, '--', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1, 'DisplayName', 'H = 5.4 T');
    
    % Add labels, legend, and title
    xlabel('Field (T)');
    ylabel('dM/dH');
    % legend('Location', 'best', 'FontSize', 8);
    title(sprintf('dM/dH for J = %.6e', testJz{1,1}(i)), 'FontSize', 10);
    
    hold off;
end

% Add overall title
sgtitle('dM/dH vs Field (T) for various B60 values', 'FontSize', 12);
%% new the same for magnetization
% Define parameters
B60Arr = [2.75e-6];%[3.4e-6, 3.35e-6, 3.3e-6,3.25e-6, 3.2e-6, 3.154e-6, 3.1e-6, 3.08e-6, 3.06e-6, 3.05e-3, 3.03e-3, 3.0e-6, 2.9e-6, 2.8e-6];
numPlots = length(B60Arr);
fnames = {'m_B60_test34e-06.h5','m_B60_test335e-06.h5','m_B60_test33e-06.h5', 'm_B60_test325e-06.h5', 'm_B60_test32e-06.h5', 'm_B60_test3154e-06.h5'...
    'm_B60_test31e-06.h5', 'm_B60_test308e-06.h5','m_B60_test306e-06.h5'...
    'm_B60_test305e-06.h5', 'm_B60_test303e-06.h5', 'm_B60_test3e-06.h5', 'm_B60_test29e-06.h5', 'm_B60_test28e-06.h5'};
% Initialize figure
fnames = {'m_B60_test3154e-06.h5','m_B60_test31e-06.h5', 'm_B60_test308e-06.h5', 'm_B60_test303e-06.h5','m_B60_test29e-06.h5','m_B60_test275e-06.h5'};
figure;
tiledlayout(ceil(length(fnames) / 2), 2, 'TileSpacing', 'Compact');

% Loop through each B60 value
for i = 1:length(fnames)
    % Construct the filename for the HDF5 file
    fname = fnames{i};
    
    % Read data from HDF5 file
    testJz = h5read(fname, '/testJz');
    H = h5read(fname, '/H');
    mag = h5read(fname, '/mag');
    Mdata2K = h5read(fname, '/Mdata2K');
    B60 = h5readatt(fname, '/', 'B60');
    
    
    % Extract 2K MPMS data
    Mdata2K_H = Mdata2K(:, 1); % Field values
    Mdata2K_M = Mdata2K(:, 2) / 1.35; % Magnetization values, scaled
    
    % Create subplot
    nexttile; 
    hold on;
    grid on;
    
    % Plot 2K MPMS data
    plot(Mdata2K_H, Mdata2K_M, 'o', 'DisplayName', '2K MPMS data');
    
    % Plot M for each test Jz value
    for j = 1:length(testJz)
        plot(H, mag(:, j), 'DisplayName', sprintf('Jz = %.3e', testJz(j)));
    end
    
    % Add vertical line at 5.4
    % xline(5.4, '--', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1, 'DisplayName', 'H = 5.4 T');
    
    % Add labels, legend, and title
    xlabel('H (T)');
    ylabel('M (\mu_B/Er)');
    % legend('Location', 'best', 'FontSize', 8);
    title(sprintf('Magnetization for B60 = %.6e', B60), 'FontSize', 10);
    xlim([0 7])
    hold off;
end

% Add overall title
sgtitle('Magnetization vs Field (T) for various B60 values', 'FontSize', 12);
