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
subplot(2,3,1)
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
subplot(2,3,4);  grid on; box on; 
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

xlim([0, 6]);
ylim([0, 7]);

title('AB plane magnetization');
xlabel('Field (T)');
ylabel('Magnetization (\mu_B/Er)');
legend show;
hold off;


% %% dmdh with data

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
colors = cool(n); % wanted inferno, can't find it :(

subplot(1,3,2); grid on; box on; 
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
    text(9, offset+0.2, labels{i}, 'FontSize', 9);
end

title('dM/dH from SCM1 \n calculated dM/dH in dotted line');
ylabel('dM/dH (arb)');
xlabel('Field (T)');
hold off;

% %% make integrated scm1 data

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
colors = cool(n);

subplot(1,3,3); grid on; box on; 
hold on;
for i = [1, 4, n]
    % Sort and normalize x and y data
    x = xArrs{i};
    y = yArrs{i};
    [x, sortIdx] = sort(x);
    y = y(sortIdx) / max(y);
    intMag = integratedMag{i} / max(integratedMag{i}); % Normalize integratedMag
    
    % Plot data
    plot(x, intMag , 'DisplayName', labels{i}, 'Color', colors(i, :));
    plot(magF, tempMag(:, i) / max(tempMag(:,i)), '--', 'Color', colors(i, :));
end

% Add labels, title, and legend
title('Integrated chi(H) \n Numerically Integrated from SCM1 Data \n Calculated Curve in Dotted');
xlabel('Field (T)');
ylabel('Magnetization (arb)');
legend('show');
hold off;
