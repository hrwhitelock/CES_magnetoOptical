% final fig generation script 12/16/24

% let's start with making the paper figures
%% load my spectroscopy data
% filename = 'spectroscopy_calculation_hopes_params_2025Jan30_02.h5';
filename = 'spectroscopy_fgr_only_high_lines2025Aug07.h5';
% filename = 'spectroscopy_fgr_all_2025Aug07.h5';
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

%% load H||a test spectroscopy data
filename = 'spectroscopy_HparA_2025jul17.h5';

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
    my_spec_data.(datasetName) = datasetData;
end
%% make colormap
yellowMap = [linspace(0, 1, 256)', linspace(0, 1, 256)', zeros(256, 1)];
colormap(yellowMap);

n = 256;  % Number of color levels

base_cmap = [
    1.0000, 1.0000, 0.8980
    1.0000, 0.9686, 0.7373
    0.9961, 0.8902, 0.5686
    0.9961, 0.8000, 0.4353
    0.9922, 0.6824, 0.3804
    0.9451, 0.5176, 0.3216
    0.8510, 0.3725, 0.2588
    0.6510, 0.2118, 0.1569
    0.4980, 0.1529, 0.1176
];

% Interpolate to desired number of colors
x = linspace(0, 1, size(base_cmap,1));
xq = linspace(0, 1, n);

oranges = interp1(x, base_cmap, xq, 'linear');

base_colors = [...
    0.556863, 0.003922, 0.321569;
    0.772549, 0.105882, 0.490196;
    0.870588, 0.466667, 0.682353;
    0.945098, 0.713725, 0.854902;
    0.992157, 0.878431, 0.937255;
    1.000000, 1.000000, 0.850980;
    0.901961, 0.960784, 0.596078;
    0.670588, 0.866667, 0.643137;
    0.400000, 0.741176, 0.674510;
    0.196078, 0.533333, 0.741176;
    0.192157, 0.211765, 0.584314];

% Interpolate to desired number of colors
x = linspace(0, 1, size(base_colors, 1));
xi = linspace(0, 1, n);

r = interp1(x, base_colors(:,1), xi, 'pchip');
g = interp1(x, base_colors(:,2), xi, 'pchip');
b = interp1(x, base_colors(:,3), xi, 'pchip');

PiYG = [r(:), g(:), b(:)];

base_cmap = [ ...
    0.2298057, 0.29871797, 0.75368315;
    0.26623388, 0.35309472, 0.80146664;
    0.30386891, 0.4065353 , 0.84495867;
    0.34280448, 0.45875775, 0.8837259 ;
    0.38301334, 0.50941904, 0.91738754;
    0.42436961, 0.55814842, 0.94561959;
    0.46666708, 0.60456257, 0.96815491;
    0.5096352 , 0.64828077, 0.98478814;
    0.55295333, 0.68892991, 0.9953756 ;
    0.59626216, 0.72614942, 0.99983624;
    0.63917622, 0.75959947, 0.99815161;
    0.68129128, 0.78896471, 0.9903632 ;
    0.72219329, 0.81395274, 0.97657471;
    0.76146477, 0.83430265, 0.95694527;
    0.79869163, 0.84978614, 0.93168855;
    0.8334661 , 0.86020798, 0.9010687 ;
    0.86539519, 0.86541021, 0.86539519;
    0.89609127, 0.84893747, 0.82088012;
    0.92290314, 0.82738471, 0.7745082 ;
    0.94647675, 0.80092744, 0.72673615;
    0.96660839, 0.76976775, 0.67800728;
    0.98324389, 0.73421376, 0.62875189;
    0.99654221, 0.69473848, 0.57937515;
    1.        , 0.65197766, 0.53026376;
    0.99581563, 0.60659675, 0.48177591;
    0.98386829, 0.55926553, 0.43424371;
    0.96412339, 0.51052391, 0.38797432;
    0.93667127, 0.46075608, 0.34326326;
    0.90170242, 0.41019373, 0.30031808;
    0.85951726, 0.35900052, 0.25930125;
    0.81046763, 0.30728987, 0.22033408;
    0.75492178, 0.25514502, 0.18352669;
    0.69326246, 0.20264226, 0.14893563;
    0.62589536, 0.14983829, 0.11659704];

% Interpolate to n colors
x = linspace(0, 1, size(base_cmap, 1));
xi = linspace(0, 1, n);

r = interp1(x, base_cmap(:,1), xi, 'pchip');
g = interp1(x, base_cmap(:,2), xi, 'pchip');
b = interp1(x, base_cmap(:,3), xi, 'pchip');

coolwarm = [r(:), g(:), b(:)];

% Anchor points for the BWR colormap from Matplotlib
base_cmap = [ ...
    0.0, 0.0, 1.0;    % blue
    1.0, 1.0, 1.0;    % white
    1.0, 0.0, 0.0];   % red

% Define interpolation scale
x = linspace(0, 1, size(base_cmap, 1));
xi = linspace(0, 1, n);

% Interpolate in RGB space
r = interp1(x, base_cmap(:,1), xi, 'linear');
g = interp1(x, base_cmap(:,2), xi, 'linear');
b = interp1(x, base_cmap(:,3), xi, 'linear');

% Combine into colormap
bwr = [r(:), g(:), b(:)];

base_cmap = [ ...
    0.0, 0.0, 0.3;   % dark blue
    0.0, 0.0, 1.0;   % blue
    1.0, 1.0, 1.0;   % white
    1.0, 0.0, 0.0;   % red
    0.5, 0.0, 0.0];  % dark red

% Corresponding positions for interpolation
x = linspace(0, 1, size(base_cmap, 1));
xi = linspace(0, 1, n);

% Interpolate each channel
r = interp1(x, base_cmap(:,1), xi, 'linear');
g = interp1(x, base_cmap(:,2), xi, 'linear');
b = interp1(x, base_cmap(:,3), xi, 'linear');

% Combine into colormap
seismic = [r(:), g(:), b(:)];

base_cmap = [ ...
    158,   1,  66;
    213,  62,  79;
    244, 109,  67;
    253, 174,  97;
    254, 224, 139;
    255, 255, 191;
    230, 245, 152;
    171, 221, 164;
    102, 194, 165;
     50, 136, 189;
     94,  79, 162] / 255;

% Positions for interpolation
x = linspace(0, 1, size(base_cmap, 1));
xi = linspace(0, 1, n);

% Interpolate RGB channels
r = interp1(x, base_cmap(:,1), xi, 'linear');
g = interp1(x, base_cmap(:,2), xi, 'linear');
b = interp1(x, base_cmap(:,3), xi, 'linear');

% Combine into colormap
spectral = [r(:), g(:), b(:)];
%% make spec without thermal lines: 
% Load data (you probably already have this loaded)
field_vals = my_spec_data.calc_field;   % e.g., size [Nfield × 1]
linesC = my_spec_data.linesB;           % e.g., size [Nfield × Nlines]

% Define energy axis
Emin = 15; Emax = 300;  % cm^-1 range
dE = 0.05;
E = Emin:dE:Emax;      % energy grid
Nfield = length(field_vals);
Nenergy = length(E);
Z = zeros(Nenergy, Nfield);

% Gaussian width in cm^-1
sigma = .8;

% Loop through each spectral line
for i = 2:121
    for f = 1:Nfield
        E0 = linesC(f, i)*8.022;  % Center of the Gaussian
        if isnan(E0), continue; end  % Skip if missing
        amp = my_spec_data.ampB(f,i);
        % Gaussian over energy axis
        G = amp.*exp(-((E - E0).^2) / (2*sigma^2));

        % Add to spectrum at field f
        Z(:, f) = Z(:, f) + G';
    end
end
% now ZF subtract
Z = Z(:,2:end)-Z(:,1); 
% Normalize (optional)
% Z = Z / max(max(Z(:)));

avg = sum(Z,2);
Z = Z-avg/99;
Z = Z / max(max(Z(:)));
% Plot
figure;
imagesc(field_vals(2:end), E, Z);
axis xy;
xlabel('Field [T]');
ylabel('Energy [cm^{-1}]');
title('Simulated 2D Spectrum from Lines');
colormap(coolwarm);
colorbar;

%% 
fig = figure;
ax1 = subplot(2,3,1);
title('B-axis IR')
hold on; 
pcolor(my_spec_data.IR_B_field,my_spec_data.IR_B_wavenums,my_spec_data.IR_dataB')
% axis xy;
shading flat
for i = 2:16%size(arrC, 2)
    plot(my_spec_data.calc_field, my_spec_data.linesB(:,i),  'r-', 'LineWidth', 1);
end
% % for i = 17:50%size(arrC, 2)
%     plot(my_spec_data.calc_field, my_spec_data.linesB(:,i),  'r--', 'LineWidth', 1);
% end
% set(ax1,'Xticklabel',[])
ylabel('Energy [cm{^-1}]')
clim([0 1])
ylim([0 100])
xlim([0 17.5])
colormap(ax1, coolwarm)
% colormap(seismic)
set(gca,'Layer','top')



%% c axis IR
figure; 
hold on; 
ax2 = subplot(2,3,2);hold on; box on; 
pcolor(my_spec_data.IR_C_field, my_spec_data.IR_C_wavenums, my_spec_data.IR_dataC')
shading flat
for i = 2:54 %size(my_spec_data.linesC, 2)
    if i <= 16
        plot(my_spec_data.calc_field, my_spec_data.linesC(:,i),  'r-', 'LineWidth',1, "DisplayName",['0 to' num2str(i-1)]);
    end
%     if i>=17 && i<31
%         plot(my_spec_data.calc_field, my_spec_data.linesC(:,i),  'r-', 'LineWidth',1, "DisplayName",['1 to' num2str(i-15)]);
%     end
%     if i>=31 && i<43
%         plot(my_spec_data.calc_field, my_spec_data.linesC(:,i),  'r-', 'LineWidth',1, "DisplayName",['2 to' num2str(i-15-14+1)])
%     end
%     if i>=43 && i<=54
%         plot(my_spec_data.calc_field, my_spec_data.linesC(:,i),  'r-', 'LineWidth',1, "DisplayName",['3 to' num2str(i-15-14-13+3)])
%     end
end
title('C-axis IR')
% ylabel('Energy [cm{^-1}]')
colormap(ax2, coolwarm)
clim([0 1])
ylim([0 100])
xlim([0 17.5])
set(gca,'Layer','top')
% 
% ax4 = subplot(2,2,4);
% hold on; box on; 
% contourf(my_spec_data.calc_field(field_idx),my_spec_data.calc_wavenums(idx), my_spec_data.simulated_IR_C(idx, field_idx), 100, 'LineStyle', 'none');
% colormap(ax4, jet)
% for i = 2:16%size(arrC, 2)
%     plot(my_spec_data.calc_field(field_idx), my_spec_data.linesC(field_idx,i), 'r--', 'LineWidth', 1);
% end
% clim([0 1])
% ylim([0 100])
% xlim([0 17.5])
% ylabel({'Data', 'Energy [cm{^-1}]'})
xlabel('H(T)')
set(gca,'Layer','top')

% ax1.Position(2)=ax3.Position(2)+ax3.Position(4)+.002; 
% ax3.Position(1)=ax1.Position(1);
% ax2.Position(1) = ax1.Position(1) +ax1.Position(3) +.002;
% 
% ax2.Position(2)=ax1.Position(2);
% ax4.Position(1)=ax2.Position(1);

% linkaxes([ax1, ax2, ax3, ax4], 'xy')
% set(ax2,'Yticklabel',[]) 
% set(ax1,'Xticklabel',[])
% set(ax2,'Xticklabel',[])
% set(ax4,'Yticklabel',[])
% 

% c axis raman

% fig = figure;
ax3 = subplot(2,3,3);box on; hold on; 
pcolor(my_spec_data.raman_field, my_spec_data.raman_wavenums, my_spec_data.ramanData'); 
shading flat; 
for i = 2:16%size(arrC, 2)
    plot(my_spec_data.calc_field, my_spec_data.linesC(:,i), 'r', 'LineWidth', 1);
end
% for i = 17:50%size(arrC, 2)
%     plot(my_spec_data.calc_field,my_spec_data.linesC(:,i),  'r--', 'LineWidth', 1);
% % end
% ylim([0 100])
% xlim([0 14])
colormap(ax3, coolwarm)
title ('C-axis Raman')
% ylabel('Energy [cm^{-1}]')
% 
% ax2 = subplot(2,1,2);
% hold on; box on; 
% contourf(my_spec_data.calc_field(field_idx),my_spec_data.calc_wavenums(idx), my_spec_data.simulated_raman(idx, field_idx), 100, 'LineStyle', 'none');
% colormap(ax2, cm)
% % let's add some lines
% for i = 2:16%size(arrC, 2)
%     plot(my_spec_data.calc_field(field_idx), my_spec_data.linesC(field_idx,i), 'r--', 'LineWidth', 1);
% end
% clim([0 1])
% xlim([0 14])
% ylim([0 100])
% ylabel('Energy [cm{^-1}]')
% xlabel('H(T)')
% 
% ax1.Position(2)=ax2.Position(2)+ax2.Position(4)+.02;
% ax2.Position(1)=ax1.Position(1);
linkaxes([ax1, ax2, ax3], 'xy')
% set(ax3,'Yticklabel',[])
clim([0 1])
ylim([0 100])
xlim([0 17.5])

ax4 = subplot(2,3,4);
% title('B-axis IR')
hold on; 
% pcolor(my_spec_data.IR_B_field,my_spec_data.IR_B_wavenums,my_spec_data.IR_dataB')
% % axis xy;
shading flat
for i = 2:16%size(arrC, 2)
    plot(my_spec_data.calc_field, my_spec_data.linesB(:,i),  'r-', 'LineWidth', 1);
end
% % for i = 17:50%size(arrC, 2)
%     plot(my_spec_data.calc_field, my_spec_data.linesB(:,i),  'r--', 'LineWidth', 1);
% end
% set(ax1,'Xticklabel',[])
ylabel('Energy [cm{^-1}]')
% clim([0 1])
% ylim([0 100])
% xlim([0 17.5])
% colormap(ax1, coolwarm)
% colormap(seismic)
set(gca,'Layer','top')

ax5 = subplot(2,3,5);box on; hold on; 
shading flat; 
for i = 2:16%size(arrC, 2)
    plot(my_spec_data.calc_field, my_spec_data.linesC(:,i), 'r', 'LineWidth', 1);
end

ax6 = subplot(2,3,6);box on; hold on; 
for i = 2:16%size(arrC, 2)
    plot(my_spec_data.calc_field, my_spec_data.linesC(:,i), 'r', 'LineWidth', 1);
end

linkaxes([ax1, ax2, ax3], ax4, ax5, ax6, 'xy')
% set(ax3,'Yticklabel',[])
% clim([0 1])
ylim([0 100])
xlim([0 17.5])


% ax2.Position(1) = ax1.Position(1)+ax1.Position(3) +.01;
% ax3.Position(1) = ax2.Position(1)+ax2.Position(3) +.01;

%% now make fig 3
% load my magnetic data
filename = 'mag_calc_mft_2025May12.h5';

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
plot(my_mag_data.magnetization20K(:,1), my_mag_data.magnetization20K(:,2)/1.37, 'o', 'DisplayName', '20K MPMS data scaled by 1.37');
plot(my_mag_data.magnetization6K(:,1), my_mag_data.magnetization6K(:,2)/1.37, 'o', 'DisplayName', '6K MPMS data scaled by 1.37');
plot(my_mag_data.magnetization2K(:,1), my_mag_data.magnetization2K(:,2)/1.37, 'o', 'DisplayName', '2K MPMS data scaled by 1.37');
plot(my_mag_data.CESMHdata(:,7)./1e4, my_mag_data.CESMHdata(:,8), 'b.', 'DisplayName', 'From Allens paper');
% 
% plot(my_mag_data.magnetization20K(:,1), my_mag_data.magnetization20K(:,2), 'x', 'DisplayName', '20K MPMS data');
% plot(my_mag_data.magnetization6K(:,1), my_mag_data.magnetization6K(:,2), 'x', 'DisplayName', '6K MPMS data');
% plot(my_mag_data.magnetization2K(:,1), my_mag_data.magnetization2K(:,2), 'x', 'DisplayName', '2K MPMS data');

% Plot MFT data
% H = horzcat(linspace(0,1,50), linspace(1.01,15, 150));
H = horzcat(linspace(0,1,100), linspace(1.01, 10, 100));
plot(my_mag_data.H, my_mag_data.tempMagC(:,11), '-', 'DisplayName', 'MFT 2K');
plot(my_mag_data.H, my_mag_data.tempMagC(:,12), '-', 'DisplayName', 'MFT 6K');
plot(my_mag_data.H, my_mag_data.tempMagC(:,13), '-', 'DisplayName', 'MFT 20K');

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
subplot(1,2,1);  grid on; box on; 
hold on;

% MPMS data
plot(my_mag_data.magnetizationAb20K(:,1), my_mag_data.magnetizationAb20K(:,2)/1.37, 'o', 'DisplayName', '20K MPMS data scaled by 1.37');
plot(my_mag_data.magnetizationAB6K(:,1), my_mag_data.magnetizationAB6K(:,2)/1.37, 'o', 'DisplayName', '6K MPMS data scaled by 1.37');
plot(my_mag_data.magnetizationAB2K(:,1), my_mag_data.magnetizationAB2K(:,2)/1.37, 'o', 'DisplayName', '2K MPMS data scaled by 1.37');


plot(my_mag_data.magnetizationAb20K(:,1), my_mag_data.magnetizationAb20K(:,2), 'x', 'DisplayName', '20K MPMS data');
plot(my_mag_data.magnetizationAB6K(:,1), my_mag_data.magnetizationAB6K(:,2), 'x', 'DisplayName', '6K MPMS data');
plot(my_mag_data.magnetizationAB2K(:,1), my_mag_data.magnetizationAB2K(:,2), 'x', 'DisplayName', '2K MPMS data');

% Allen's paper data
plot(my_mag_data.CESMHdata(:, 1) / 1e4, my_mag_data.CESMHdata(:, 2), 'b.', 'DisplayName', 'From Allen''s paper');

% MFT and calculated magnetization
plot(my_mag_data.H, my_mag_data.tempMagB(:,11), '-', 'DisplayName', 'MFT 2K');
plot(my_mag_data.H, my_mag_data.tempMagB(:,12), '-', 'DisplayName', 'MFT 6K');
plot(my_mag_data.H, my_mag_data.tempMagB(:,13), '-', 'DisplayName', 'MFT 20K');

xlim([0, 6]);
ylim([0, 7]);

title('AB plane magnetization');
xlabel('Field (T)');
ylabel('Magnetization (\mu_B/Er)');
legend show;
hold off;

%% dmdh with data
figure;
labels = my_mag_data.dmdhLabels;
n = length(labels);
colors = abyss(n);

H_vals = my_mag_data.H;
temps = my_mag_data.temps;

allenT = chi_from_allen.temperature;
allenChi = chi_from_allen.chi_prime;

hold on; box on; grid on;

for i = 1:n
    % === Load experimental data ===
    H_exp = my_mag_data.dmdhField{i}; 
    dmdh_exp = my_mag_data.dmdhData{i};

    % Sort and truncate to H >= 0
    [H_exp, sortIdx] = sort(H_exp);
    dmdh_exp = dmdh_exp(sortIdx);
    mask_exp = H_exp >= 0;
    H_exp = H_exp(mask_exp);
    dmdh_exp = dmdh_exp(mask_exp);

    % Normalize and scale
    area_exp = trapz(H_exp, dmdh_exp);
    dmdh_exp = dmdh_exp / area_exp;

    [~, t_idx] = min(abs(allenT - temps(i)));
    scale = allenChi(t_idx);
    dmdh_exp = dmdh_exp * scale;

    % Plot experimental
    plot(H_exp, dmdh_exp, '-', ...
        'Color', colors(i, :), ...
        'DisplayName', [labels{i}, ' (exp)']);

    % === Load and truncate calculated data to ≤12 T ===
    H_calc = H_vals;
    dmdh_calc = my_mag_data.dmdhC(:, i);
    mask_calc = (H_calc >= 0) & (H_calc <= 12);
    H_calc = H_calc(mask_calc);
    dmdh_calc = dmdh_calc(mask_calc);

    % Normalize and scale (after truncation)
    area_calc = trapz(H_calc, dmdh_calc);
    dmdh_calc = dmdh_calc / area_calc;
    dmdh_calc = dmdh_calc * scale;

    % Plot calculated
    plot(H_calc, dmdh_calc, '--', ...
        'Color', colors(i, :), ...
        'DisplayName', [labels{i}, ' (calc ≤12T)']);
end

xlabel('Field (T)');
ylabel('dM/dH (normalized area × χ′)');
title('dM/dH Curves Truncated at 12T, Area-Normalized and χ′-Scaled');
legend('show');
hold off;

%% make integrated data
figure;
labels = my_mag_data.dmdhLabels;
n = length(labels);
colors = abyss(n);

H_vals = my_mag_data.H;
temps = my_mag_data.temps;

allenT = chi_from_allen.temperature;
allenChi = chi_from_allen.chi_prime;

hold on; box on; grid on;

for i = 1:n
    % === Experimental Integrated M(H) ===
    H_exp = my_mag_data.dmdhField{i};
    dmdh_exp = my_mag_data.dmdhData{i};

    % Sort and keep H ≥ 0
    [H_exp, sortIdx] = sort(H_exp);
    dmdh_exp = dmdh_exp(sortIdx);
    mask = H_exp >= 0;
    H_exp = H_exp(mask);
    dmdh_exp = dmdh_exp(mask);

    % Normalize dM/dH area to 1
    area_exp = trapz(H_exp, dmdh_exp);
    dmdh_exp = dmdh_exp / area_exp;

    % Integrate to get M(H)
    M_exp = cumtrapz(H_exp, dmdh_exp);

    % Scale by Allen χ'
    [~, t_idx] = min(abs(allenT - temps(i)));
    scale = allenChi(t_idx);
    M_exp = M_exp * scale;

    % Plot experimental integrated
    plot(H_exp, M_exp, '-', ...
        'DisplayName', [labels{i}, ' (exp)'], ...
        'Color', colors(i, :));

    % === Calculated Integrated M(H), Truncated to ≤ 12 T ===
    H_calc = H_vals;
    dmdh_calc = my_mag_data.dmdhC(:, i);

    % Keep H ∈ [0, 12] T
    mask = (H_calc >= 0) & (H_calc <= 12);
    H_calc = H_calc(mask);
    dmdh_calc = dmdh_calc(mask);

    % Normalize dM/dH area to 1
    area_calc = trapz(H_calc, dmdh_calc);
    dmdh_calc = dmdh_calc / area_calc;

    % Integrate to get M(H)
    M_calc = cumtrapz(H_calc, dmdh_calc);

    % Scale by Allen χ'
    M_calc = M_calc * scale;

    % Plot calculated integrated
    plot(H_calc, M_calc, '--', ...
        'DisplayName', [labels{i}, ' (calc ≤12T)'], ...
        'Color', colors(i, :));
end

xlabel('Field (T)');
ylabel('Magnetization (arb, scaled)');
title('Integrated M(H): Area-Normalized and χ′(T)-Scaled');
legend('show');
hold off;

%% try different normalization
figure;
labels = my_mag_data.dmdhLabels;
n = length(labels);
colors = abyss(n);

H_vals = my_mag_data.H;
temps = my_mag_data.temps;

allenT = chi_from_allen.temperature;
allenChi = chi_from_allen.chi_prime;

hold on; box on; grid on;

for i = 1:n
    % === Experimental dM/dH ===
    H_exp = my_mag_data.dmdhField{i};
    dmdh_exp = my_mag_data.dmdhData{i};

    % Truncate to H >= 0
    [H_exp, sortIdx] = sort(H_exp);
    dmdh_exp = dmdh_exp(sortIdx);
    mask_exp = H_exp >= 0;
    H_exp = H_exp(mask_exp);
    dmdh_exp = dmdh_exp(mask_exp);

    % Normalize area
    area_exp = trapz(H_exp, dmdh_exp);
    dmdh_exp_norm = dmdh_exp / area_exp;

    % Integrate
    M_exp = cumtrapz(H_exp, dmdh_exp_norm);

    % Scale by Allen chi
    [~, t_idx] = min(abs(allenT - temps(i)));
    scale_exp = allenChi(t_idx);
    M_exp = M_exp * scale_exp;

    % Plot experimental result
    plot(H_exp, M_exp, '-', ...
        'DisplayName', [labels{i}, ' (exp)'], ...
        'Color', colors(i, :));

    % === Calculated dM/dH ===
    H_calc = H_vals;
    dmdh_calc = my_mag_data.dmdhC(:, i);

    % Truncate to H in [0, 12]
    mask_calc = (H_calc >= 0) & (H_calc <= 12);
    H_calc = H_calc(mask_calc);
    dmdh_calc = dmdh_calc(mask_calc);

    %normalize area to 1

    % Find dM/dH at 12T in both calc and exp
    target_H = 12;
    % Interpolate experimental dM/dH at 12T
    dmdh_exp_at_12T = Interp1NonUnique(H_exp, dmdh_exp, target_H);
    dmdh_calc_at_12T = interp1(H_calc, dmdh_calc, target_H, 'linear', 'extrap');

    % Compute offset
    offset = dmdh_exp_at_12T - dmdh_calc_at_12T;
    dmdh_calc_shifted = dmdh_calc + offset;

    % Normalize area
    area_calc = trapz(H_calc, dmdh_calc_shifted);
    dmdh_calc_norm = dmdh_calc_shifted / area_calc;

    % Integrate
    M_calc = cumtrapz(H_calc, dmdh_calc_norm);
    M_calc = M_calc*scale_exp;

    % Plot calculated result
    plot(H_calc, M_calc, '--', ...
        'DisplayName', [labels{i}, ' (calc, offset+area norm)'], ...
        'Color', colors(i, :));
end

xlabel('Field (T)');
ylabel('Magnetization (arb)');
title('M(H): Experimental Scaled by Allen χ′, Theory Offset & Area-Normalized');
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
filename = 'sus_mft_0p48ueVJJz_2025May11.h5';

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
    [sortedTemp, idx] = sort(my_sus_data.temps); 
    sus = my_sus_data.susC(:,i); 
    sortedSus = sus(idx);
    plot(sortedTemp, sortedSus, 'DisplayName', [num2str(my_sus_data.fieldVals(i)), 'T, mean field'])
end
% add data
for i= 1:length(my_sus_data.data_sus_C)
    plot(my_sus_data.data_temps_C{i}, 1./(my_sus_data.data_sus_C{i}.*1.37), 'b.', 'DisplayName', [my_sus_data.clabels{i}, 'scaled by 1.37']); 

end

for i= 1:length(my_sus_data.data_sus_C)
    plot(my_sus_data.data_temps_C{i}, 1./(my_sus_data.data_sus_C{i}), 'bx', 'DisplayName', [my_sus_data.clabels{i}, ' no scaling']); 

end
title('Susceptibility, c-axis, my params')
legend(); 
xlabel('Temperature [K]'); 
ylabel('\chi')
%% ab plane
figure; hold on; box on; grid on; 
for i = 1: length(my_sus_data.susB(1,:))
    [sortedTemp, idx] = sort(my_sus_data.temps); 
    sus = my_sus_data.susB(:,i); 
    sortedSus = sus(idx);
    plot(sortedTemp, sortedSus, 'DisplayName', [num2str(my_sus_data.fieldVals(i)), 'T, mean field'])
end
% add data
for i= 1:length(my_sus_data.data_sus_AB)
    plot(my_sus_data.data_temps_AB{i}, 1./(my_sus_data.data_sus_AB{i}*1.35), 'b.', 'DisplayName', [my_sus_data.blabels{i}, 'scaled by 1.37']); 

end


for i= 1:length(my_sus_data.data_sus_AB)
    plot(my_sus_data.data_temps_AB{i}, 1./(my_sus_data.data_sus_AB{i}), 'b.', 'DisplayName', [my_sus_data.blabels{i}, 'no scaling']); 

end
title('Susceptibility, ab plane, my params')
legend(); 
xlabel('Temperature [K]'); 
ylabel('\chi')


%% low temp stuff I forgot about oops
figure; hold on; box on; grid on; 
% for i = 1: length(my_sus_data.lowTempSusC(1,:))
%     plot(linspace(0,1,100), my_sus_data.lowTempSusC(i,:), 'DisplayName', my_sus_data.lowTempLabels{i})
% end
% add data
for i= 1:length(my_sus_data.data_lowTemps)
    plot(my_sus_data.data_lowTemps{i, 1}, 1e5*my_sus_data.data_lowTempSus{i,1}/1.85, 'DisplayName', my_sus_data.lowTempLabels{i}); 

end
title()
legend(); 
xlabel('Temperature [K]'); 
ylabel('\chi (arb)')
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

%%
% load my magnetic data
filename = 'mag_calc_mft_2025May09.h5';

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

%% make temp dependence for M_c vs H
figure;
hold on;
grid on;


% Plot calculated magnetization curves
numCurves = size(my_mag_data.tempMagB, 2);
colors = jet(numCurves); % wanted inferno, can't find it :
for i = 1:numCurves
    plot(my_mag_data.H, my_mag_data.tempMagB(:, i), 'DisplayName', num2str(my_mag_data.temps(i)), 'color', colors(i, :));
end
legend()
xlabel('H[T]')
ylabel('M \mu_B/Er')
title('')

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

%% make fig 1
% Load data from HDF5
fileName = 'EvsH_mft_T10K_1015May12.h5';
ZFevals = h5read(fileName, '/ZFevals');
ABevals = h5read(fileName, '/ABevals');
Cevals = h5read(fileName, '/Cevals');
ABevals_nomft = h5read(fileName, '/ABevals_nomft');
Cevals_nomft = h5read(fileName, '/Cevals_nomft');
field = h5read(fileName, '/field');
B20 = h5readatt(fileName, '/', 'B20');
B40 = h5readatt(fileName, '/', 'B40');
B43 = h5readatt(fileName, '/', 'B43');
B60 = h5readatt(fileName, '/', 'B60');
B63 = h5readatt(fileName, '/', 'B63');
B66 = h5readatt(fileName, '/', 'B66');


% fileName = 'EvsH_noNorm_allens_params_2025Jan14.h5';
% ZFevals_allen = h5read(fileName, '/ZFevals');
% ABevals_allen = h5read(fileName, '/ABevals');
% Cevals_allen = h5read(fileName, '/Cevals');
% ABevals_allen_nomft = h5read(fileName, '/ABevals_nomft');
% Cevals_allen_nomft = h5read(fileName, '/Cevals_nomft');
% field = h5read(fileName, '/field');
% B20_allen = h5readatt(fileName, '/', 'B20');
% % Create the figure with subplots
figure;
hold on;
tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% Plot ZFevals (H = 0)
nexttile;hold on;
for i = 1:size(ZFevals)
    yline(ZFevals(i)+abs(ZFevals(1)), 'LineWidth', 1.2);
    % yline(ZFevals_allen(i)+abs(ZFevals_allen(1)), 'LineWidth', 1.2, 'LineStyle','--')
end
title('H = 0');
% xlabel('Field (T)');
ylabel('Energy');

% Plot ABevals (H || b)
nexttile;hold on;
hold on;
for i = 1:size(ABevals, 2)
    plot(field, ABevals(:, i)-ABevals(1, 1), 'color', 'black','LineWidth', 1.2, 'DisplayName', 'my params mft');
    % plot(field, ABevals_nomft(:, i)-ABevals_nomft(:, 1), 'color', 'red','LineStyle', '--','LineWidth', 1.2, 'DisplayName', 'my params no mft');
    % plot(field, ABevals_allen(:, i)-ABevals_allen(:, 1),'color', 'cyan', 'LineWidth', 1.2, 'LineStyle', ':', 'DisplayName', 'allen params mft');
    % plot(field, ABevals_allen_nomft(:, i)+abs(ZFevals_allen(1)),'color', 'blue', 'LineWidth', 1.2, 'LineStyle', '-.', 'DisplayName', 'allen params no mft');
end
title('H || b');
xlabel('Field (T)');

% Plot Cevals (H || c)
nexttile;
hold on;
for i = 1:size(Cevals, 2)
    plot(field, Cevals(:, i)-Cevals(1, 1), 'color', 'black', 'LineWidth', 1.2, 'DisplayName', 'my params mft');
    % plot(field, Cevals_nomft(:, i)-Cevals_nomft(:, 1), 'color', 'red', 'LineStyle', '--','LineWidth', 1.2, 'DisplayName', 'my params no mft');
    % plot(field, Cevals_allen(:, i)-Cevals_nomft(:, 1), 'color', 'cyan','LineWidth', 1.2, 'LineStyle', ':', 'DisplayName', 'allen params mft');
    % plot(field, Cevals_allen_nomft(:, i)+abs(ZFevals_allen(1)), 'LineWidth', 1.2, 'LineStyle', '--', 'DisplayName', 'allen params mft');
end
title('H || c');
xlabel('Field (T)');

% Set shared properties
for ax = 1:3
    nexttile(ax);
    xlim([0, 100]);
    ylim([-50, 50]);
end

% % Add overall title
% sgtitle(sprintf('B20: %.2f, B40: %.2f, B43: %.2f, B60: %.2f, B63: %.2f, B66: %.2f', ...
%     B20, B40, B43, B60, B63, B66), 'FontSize', 14);

%% checking J with different temperatures

% Load data from HDF5
fileName = 'jjz_temperature_test_dec_20_model.h5';
H = h5read(fileName, '/H');
Mdata2K = h5read(fileName, '/Mdata2K');
temps = h5read(fileName, '/temps');

% Extract groups for J values
info = h5info(fileName);
J_groups = {info.Groups.Name};  % Names of groups (e.g., '/J_-2.11e-3')

% Create one figure with subplots
figure;
numPlots = length(J_groups);  % Number of J values
cols = ceil(sqrt(numPlots));  % Number of columns for subplots
rows = ceil(numPlots / cols); % Number of rows for subplots

for i = 1:numPlots
    groupName = J_groups{i};
    J = h5readatt(fileName, groupName, 'J');
    mag = h5read(fileName, [groupName '/mag']);
    
    % Create subplot
    subplot(rows, cols, i);
    hold on;
    plot(Mdata2K(:, 1), Mdata2K(:, 2) / 1.37, 'bo', 'DisplayName', '2K MPMS Data');
    for j = 1:size(mag, 2)
        plot(H, mag(:, j), 'DisplayName', sprintf('Temp = %.1f K', temps(j)));
    end
    xlim([0 8])
    ylim([0 8])
    hold off;
    
    % Add labels, legend, and title
    xlabel('H [T]');
    ylabel('M [\mu_B]');
    title(sprintf('J = %.5f', J));
    grid on;
    if i == numPlots  % Add legend to the last subplot
        legend show;
    end
end

% Adjust layout for better spacing
sgtitle('Magnetization Curves for Different J Values');

%% plot by temperature instead of J
fileName = 'jjz_temperature_test_dec_20_model.h5';
H = h5read(fileName, '/H');
Mdata2K = h5read(fileName, '/Mdata2K');
temps = h5read(fileName, '/temps');

% Extract groups for J values
info = h5info(fileName);
J_groups = {info.Groups.Name};  % Names of groups (e.g., '/J_-2.11e-3')

% Preallocate for all J data
numTemps = length(temps);
numJ = length(J_groups);
allMag = zeros(numTemps, length(H), numJ);  % [Temperature x Field x J]

% Loop through each J group to read data
for jIdx = 1:numJ
    groupName = J_groups{jIdx};
    mag = h5read(fileName, [groupName '/mag']);  % [Temperature x Field]
    allMag(:, :, jIdx) = mag';
end

% Create one figure with subplots organized by temperature
figure;
cols = ceil(sqrt(numTemps));  % Number of columns for subplots
rows = ceil(numTemps / cols); % Number of rows for subplots

for tIdx = 1:numTemps
    % Create subplot for each temperature
    subplot(rows, cols, tIdx);
    hold on;
    plot(Mdata2K(:, 1), Mdata2K(:, 2) / 1.37, 'bo', 'DisplayName', '2K MPMS Data');
    
    % Plot magnetization for each J value
    for jIdx = 1:numJ
        groupName = J_groups{jIdx};
        J = h5readatt(fileName, groupName, 'J');
        plot(H, squeeze(allMag(tIdx, :, jIdx)), ...
             'DisplayName', sprintf('J = %.5f', J));
    end
    hold off;
    
    % Add labels, legend, and title
    xlabel('H [T]');
    ylabel('M [\mu_B]');
    title(sprintf('Temp = %.1f K', temps(tIdx)));
    grid on;
    if tIdx == numTemps  % Add legend to the last subplot
        legend show;
    end
end

% Adjust layout for better spacing
sgtitle('Magnetization Curves for dec 20 model');

%% make spec line fig
% Define the temperatures and corresponding file names
temps = [ 10];
colors = lines(length(temps));  % Distinct colors for each temp

figure;
tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% Plot H = 0 ZFevals from the first temperature
fileName = sprintf('EvsH_mft_T%dK_1015May12.h5', temps(1));
ZFevals = h5read(fileName, '/ZFevals');
nexttile;
hold on;
for i = 1:numel(ZFevals)
    yline(ZFevals(i)-ZFevals(1), 'LineWidth', 1.2);
end
title('H = 0');
ylabel('Energy');

% Loop over temps to read and plot ABevals and Cevals
for ti = 1:length(temps)
    T = temps(ti);
    fileName = sprintf('EvsH_mft_T%dK_1015May12.h5', T);
    field = h5read(fileName, '/field');
    ABevals = h5read(fileName, '/ABevals');
    Cevals = h5read(fileName, '/Cevals');

    % H || b
    nexttile(2);
    hold on;
    for i = 1:size(ABevals, 2)
        plot(field, ABevals(:, i) - ABevals(:, 1), 'Color', colors(ti,:), 'DisplayName',['to ', num2str(i)]);
    end
    title('H || b');
    xlabel('Field (T)');

    % H || c
    nexttile(3);
    hold on;
    for i = 1:size(Cevals, 2)
        plot(field, Cevals(:, i) - Cevals(:, 1), 'Color', colors(ti,:), 'DisplayName', ['to ', num2str(i)]);
    end
    title('H || c');
    xlabel('Field (T)');
end

% Set shared plot properties
for ax = 1:3
    nexttile(ax);
    xlim([0 10]);
    ylim([-5 30]);
    grid on;
end

% Add legends only once
nexttile(2); legend('Location', 'best');


%% let's plot the magnetotropic coeff

% Load data
filename = 'magnetotropic_phi_90_18jun2025.h5';

t_arr = h5read(filename, '/t_arr');
Hval = h5read(filename, '/Hval');
theta = h5read(filename, '/theta');
k_temp = h5read(filename, '/k_temp'); % size: [length(theta), length(Hval), length(t_arr)]
% NOTE: HDF5 saves in row-major, so permute dimensions to match Python's [t, H, theta]
k_temp = k_temp.r;
k_temp = permute(k_temp, [3, 2, 1]);

% Plot
nTemps = length(t_arr);
nH = length(Hval);
colors = turbo(nH);

figure('DefaultAxesFontSize',12);

tiledlayout(ceil(nTemps/2), 2, 'TileSpacing','compact');

for t_idx = 1:nTemps
    nexttile;
    
    deg = theta * 180 / pi;
    
    for h_idx = 1:nH
        plot(deg, squeeze(k_temp(t_idx, h_idx, :)), 'Color', colors(h_idx,:), 'LineWidth', 1.2);
        hold on;
    end
    
    title(sprintf('T = %.1f K', t_arr(t_idx)));
    xlabel('angle [deg]');
    ylabel('k');
    grid on;
    
    if t_idx == 1
        legend(arrayfun(@(hval) sprintf('%.1f T', hval), Hval, 'UniformOutput', false));
    end
end

%%
filename = 'spectroscopy_fgr_all_2025Jul30.h5';

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
% Load data from CSV files
% ramanData = readmatrix('ramanData.csv');  % assuming saved separately
% binaryAmp = readmatrix('binary_amp.csv');
ampC = my_spec_data.ampC;%
% ampC = readmatrix('amplitude.csv');
arrC = my_spec_data.linesC;                       %
% arrC = readmatrix('energies.csv');
calc_field = my_spec_data.calc_field;%
% calc_field = readmatrix('calc_field.csv');

% Generate field and wavenumber axes
% [n_wavenum, n_field] = size(ramanData);
wavenums = linspace(0, 120, length(ampC));  % adjust as appropriate
field = calc_field;        % adjust as appropriate

% Create figure and contour plot
figure;
hold on;
% contourf(field, wavenums, ramanData, 100, 'LineColor', 'none');
% colormap('Oranges');
% colorbar;
% xlim([0, 14]);
% ylim([0, 120]);
xlabel('Field (T)');
ylabel('Energy (cm^{-1})');

% Overlay colored lines where binaryAmp == 1
n_lines = size(arrC, 1);
for i = 2:55
    mask = ampC(:, i)>0;
    if any(mask)
        if i<17
            plot(calc_field(mask), arrC(mask, i)*8.022, 'b-', 'LineWidth', 1.5, 'DisplayName', ['0 to' num2str(i-1)]); % 'c' = cyan
        end
        if i>=17 && i<32
            plot(calc_field(mask), arrC(mask, i)*8.022, 'b--', 'LineWidth', 1.5, 'DisplayName',['1 to' num2str(i-15)] ); % 'c' = cyan
        end
        if i>=31 && i<44
            plot(calc_field(mask), arrC(mask, i)*8.022, 'b:', 'LineWidth', 1.5, 'DisplayName',['2 to' num2str(i-15-14+1)] ); % 'c' = cyan
        end

        if i>=43 && i<54
            plot(calc_field(mask), arrC(mask, i)*8.022, 'b:', 'LineWidth', 1.5, 'DisplayName',['3 to' num2str(i-15-14-13+3)] ); % 'c' = cyan
        end
    end
end

title('H||c, selection rulesmasked for amplitude >1e-2  computational error means we actually never get 0');
hold off;
%%
figure(); 
contourf(calc_field, my_spec_data.calc_wavenums, my_spec_data.simulated_IR_B./max(my_spec_data.simulated_IR_B), 100, 'linestyle', 'none');
colormap(coolwarm)

%%

simdata = (my_spec_data.simulated_C(:,2:end)-my_spec_data.simulated_C(:,1)); 
% simdata = simdata./max(max(simdata));
avgdata = zeros([length(simdata(:,1)),1]);
for i =1:length(simdata(1,:)) 
    avgdata = avgdata + simdata(:,i);
end
simdata = simdata - avgdata./length(simdata(1,:));
simdata = simdata./max(max(simdata));
figure(); 
contourf(calc_field(2:end), my_spec_data.calc_wavenums, simdata, 100, 'linestyle', 'none');
colormap(coolwarm)
xlabel('H [T]')
ylabel('E [cm^{-1}]')
clim([-.1, .1])
%% make sim figures
figure; hold on; box on;
contourf(my_spec_data.calc_field, my_spec_data.calc_wavenums, my_spec_data.simulated_IR_B, 100, 'linestyle', 'none')

%% plot IR data to get phonons
bfname = '/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/Maglab_IRMeasurement_June2022/ProcessingCode/Load1_TrimData/P2_CsEr_100-FIR_RAWAVG.dat';
cfname = '/Users/hopeless/Desktop/LeeLab/data/CsErSe2_data/Maglab_IRMeasurement_June2022/ProcessingCode/Load2_TrimData/P3_CsEr_100_RAWAVG.dat';


%%
% figure; hold on; box on; grid on; 
for i = 1:83
    plot(wavenumber, B_IR_data(:,i)+i*.0)
end
plot(wav, pinhole, 'mo--')

%% 
% take ZF spec, subtract the pinhole
i = 1;
figure; hold on; box on; grid on; 
plot(wavenumber, B_IR_data(:,i)./max(B_IR_data(:,i)));
plot(wav, pinhole./max(pinhole), 'mo--')
% plot(wavenumber, B_IR_data(:,i)./interp_pinhole)
%% make average spectrum 
data = zeros(length(B_IR_data), 1);
for i = 1:length(B_IR_data)
    % disp(i)
    data(i) = sum(B_IR_data(i,:),'omitnan')/length(~isnan(B_IR_data(i,:)))+ sum(C_IR_data(i,:),'omitnan')/length(~isnan(C_IR_data(i,:)));
end 
%%
figure; hold on; box on; grid on; 
plot(wavenumber, data./max(data))
plot(wav, pinhole./6.2, 'm--')
%% interp, divide by pinhole

% interp_data = interp1(wavenumber, data/max(data), linspace(min(wavenumber),max(wavenumber), ))
interp_pinhole = interp1(wav, pinhole./6.2, wavenumber); 
data = data./max(data);
%%
div_data = data./interp_pinhole; 

%%
figure; hold on; box on; grid on; 
plot(wavenumber, 1-(div_data./max(div_data)))

%% 
usable = div_data(400:end); 
wave = wavenumber(400:end);

%%
figure; hold on; box on; grid on; 
plot(wave, 1-(usable./max(usable)))

%% rq, let's look at just ZF

simAB = h5read('simAB_data.h5', '/simAB');
field = h5read('simAB_data.h5', '/field');
wavenum = h5read('simAB_data.h5', '/wavenum');

% Create the contour plot
figure;
contourf(field, wavenum, simAB, 100, 'LineColor', 'none');  % Transpose for orientation
xlabel('H [T]');
ylabel('E [cm^{-1}]');
title('Simulated Zeeman Spectrum (AB Plane)');
colorbar;

%% 
figure; hold on; box on; grid on; 
plot(wavenumber, ir(:,1), "DisplayName", num2str(field(1))); 
plot(wavenumber, ir(:,15), "DisplayName", num2str(field(15))); 
plot(wavenumber, ir(:,43), "DisplayName", num2str(field(43))); 
plot(wavenumber, ir(:,71), "DisplayName", num2str(field(71))); 
legend()
xlabel('E [cm^{-1}]')
ylabel('Transmission [arb]')

%% subtract zf
figure; hold on; box on; grid on; 
% plot(wavenumber, ir(:,1), "DisplayName", num2str(field(1))); 
plot(wavenumber, ir(:,15)-ir(:,1), "DisplayName", num2str(field(15))); 
plot(wavenumber, ir(:,43)-ir(:,1), "DisplayName", num2str(field(43))); 
plot(wavenumber, ir(:,71)-ir(:,1), "DisplayName", num2str(field(71))); 
legend()
xlabel('E [cm^{-1}]')
ylabel('Transmission [arb]')

%% figure
figure; box on;
contourf(field(2:71), wavenumber, -1*(ir(:,2:71)-ir(:,1)),100, 'linestyle', 'none')
xlabel('H[T]')
ylabel('E [cm^{-1}]')
clim([0 1])
%% divide
figure; hold on; box on; grid on; 
% plot(wavenumber, ir(:,1), "DisplayName", num2str(field(1))); 
plot(wavenumber, ir(:,15)./ir(:,72), "DisplayName", num2str(field(15))); 
plot(wavenumber, ir(:,43)./ir(:,72), "DisplayName", num2str(field(43))); 
plot(wavenumber, ir(:,71)./ir(:,72), "DisplayName", num2str(field(71))); 
legend()
xlabel('E [cm^{-1}]')
ylabel('Transmission [arb]')



%% figure
figure; box on;
contourf(field(1:71), wavenumber, -1-(ir(:,1:71)./ir(:,72)),100, 'linestyle', 'none')
xlabel('H[T]')
ylabel('E [cm^{-1}]')
% clim([0 1])

%% figure
figure; box on;
contourf(field_b(2:40), wavenumber, 1-(ir_b(:,2:40)-ir_b(:,1))+1,256, 'linestyle', 'none')
xlabel('H[T]')
ylabel('E [cm^{-1}]')
% clim([0 1])