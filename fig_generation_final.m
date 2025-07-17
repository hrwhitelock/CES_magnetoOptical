% final fig generation script 12/16/24

% let's start with making the paper figures
%% load my spectroscopy data
filename = 'spectroscopy_calculation_hopes_params_2025Jan30_02.h5';

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
yellowMap = [linspace(0, 255, 256)', linspace(0, 255, 256)', zeros(256, 1)];
colormap(yellowMap);

%% make spec without thermal lines: 
% Load data (you probably already have this loaded)
field_vals = my_spec_data.calc_field;   % e.g., size [Nfield × 1]
linesB = my_spec_data.linesB;           % e.g., size [Nfield × Nlines]

% Define energy axis
Emin = 0; Emax = 100;  % cm^-1 range
dE = 0.05;
E = Emin:dE:Emax;      % energy grid
Nfield = length(field_vals);
Nenergy = length(E);
Z = zeros(Nenergy, Nfield);

% Gaussian width in cm^-1
sigma = 0.5;

% Loop through each spectral line
for i = 2:16
    for f = 1:Nfield
        E0 = linesB(f, i);  % Center of the Gaussian
        if isnan(E0), continue; end  % Skip if missing

        % Gaussian over energy axis
        G = exp(-((E - E0).^2) / (2*sigma^2));

        % Add to spectrum at field f
        Z(:, f) = Z(:, f) + G';
    end
end

% Normalize (optional)
Z = Z / max(Z(:));

% Plot
figure;
imagesc(field_vals, E, Z);
axis xy;
xlabel('Field [T]');
ylabel('Energy [cm^{-1}]');
title('Simulated 2D Spectrum from Lines');
colormap(jet);
colorbar;

%% 
fig = figure;
ax1 = subplot(2,2,1);
title('B-axis IR')
hold on; 
pcolor(my_spec_data.IR_B_field,my_spec_data.IR_B_wavenums,my_spec_data.IR_dataB')
% axis xy;
shading flat
for i = 2:16%size(arrC, 2)
    plot(my_spec_data.calc_field, my_spec_data.linesB(:,i),  'r--', 'LineWidth', 1);
end
set(ax1,'Xticklabel',[])
ylabel({'data', 'Energy [cm{^-1}]'})
clim([0 1])
ylim([0 100])
xlim([0 17.5])
colormap(ax1, jet)
set(gca,'Layer','top')

ax3 = subplot(2,2,3);
hold on; 
idx = 1:1:length(my_spec_data.calc_wavenums);
field_idx = 1:1:length(my_spec_data.calc_field); 
contourf(my_spec_data.calc_field(field_idx),my_spec_data.calc_wavenums(idx), my_spec_data.simulated_IR_B(idx, field_idx), 100,'LineStyle', 'none');
% contourf(np.linspace(0,18.5),np.linspace(0, 100, 10000), Vq)
colormap(ax3, jet)
clim([0 1])
% xlabel({'calculation', 'Energy [cm{^-1}]'})
% ylabel('H(T)')

for i = 33:50%size(arrC, 2)
    plot(my_spec_data.calc_field, my_spec_data.linesB(:,i), 'r--', 'LineWidth', 1);
end

clim([0 1])
ylim([0 100])
xlim([0 17.5])
set(gca,'Layer','top')
% 

% linkaxes([ax1, ax2], 'xy')

% %% c axis IR
% figure; hold on; 
ax2 = subplot(2,2,2);hold on; box on; 
pcolor(my_spec_data.IR_C_field, my_spec_data.IR_C_wavenums, my_spec_data.IR_dataC')
shading flat
for i = 2:16 %size(my_spec_data.linesC, 2)
    plot(my_spec_data.calc_field(field_idx), my_spec_data.linesC(field_idx,i),  'r--', 'LineWidth',1);
end
title('C-axis IR')
ylabel('Energy [cm{^-1}]')
colormap(ax2, jet)
clim([0 1])
ylim([0 100])
xlim([0 17.5])
set(gca,'Layer','top')
% 
ax4 = subplot(2,2,4);
hold on; box on; 
contourf(my_spec_data.calc_field(field_idx),my_spec_data.calc_wavenums(idx), my_spec_data.simulated_IR_C(idx, field_idx), 100, 'LineStyle', 'none');
colormap(ax4, jet)
for i = 2:16%size(arrC, 2)
    plot(my_spec_data.calc_field(field_idx), my_spec_data.linesC(field_idx,i), 'r--', 'LineWidth', 1);
end
clim([0 1])
ylim([0 100])
xlim([0 17.5])
ylabel({'Data', 'Energy [cm{^-1}]'})
xlabel('H(T)')
set(gca,'Layer','top')

ax1.Position(2)=ax3.Position(2)+ax3.Position(4)+.002; 
ax3.Position(1)=ax1.Position(1);
ax2.Position(1) = ax1.Position(1) +ax1.Position(3) +.002;

ax2.Position(2)=ax1.Position(2);
ax4.Position(1)=ax2.Position(1);

linkaxes([ax1, ax2, ax3, ax4], 'xy')
set(ax2,'Yticklabel',[]) 
set(ax1,'Xticklabel',[])
set(ax2,'Xticklabel',[])
set(ax4,'Yticklabel',[])


%% c axis raman

fig = figure;
ax1 = subplot(2,1,1);box on; hold on; 
pcolor(my_spec_data.raman_field, my_spec_data.raman_wavenums, my_spec_data.ramanData'); 
shading flat; 
for i = 33:45%size(arrC, 2)
    plot(my_spec_data.calc_field(field_idx), my_spec_data.linesC(field_idx,i), 'b--', 'LineWidth', 1);
end
% for i = 17:50%size(arrC, 2)
%     plot(my_spec_data.calc_field,my_spec_data.linesC(:,i),  'r--', 'LineWidth', 1);
% end
ylim([0 100])
xlim([0 14])
colormap(ax1, cm)
title ('C-axis Raman')
ylabel('Energy [cm^{-1}]')

ax2 = subplot(2,1,2);
hold on; box on; 
contourf(my_spec_data.calc_field(field_idx),my_spec_data.calc_wavenums(idx), my_spec_data.simulated_raman(idx, field_idx), 100, 'LineStyle', 'none');
colormap(ax2, cm)
% let's add some lines
for i = 2:16%size(arrC, 2)
    plot(my_spec_data.calc_field(field_idx), my_spec_data.linesC(field_idx,i), 'r--', 'LineWidth', 1);
end
clim([0 1])
xlim([0 14])
ylim([0 100])
ylabel('Energy [cm{^-1}]')
xlabel('H(T)')

ax1.Position(2)=ax2.Position(2)+ax2.Position(4)+.02;
ax2.Position(1)=ax1.Position(1);
linkaxes([ax1, ax2], 'xy')
set(ax1,'Xticklabel',[])

% 
% ax3.Position(2)=ax6.Position(2)+ax6.Position(4)+.02;
% ax1.Position(2)=ax4.Position(2)+ax4.Position(4)+.02;
% % 
% % ax4.Position(2)=ax6.Position(2)+ax6.Position(4);
% ax2.Position(2)=ax5.Position(2)+ax5.Position(4)+.02;
% % 
% % 
% % ax4.Position(1)=ax3.Position(1)+ax3.Position(3);
% ax2.Position(1)=ax1.Position(1)+ax1.Position(3);
% ax5.Position(1)=ax4.Position(1)+ax4.Position(3);
% linkaxes([ax1, ax2, ax3, ax4, ax5, ax6], 'xy')
% % ylabel('Energy [cm^{-1}]')
% % xlabel('H [t]')
% set(ax2,'Yticklabel',[]) 
% set(ax1,'Xticklabel',[])
% set(ax2,'Xticklabel',[])
% 
% set(ax5,'Yticklabel',[]) 
% % set(ax5,'Xticklabel',[])
% 
% 
% set(ax3,'Xticklabel',[])

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
fileName = 'EvsH_noNorm_hopes_params_2025Jan14.h5';
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


fileName = 'EvsH_noNorm_allens_params_2025Jan14.h5';
ZFevals_allen = h5read(fileName, '/ZFevals');
ABevals_allen = h5read(fileName, '/ABevals');
Cevals_allen = h5read(fileName, '/Cevals');
ABevals_allen_nomft = h5read(fileName, '/ABevals_nomft');
Cevals_allen_nomft = h5read(fileName, '/Cevals_nomft');
field = h5read(fileName, '/field');
B20_allen = h5readatt(fileName, '/', 'B20');
% Create the figure with subplots
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
    plot(field, ABevals(:, i)-ABevals(:, 1), 'color', 'black','LineWidth', 1.2, 'DisplayName', 'my params mft');
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
    plot(field, Cevals(:, i)-Cevals(:, 1), 'color', 'black', 'LineWidth', 1.2, 'DisplayName', 'my params mft');
    % plot(field, Cevals_nomft(:, i)-Cevals_nomft(:, 1), 'color', 'red', 'LineStyle', '--','LineWidth', 1.2, 'DisplayName', 'my params no mft');
    % plot(field, Cevals_allen(:, i)-Cevals_nomft(:, 1), 'color', 'cyan','LineWidth', 1.2, 'LineStyle', ':', 'DisplayName', 'allen params mft');
    % plot(field, Cevals_allen_nomft(:, i)+abs(ZFevals_allen(1)), 'LineWidth', 1.2, 'LineStyle', '--', 'DisplayName', 'allen params mft');
end
title('H || c');
xlabel('Field (T)');

% Set shared properties
for ax = 1:3
    nexttile(ax);
    xlim([0, 10]);
    ylim([-5, 30]);
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
temps = [1, 2, 4, 6, 10, 20];
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
        plot(field, ABevals(:, i) - ABevals(:, 1), 'Color', colors(ti,:), 'DisplayName', sprintf('%d K', T));
    end
    title('H || b');
    xlabel('Field (T)');

    % H || c
    nexttile(3);
    hold on;
    for i = 1:size(Cevals, 2)
        plot(field, Cevals(:, i) - Cevals(:, 1), 'Color', colors(ti,:), 'DisplayName', sprintf('%d K', T));
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
