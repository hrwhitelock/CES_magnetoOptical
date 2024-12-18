%% first we plot jst the lines with no data
% Load data from HDF5
fileName = 'zeeman_split_lines.h5';
fieldArr = h5read(fileName, '/fieldArr');
arrC_neg = h5read(fileName, '/arrC_neg');
arrC_zero = h5read(fileName, '/arrC_zero');
B20 = h5readatt(fileName, '/', 'B20');
B40 = h5readatt(fileName, '/', 'B40');
B43 = h5readatt(fileName, '/', 'B43');
B60 = h5readatt(fileName, '/', 'B60');
B63 = h5readatt(fileName, '/', 'B63');
B66 = h5readatt(fileName, '/', 'B66');

% Plot configurations
figure;
hold on;
grid on;
xlim([0, 14]);
ylim([0, 120]);

% Add negative Jz splitting data

% now add data
filename = 'CsErSe2_IR_c_axis_line_calculation.h5';

data = h5read(filename, '/data');
wavenums = h5read(filename, '/wavenums');
fields = str2double(h5read(filename, '/fields'));

contourf(fields, wavenums, data', 100, 'LineStyle', 'none')
for i = 1:length(arrC_neg(1,:))
    if i <= 16
        plot(fieldArr, arrC_neg(:, i), 'r', 'DisplayName', sprintf('Line %d (Jz = -2.53e-3)', i), 'LineWidth', 0.8, 'Color', [1, 0, 0, 0.7]);
    end
end
xlabel('H [T]')
xlim([0 18])
ylabel('Energy [cm^{-1}]')
title('H || c Allens params, JJz = 2.53ueV')
%% Add zero Jz splitting data
figure(); 
hold on; 
box on; 
contourf(fields, wavenums, data', 100, 'LineStyle', 'none')
for i = 1:1:length(arrC_zero(1,:))
    if i <= 16
        plot(fieldArr, arrC_zero(:, i), 'r--', 'DisplayName', sprintf('Line %d (Jz = 0)', i), 'LineWidth', 0.8, 'Color', [1, 0, 0, 0.7]);
    end
end

% Add title and labels
title('H || c Allens params, no mft')
xlabel('Field (T)', 'FontSize', 10);
ylabel('Energy (cm^{-1})', 'FontSize', 10);
ylim([0 120])
hold off;


%% now we do same for AB plane
%% first we plot jst the lines with no data
% Load data from HDF5
fileName = 'zeeman_split_lines_AB.h5';
fieldArr = h5read(fileName, '/fieldArr');
arrC_neg = h5read(fileName, '/arrB_neg');
arrC_zero = h5read(fileName, '/arrB_zero');
B20 = h5readatt(fileName, '/', 'B20');
B40 = h5readatt(fileName, '/', 'B40');
B43 = h5readatt(fileName, '/', 'B43');
B60 = h5readatt(fileName, '/', 'B60');
B63 = h5readatt(fileName, '/', 'B63');
B66 = h5readatt(fileName, '/', 'B66');
% Jperp = h5readatt(fileName, '/', 'Jperp')
% Plot configurations
figure;
hold on;
grid on;
xlim([0, 14]);
ylim([0, 120]);

% Add negative Jz splitting data

% now add data
filename = 'CsErSe2_IR_b_axis_simulation.h5';

data = h5read(filename, '/data');
wavenums = h5read(filename, '/wavenums');
fields = str2double(h5read(filename, '/fields'));


% contourf(fields, wavenums, data', 100, 'LineStyle', 'none')
for i = 1:length(arrC_neg(1,:))
    if i <= 16
        plot(fieldArr, arrC_neg(:, i), 'r', 'DisplayName', sprintf('Line %d (Jz = -2.53e-3)', i), 'LineWidth', 0.8, 'Color', [1, 0, 0, 0.7]);
    end
end
xlabel('H [T]')
xlim([0 18])
ylabel('Energy [cm^{-1}]')
title('H || c Allens params, JJperp = -.575ueV')
% %% Add no mft
% figure(); 
hold on; 
box on; 
% contourf(fields, wavenums, data', 100, 'LineStyle', 'none')
for i = 1:1:length(arrC_zero(1,:))
    if i <= 16
        plot(fieldArr, arrC_zero(:, i), 'r--', 'DisplayName', sprintf('Line %d (Jz = 0)', i), 'LineWidth', 0.8, 'Color', [1, 0, 0, 0.7]);
    end
end

% Add title and labels
title('H || b Allens params, no mft')
xlabel('Field (T)', 'FontSize', 10);
ylabel('Energy (cm^{-1})', 'FontSize', 10);
ylim([0 120])
hold off;