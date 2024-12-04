
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% starting with raman data
%% Load HDF5 data
filename = 'CsErSe2_raman_line_calculation.h5';

ramanData = h5read(filename, '/ramanData');
wavenums = h5read(filename, '/wavenums');
fields = str2double(h5read(filename, '/fields'));
waveArr = h5read(filename, '/waveArr');
fieldArr = h5read(filename, '/fieldArr');
ampC = h5read(filename, '/ampC');
arrC = h5read(filename, '/arrC');

B20 = h5readatt(filename, '/', 'B20');
B40 = h5readatt(filename, '/', 'B40');
B43 = h5readatt(filename, '/', 'B43');
B60 = h5readatt(filename, '/', 'B60');
B63 = h5readatt(filename, '/', 'B63');
B66 = h5readatt(filename, '/', 'B66');

%% Plot 1: Raman Data (0-100 cm^-1)
figure;
hold on;
colors = flipud(cool(length(fields)));
j= 0;
for i = length(fields):-1:1
    plot(wavenums, ramanData(i, :) + j/3, 'Color', colors(i, :));
    if i == length(fields)
        text(90, j/3 + 0.15, 'H = 0T', 'FontSize', 9);
    elseif i == 1
        text(90, j/3 + 0.15, ['H = ', num2str(field(i)), 'T'], 'FontSize', 9);
    end
    j = j+1; 
end
xlabel('Energy (cm^{-1})');
ylabel('Intensity (arb)');
title('Raman Data, full range');
% xlim([0, 100]);
hold off;

%% Plot 2: Zeeman Splitting Lines Overlaid on waterfall plot
figure;
hold on;
j= 0;
for i = length(fields):-1:1
    plot(wavenums, ramanData(i, :) + fields(i), 'Color', colors(i, :));
    if i == length(fields)
        text(90, fields(i) + 0.15, 'H = 0T', 'FontSize', 9);
    elseif i == 1
        text(90, fields(i) + 0.15, ['H = ', num2str(field(i)), 'T'], 'FontSize', 9);
    end
    j = j+1; 
end
for i = 1:17%size(arrC, 2)
    plot(arrC(:, i), fieldArr, 'r', 'LineWidth', 0.7);
end
for i = 17:size(arrC,2)
    plot( arrC(:, i),fieldArr, 'r--', 'LineWidth', 0.7);
end
colorbar;
xlabel('Field (T)');
ylabel('Energy (cm^{-1})');
title(['Zeeman Splitting with Overlay (B20: ', num2str(B20), ...
    ', B40: ', num2str(B40), ', B43: ', num2str(B43), ...
    ', B60: ', num2str(B60), ', B63: ', num2str(B63), ...
    ', B66: ', num2str(B66), ')', 'I havent figured out how to make fading lines in matlab yet']);
hold off;

%% Plot 3: Simulated Data
%Load HDF5 data
filename = 'CsErSe2_raman_simulation.h5';

ramanData = h5read(filename, '/ramanData');
wavenums = h5read(filename, '/wavenums');
fields = str2double(h5read(filename, '/fields'));
waveArr = h5read(filename, '/waveArr');
fieldArr = h5read(filename, '/fieldArr');
simulatedRaman = h5read(filename, '/simulatedData'); 

B20 = h5readatt(filename, '/', 'B20');
B40 = h5readatt(filename, '/', 'B40');
B43 = h5readatt(filename, '/', 'B43');
B60 = h5readatt(filename, '/', 'B60');
B63 = h5readatt(filename, '/', 'B63');
B66 = h5readatt(filename, '/', 'B66');

fig = figure;
subplot(1,2,1)
contourf(fields, wavenums, ramanData', 100, 'LineStyle', 'none')
ylim([0 120])
xlabel('Field (T)')
ylabel('Energy (cm^{-1})');
title('Raman Data')
subplot(1,2,2)
contourf(fieldArr, waveArr, simulatedRaman, 100, 'LineStyle', 'none');
colormap("jet")
colorbar;
xlabel('Field (T)');
ylabel('Energy (cm^{-1})');
title('Simulated Data');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% next do IR c axis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load from hdf5
filename = 'CsErSe2_IR_c_axis_line_calculation.h5';

data = h5read(filename, '/data');
wavenums = h5read(filename, '/wavenums');
fields = str2double(h5read(filename, '/fields'));
waveArr = h5read(filename, '/waveArr');
fieldArr = h5read(filename, '/fieldArr');
ampC = h5read(filename, '/ampC');
arrC = h5read(filename, '/arrC');

B20 = h5readatt(filename, '/', 'B20');
B40 = h5readatt(filename, '/', 'B40');
B43 = h5readatt(filename, '/', 'B43');
B60 = h5readatt(filename, '/', 'B60');
B63 = h5readatt(filename, '/', 'B63');
B66 = h5readatt(filename, '/', 'B66');

%% Plot 1: Raman Data (0-100 cm^-1)
figure;
hold on;
colors = flipud(cool(length(fields)));
for i = length(fields):-1:1
    plot(wavenums, data(i, :) + fields(i), 'Color', colors(i, :));
    if i == length(fields)
        text(90, fields(i) + 0.15, 'H = 0T', 'FontSize', 9);
    elseif i == 1
        text(90, fields(i) + 0.15, ['H = ', num2str(field(i)), 'T'], 'FontSize', 9);
    end 
end
xlabel('Energy (cm^{-1})');
ylabel('Intensity (arb)');
title('IR Data, c-axis, full range');
% xlim([0, 100]);
hold off;

%% Plot 2: Zeeman Splitting Lines Overlaid on waterfall plot
figure;
hold on;
j= 0;
for i = length(fields):-1:1
    plot(wavenums, data(i, :) + fields(i), 'Color', colors(i, :));
    if i == length(fields)
        text(90, fields(i) + 0.15, 'H = 0T', 'FontSize', 9);
    elseif i == 1
        text(90, fields(i) + 0.15, ['H = ', num2str(field(i)), 'T'], 'FontSize', 9);
    end
    j = j+1; 
end
for i = 1:17%size(arrC, 2)
    plot(arrC(:, i), fieldArr, 'r', 'LineWidth', 0.7);
end
for i = 17:size(arrC,2)
    plot( arrC(:, i),fieldArr, 'r--', 'LineWidth', 0.7);
end
% colorbar;
xlabel('Field (T)');
ylabel('Energy (cm^{-1})');
title(['IR c-axis with Overlay (B20: ', num2str(B20), ...
    ', B40: ', num2str(B40), ', B43: ', num2str(B43), ...
    ', B60: ', num2str(B60), ', B63: ', num2str(B63), ...
    ', B66: ', num2str(B66), ')', 'I havent figured out how to make fading lines in matlab yet']);
ylim([0 18])
hold off;

%% Plot 3: Simulated Data
%Load HDF5 data
filename = 'CsErSe2_IR_c_axis_simulation.h5';

data = h5read(filename, '/data');
wavenums = h5read(filename, '/wavenums');
fields = str2double(h5read(filename, '/fields'));
waveArr = h5read(filename, '/waveArr');
fieldArr = h5read(filename, '/fieldArr');
simulatedIR = h5read(filename, '/simulatedData'); 

B20 = h5readatt(filename, '/', 'B20');
B40 = h5readatt(filename, '/', 'B40');
B43 = h5readatt(filename, '/', 'B43');
B60 = h5readatt(filename, '/', 'B60');
B63 = h5readatt(filename, '/', 'B63');
B66 = h5readatt(filename, '/', 'B66');

fig = figure;
subplot(1,2,1)
contourf(fields, wavenums, data', 100, 'LineStyle', 'none')
ylim([0 120])
clim([0 1])
xlabel('Field (T)')
ylabel('Energy (cm^{-1})');
title('IR c-axis Data')
subplot(1,2,2)
contourf(fieldArr, waveArr, simulatedIR, 100, 'LineStyle', 'none');
colormap("jet")
clim([0 1])
colorbar;
xlabel('Field (T)');
ylabel('Energy (cm^{-1})');
title('Simulated Data');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% next do IR c axis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load from hdf5
filename = 'CsErSe2_IR_b_axis_line_calculation.h5';

data = h5read(filename, '/data');
wavenums = h5read(filename, '/wavenums');
fields = str2double(h5read(filename, '/fields'));
waveArr = h5read(filename, '/waveArr');
fieldArr = h5read(filename, '/fieldArr');
ampAB = h5read(filename, '/ampAB');
arrAB = h5read(filename, '/arrAB');

B20 = h5readatt(filename, '/', 'B20');
B40 = h5readatt(filename, '/', 'B40');
B43 = h5readatt(filename, '/', 'B43');
B60 = h5readatt(filename, '/', 'B60');
B63 = h5readatt(filename, '/', 'B63');
B66 = h5readatt(filename, '/', 'B66');

%% Plot 1: waterfall Data (0-100 cm^-1)
figure;
hold on;
colors = flipud(cool(length(fields)));
for i = length(fields):-1:1
    plot(wavenums, data(i, :) + fields(i), 'Color', colors(i, :));
    if i == length(fields)
        text(90, fields(i) + 0.15, ['H = ', num2str(fields(i)), 'T'], 'FontSize', 9);
    elseif i == 1
        text(90, fields(i) + 0.15, ['H = ', num2str(fields(i)), 'T'], 'FontSize', 9);
    end 
end
xlabel('Energy (cm^{-1})');
ylabel('Intensity (arb)');
title('IR Data, b-axis, full range');
% xlim([0, 100]);
hold off;

%% Plot 2: Zeeman Splitting Lines Overlaid on waterfall plot
figure;
hold on;
j= 0;
for i = length(fields):-1:1
    plot(wavenums, data(i, :) + fields(i), 'Color', colors(i, :));
    if i == length(fields)
        text(90, fields(i) + 0.15, ['H = ', num2str(fields(i)), 'T'], 'FontSize', 9);
    elseif i == 1
        text(90, fields(i) + 0.15, ['H = ', num2str(fields(i)), 'T'], 'FontSize', 9);
    end
    j = j+1; 
end
for i = 1:17%size(arrC, 2)
    plot(arrAB(:, i), fieldArr, 'r', 'LineWidth', 0.7);
end
for i = 17:size(arrAB,2)
    plot( arrAB(:, i),fieldArr, 'r--', 'LineWidth', 0.7);
end
% colorbar;
xlabel('Field (T)');
ylabel('Energy (cm^{-1})');
title(['IR b-axis with Overlay (B20: ', num2str(B20), ...
    ', B40: ', num2str(B40), ', B43: ', num2str(B43), ...
    ', B60: ', num2str(B60), ', B63: ', num2str(B63), ...
    ', B66: ', num2str(B66), ')', 'I havent figured out how to make fading lines in matlab yet']);
ylim([0 18])
hold off;

%% Plot 3: Simulated Data
%Load HDF5 data
filename = 'CsErSe2_IR_b_axis_simulation.h5';

data = h5read(filename, '/data');
wavenums = h5read(filename, '/wavenums');
fields = str2double(h5read(filename, '/fields'));
waveArr = h5read(filename, '/waveArr');
fieldArr = h5read(filename, '/fieldArr');
simulatedIR = h5read(filename, '/simulatedData'); 

B20 = h5readatt(filename, '/', 'B20');
B40 = h5readatt(filename, '/', 'B40');
B43 = h5readatt(filename, '/', 'B43');
B60 = h5readatt(filename, '/', 'B60');
B63 = h5readatt(filename, '/', 'B63');
B66 = h5readatt(filename, '/', 'B66');

fig = figure;
subplot(1,2,1)
contourf(fields, wavenums, data', 100, 'LineStyle', 'none')
ylim([0 120])
clim([0 1])
xlabel('Field (T)')
ylabel('Energy (cm^{-1})');
title('IR b-axis Data')
subplot(1,2,2)
contourf(fieldArr, waveArr, simulatedIR, 100, 'LineStyle', 'none');
colormap("jet")
clim([0 1])
colorbar;
xlabel('Field (T)');
ylabel('Energy (cm^{-1})');
title('Simulated Data');