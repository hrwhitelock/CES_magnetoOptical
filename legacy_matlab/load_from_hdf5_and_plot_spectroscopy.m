
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
for i = length(fields):-5:1
    plot(wavenums, ramanData(i, :) + j/3, 'Color', colors(i, :));
    if i == length(fields)
        text(90, j/3 + 0.15, 'H = 0T', 'FontSize', 9);
    elseif i == 1
        text(90, j/3 + 0.15, 'H = 14T', 'FontSize', 9);
    end
    j = j+1; 
end
xlabel('Energy (cm^{-1})');
ylabel('Intensity (arb)');
title('Raman Data, full range');
% xlim([0, 100]);
hold off;
% add raman lines
xline(46.1, 'g--') %E2g
xline(56.1, 'r--') %E2g
xline(91.1, 'r--') %E2g
xline(119.1, 'r--') %E2g
xline(125.9, 'g--') %E2g
xline(127.3, 'g--') %E2g
xline(164.8, 'g--') %E2g
xline(168.8, 'r--') %E2g


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
% next do IR b axis
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% let's make the six panelled figure 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% let's go AB plane IR, C-axis IR, c-axis Raman
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
ax1 = subplot(2,3,1);
contourf(fields, wavenums, data', 100, 'LineStyle', 'none')
ylim([0 120])
clim([0 1])
colormap(ax1, flipud(parula))
ax4 = subplot(2,3,4);
hold on; 
contourf(fieldArr, waveArr, simulatedIR, 100, 'LineStyle', 'none');
colormap(ax4, flipud(parula))
clim([0 1])
% add calculated lines

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
for i = 1:17%size(arrC, 2)
    plot(fieldArr, arrAB(:, i),  'r:', 'LineWidth', 0.5);
end
% for i = 17:size(arrAB,2)
%     plot( fieldArr,arrAB(:, i), 'r--', 'LineWidth', 0.7);
% end
ylim([0,100])
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

% fig = figure;
ax2 = subplot(2,3,2);
contourf(fields, wavenums, data', 100, 'LineStyle', 'none')
colormap(ax2, flipud(parula))
ylim([0 120])
clim([0 1])

ax5 = subplot(2,3,5);
hold on; 
contourf(fieldArr, waveArr, simulatedIR, 100, 'LineStyle', 'none');
colormap(ax5, flipud(parula))
clim([0 1])
% add lines 
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

for i = 1:17%size(arrC, 2)
    plot(fieldArr, arrC(:, i),  'r:', 'LineWidth', 0.5);
end


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

% fig = figure;
ax3 = subplot(2,3,3);
contourf(fields, wavenums, ramanData', 100, 'LineStyle', 'none')
ylim([0 120])
xlim([0 17.5])
colormap(ax3, 'cool')

ax6 = subplot(2,3,6);
hold on; 
contourf(fieldArr, waveArr, simulatedRaman, 100, 'LineStyle', 'none');
colormap(ax6, 'cool')
% let's add some lines
clim([0 1])

xlim([0 17.5])
ylim([0 100])

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

for i = 1:17%size(arrC, 2)
    plot(fieldArr, arrC(:, i),  'k:', 'LineWidth', 0.5);
end

ax3.Position(2)=ax6.Position(2)+ax6.Position(4);
ax1.Position(2)=ax4.Position(2)+ax4.Position(4);
% 
% ax4.Position(2)=ax6.Position(2)+ax6.Position(4);
ax2.Position(2)=ax5.Position(2)+ax5.Position(4);
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
% set(ax3,'Yticklabel',[])

% set(ax6,'Yticklabel',[])
