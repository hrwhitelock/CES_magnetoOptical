function selection_slider()
    % Load HDF5 data
    filename = 'spectroscopy_fgr_all_2025Aug07.h5';
    info = h5info(filename, '/');
    my_spec_data = struct; 
    for i = 1:length(info.Datasets)
        datasetName = info.Datasets(i).Name;
        datasetData = h5read(filename, ['/', datasetName]);
        my_spec_data.(datasetName) = datasetData;
    end

    ampC = my_spec_data.ampC;
    arrC = my_spec_data.linesC;
    calc_field = my_spec_data.calc_field;

    wavenums = linspace(0, 120, length(ampC));

    % Create figure and UI components
    f = figure('Name','Spectroscopy Mask GUI','Position',[100 100 800 600]);

    ax = axes('Parent', f, 'Position', [0.1 0.3 0.85 0.65]);
    xlabel(ax, 'Field (T)');
    ylabel(ax, 'Energy (cm^{-1})');

    % Log range settings (log10)
    logMin = -10;  % 10^-6
    logMax = 1;  % 10^-2
    initLog = -4; % initial = 10^-4

    % Slider for threshold (linear between logMin and logMax)
    uicontrol('Style','text','Parent',f,'String','Mask Threshold (log)',...
              'Units','normalized','Position',[0.1 0.15 0.2 0.05]);
    slider = uicontrol('Style','slider','Parent',f,...
              'Units','normalized','Position',[0.1 0.1 0.7 0.05],...
              'Min', logMin, 'Max', logMax, 'Value', initLog,...
              'Callback', @(src,~) updatePlot(10.^src.Value));

    % Display current threshold value
    threshLabel = uicontrol('Style','text','Parent',f,'Units','normalized',...
                            'Position',[0.82 0.1 0.1 0.05],...
                            'String',num2str(slider.Value,'%.1e'));

    % Initial plot
    updatePlot(slider.Value);

    % --- Nested function to update plot based on threshold ---
    function updatePlot(thresh)
        cla(ax);
        hold(ax,'on');
        n_lines = size(arrC, 1);

        for i = 2:55
            mask = ampC(:, i) > thresh;
            if any(mask)
                if i < 17
                    plot(ax, calc_field(mask), arrC(mask, i)*8.022, 'b-', 'LineWidth', 1.5);
                elseif i >= 17 && i < 32
                    plot(ax, calc_field(mask), arrC(mask, i)*8.022, 'b--', 'LineWidth', 1.5);
                elseif i >= 31 && i < 44
                    plot(ax, calc_field(mask), arrC(mask, i)*8.022, 'b:', 'LineWidth', 1.5);
                elseif i >= 43 && i < 54
                    plot(ax, calc_field(mask), arrC(mask, i)*8.022, 'b-.', 'LineWidth', 1.5);
                end
            end
        end

        title(ax, ['H||c, Mask Threshold = ' num2str(thresh,'%.1e')]);
        hold(ax,'off');
        set(threshLabel,'String',num2str(thresh,'%.1e'));
    end
end