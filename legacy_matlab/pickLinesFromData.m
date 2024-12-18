% pick data from plot

%% first, draw surface + picked peaks

figure(); 
hold on; 
[X,Y] = meshgrid(dataStruct.wavenumbers, dataStruct.avgField);
c = contourf(X,Y, (2-dataStruct.avgData0field)*2000, 150, 'Linestyle', 'none');
view(2)
% clim([2500, 3500])
xlim([4 120])
colormap(cool);
title('CsErSe2 zf subtraction')
xlabel('wavenumber (units???)')
ylabel('field')
scatter(dataStruct.wavenumbers, dataStruct.pks2d, 'black')

%% now we want to loop through the lines

numLines = 1;% total number of lines I'm trying to sort from the plot
% for the CES case we're going to jsut go with the number of lines ian fit
wave = dataStruct.pksWave;
field = dataStruct.pksField; 
lines = []; 
for i= 1:numLines % we want to skip the gs because we necessarily do not observe it in IR or raman
    % do picking
    % write this so we stop picking for each line on key press
    disp(i-1) % really hate 1 index
    line = []; 
    mymsg = msgbox('nextline?');
    while ishandle(mymsg)
        [x,y]=ginput; 
    end
    for i = 1:length(x)
        % find which data points we've selected from the peak finding arr
        % this is not the most efficient but whatever
        for j =1:length(wave)
            if abs(x(i)-wave(j))<1 && abs(y(i)-field(j))<.5
                line = [line, [wave(j); field(j)]];
            end
        end
    end
    lines = [lines; line]; 
end