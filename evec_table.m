%% Load the field values and eigenvectors
H = h5read('moo_3d_result.h5', '/z/T_0/Heff');  % [Nh]
evecs = h5read('moo_3d_result.h5', '/z/T_0/evecs');  % size: [16, 16, Nh]

% Ensure eigenvectors are [spin, state, field]
if size(evecs,1) ~= 16
    evecs = permute(evecs, [2 1 3]);
end

% Spin state labels
Jz_states = (15:-2:-15)';  % Column vector of spin projections

% Find indices of H closest to every 0.5 T step
H_sampled = 0:0.5:9;
[H_indices, unique_idx] = unique(arrayfun(@(hval) ...
    find(abs(H - hval) == min(abs(H - hval)), 1), H_sampled), 'stable');

% Loop over each eigenstate and generate a table
for state_idx = 1:16
    % Preallocate matrix 
    weights_matrix = zeros(16, numel(H_indices));

    for j = 1:numel(H_indices)
        h_idx = H_indices(j);
        vec = evecs.r(:, state_idx, h_idx);        % 16-element complex vector
        weights = abs(vec).^2;                   % Weight of each spin component
        weights_matrix(:, j) = weights;          % Store in column
    end

    % Create and display table
    T = array2table(weights_matrix, ...
        'RowNames', cellstr(string(Jz_states)), ...
        'VariableNames', cellstr("H_" + string(round(H(H_indices),1))));

    disp(['Eigenvector ', num2str(state_idx)]);
    disp(T);

    % Optionally save to file:
    writetable(T, sprintf('Eigenstate_%02d_weights.csv', state_idx), 'WriteRowNames', true);
end

%%
nStates = 16;
eigenTables = cell(nStates, 1);
 
for state_idx = 1:nStates
    filename = sprintf('Eigenstate_%02.0f_weights.csv', state_idx);
    eigenTables{state_idx} = readtable(filename, 'ReadRowNames', true);
end

 