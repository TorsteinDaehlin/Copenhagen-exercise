function CollateDisaggregated()

% Load participant characteristics and get index separating participants by
% sex
[pchar_file, pchar_folder] = uigetfile('..\*xlsx', 'Select participant characteristics');
if ~isequal(pchar_folder, 0) || ~isequal(pchar_file, 0)
    pchar = readtable(fullfile(pchar_folder, pchar_file));
else
    warning('Please select participant characteristics file');
end
fprintf('Data loaded successfully from %s\n', fullfile(pchar_folder, pchar_file));
pchar.Gender = categorical(pchar.Gender, {'W', 'M'}, {'F', 'M'});
pchar = rmmissing(pchar);
idx = ismember(pchar.Gender, 'F');

% Get data path
[data_file, data_folder] = uigetfile('..\*.xlsx', 'Select data');
if isequal(data_folder, 0) || isequal(data_file, 0)
    warning('Please select data path');
end

% Define condition names
trials = {'C_A_means', 'Knee'; 'C_B_means', 'Leg'; 'C_C_means', 'Ankle'};

for i = 1:size(trials, 1)
    % Load data
    data = readtable(fullfile(data_folder, data_file), 'Sheet', trials{i, 1});
    fprintf('Data loaded successfully from %s, sheet %s\n', ...
        fullfile(data_folder, data_file), trials{i, 1});
    
    % Extract outcomes of interest
    if i == 1
        F_data = [table(trials(i, 2), 'VariableNames', {'Support_location'}) ...
            mean(data(idx, 2:end))];
        M_data = [table(trials(i, 2), 'VariableNames', {'Support_location'}) ...
            mean(data(~idx, 2:end))];
    else
        F_data = [F_data; table(trials(i, 2), 'VariableNames', {'Support_location'}) ...
            mean(data(idx, 2:end))];
        M_data = [M_data; table(trials(i, 2), 'VariableNames', {'Support_location'}) ...
            mean(data(~idx, 2:end))];
    end
end

% Write output
[out_file, out_folder] = uiputfile('*.xlsx', 'Save output as');
writetable(F_data, fullfile(out_folder, out_file), 'Sheet', 'Females');
writetable(F_data, fullfile(out_folder, out_file), 'Sheet', 'Males');

end