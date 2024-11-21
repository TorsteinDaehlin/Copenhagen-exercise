function ExportResults(R, ts, dst_path)

% Get trial names
trial_names = fieldnames(R);

% Loop over trials
for i = 1:length(trial_names)
    % Convert struct to table
    tbls.(trial_names{i}) = struct2table(R.(trial_names{i}));
    tbls.(trial_names{i}) = movevars(tbls.(trial_names{i}), 'id', 'Before', 1);

    % Write tables to excel spreadsheet
    writetable(tbls(trial_names{i}), fullfile(dst_path, 'CPH_results.xlsx'), ...
        'Sheet', trial_names{i});
end

% Save .mat files
save(fullfile(dst_path, 'CPH_results.mat'), 'tbls');
save(fullfile(dst_path, 'CPH_timeseries.mat'), 'ts');

end