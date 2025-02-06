function CollateResults()

% Load data
[in_file, in_dir] = uigetfile('..\..\*.mat', 'Select .mat results file');
load(fullfile(in_dir, in_file), 'tbls');
idx_file = strtok(in_file, '.');
valid_idx = readtable(fullfile(in_dir, [idx_file '.xlsx']), 'Sheet', 'valid_reps');
idx_header = valid_idx.Properties.VariableNames;

tbl_names = fieldnames(tbls);
for i = tbl_names'
    % Get variable names of table
    var_names = tbls.(i{:}).Properties.VariableNames;
    unique_names = var_names(~contains(var_names, {'id'}));
    unique_names = unique(regexprep(unique_names, '.$', '', 'lineanchors'));

    % Get index matrix for current condition
    cols = contains(idx_header, i);
    isValid = logical(table2array(valid_idx(:, cols)));

    for j = unique_names
        mat = table2array(tbls.(i{:})(:, contains(var_names, j)));
        for k = 1:size(mat, 1)
            s_means.(j{:})(k, 1) = mean(mat(k, isValid(k, :)), 2);
            s_sds.(j{:})(k, 1) = std(mat(k, isValid(k, :)), [], 2);
        end
    end

    % Export data
    T_means = [tbls.(i{:})(:, 'id') struct2table(s_means)];
    T_sds = [tbls.(i{:})(:, 'id') struct2table(s_sds)];

    writetable(T_means, fullfile(in_dir, 'CPH_results_collated.xlsx'), 'Sheet', [i{:} '_means']);
    writetable(T_sds, fullfile(in_dir, 'CPH_results_collated.xlsx'), 'Sheet', [i{:} '_SDs']);
end
end