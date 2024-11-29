function ExportResults(tbls, R, ts, dst_path, subj_list)

% Get trial names
trial_names = fieldnames(R);

% Loop over trials
for i = 1:length(trial_names)
    % Convert struct to table and append
    tbls.(trial_names{i}) = [tbls.(trial_names{i}) ; struct2table(R.(trial_names{i}))];
    tbls.(trial_names{i}) = movevars(tbls.(trial_names{i}), 'id', 'Before', 1);

    % Remove duplicates
    tbls.(trial_names{i}) = unique(tbls.(trial_names{i}), 'stable');
    tbls.(trial_names{i}) = RemoveDuplicates(tbls.(trial_names{i}), subj_list);

    % Sort table
    tbls.(trial_names{i}) = sortrows(tbls.(trial_names{i}), 'id');

    % Write tables to excel spreadsheet
    writetable(tbls.(trial_names{i}), fullfile(dst_path, 'CPH_results_all_reps.xlsx'), ...
        'Sheet', trial_names{i});
end

% Save .mat files
save(fullfile(dst_path, 'CPH_results_all_reps.mat'), 'tbls');
save(fullfile(dst_path, 'CPH_timeseries.mat'), 'ts');

end

function tbl = RemoveDuplicates(tbl, subj_list)

for i = 1:length(subj_list)
    tf = strcmp(subj_list(i), tbl.id);
    if sum(tf) > 1
        % Get duplicate rows
        duplicate_rows = tbl(tf,:);
       
        if sum(all(ismissing(duplicate_rows(:,2:end)))) == 1 % Retain row without missing values
            new_row = rmmissing(duplicate_rows,'DataVariables',2:size(duplicate_rows,2));
        else % Retain newest row
            new_row = duplicate_rows(end,:);
        end

        % Remove old rows and add new
        tbl(tf, :) = [];
        tbl = [tbl; new_row];

        % Recursively call RetainUnique
        tbl = RemoveDuplicates(tbl, subj_list);
    end
end
end