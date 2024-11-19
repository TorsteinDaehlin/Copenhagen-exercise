function PlotKinematicsChecks(t, kinematics, roi, subj, dynamic_nr)

% Get kinematics of interest
var_names = fieldnames(kinematics);
var_names(contains(var_names, {'dynamic_lcs', 'jc'})) = [];

for i = 1:length(var_names)
    % Create axis label
    lbl = replace(var_names{i}, '_', ' ');
    lbl = [upper(lbl(1)) lower(lbl(2:end))];

    % Plot
    fig = DrawPlot(t, kinematics.(var_names{i}), roi, subj, dynamic_nr, lbl);

    % Save figure
    exportgraphics(fig, ...
    fullfile(subj.check_path, [subj.id '_kinematics_' var_names{i} '_' num2str(dynamic_nr) '.jpg']), ...
    'BackgroundColor','white');

    % Close
    close(fig);
end


end

function fig = DrawPlot(t, var, roi, subj, dynamic_nr, lbl)

% Create data check figure
cmap = parula(9);
cmap = flipud([cmap(1,:); cmap(4, :); cmap(8,:)]);

% Plot segment COM position
fig = figure('Name', [subj.id ' - trial ' subj.move_name{dynamic_nr}], ...
    'WindowState','maximized');

% Get segment names and numbers of plots
field_names = fieldnames(var);
n_seg = length(field_names);
n_roi = size(roi, 1);

% Loop over plot windows
for i = 1:n_roi
    for j = (1:n_roi:n_seg*n_roi)-1
        subplot(n_seg, n_roi, i+j);
        if isstruct(var.(field_names{(j/n_roi) + 1}))
            p = plot(t(roi(i,1):roi(i,2)), var.(field_names{(j/n_roi) + 1}).right(roi(i,1):roi(i,2),:), 'LineWidth', 1);
        else
            p = plot(t(roi(i,1):roi(i,2)), var.(field_names{(j/n_roi) + 1})(roi(i,1):roi(i,2),:), 'LineWidth', 1);
        end
        p(1).Color = cmap(1, :);
        p(2).Color = cmap(2, :);
        p(3).Color = cmap(3, :);
        title([replace(field_names{(j/n_roi) + 1}, '_', ' ') ' - rep ' num2str(i)]);
        ylabel(lbl, 'Interpreter', 'latex');
        xlabel('Time (s)', 'Interpreter', 'latex');
    end
end
legend(p, {'X', 'Y', 'Z'}, 'Position', [0.45, 0.02, 0.1, 0.02], ...
    'Orientation', 'horizontal');
end