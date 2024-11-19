function PlotNjmChecks(t, njm, grf, roi, subj, dynamic_nr)

% Create colour palette
cmap = parula(9);
cmap = flipud([cmap(1,:); cmap(4, :); cmap(8,:)]);

% Create njm figure
fig = figure('Name', [subj.id ' - trial ' subj.move_name{dynamic_nr}], ...
    'WindowState', 'maximized');

% Get joint names
joint_names = fieldnames(njm);

% Loop over regions of interst
n_roi = size(roi, 1);
for i = 1:n_roi
    for j = (1:n_roi:length(joint_names)*n_roi)-1
        subplot(length(joint_names), n_roi, i+j);
        p = plot(t(roi(i,1):roi(i,2)), njm.(joint_names{(j/n_roi) + 1})(roi(i,1):roi(i,2),:)./subj.mass, ...
            'LineWidth', 1);
        p(1).Color = cmap(1,:);
        p(2).Color = cmap(2,:);
        p(3).Color = cmap(3,:);
        title([replace(joint_names{(j/n_roi) + 1}, '_', ' ') ' - rep ' num2str(i)]);
        xlabel('Time (s)', 'Interpreter', 'latex');
        ylabel('NJM ($\frac{N \cdot m}{kg}$)', 'Interpreter', 'latex');
    end
end
legend(p, {'X', 'Y', 'Z'}, 'Position', [0.45, 0.02, 0.1, 0.02], ...
    'Orientation', 'horizontal');

% Save NJM figure
exportgraphics(fig, ...
    fullfile(subj.check_path, [subj.id '_NJM_' num2str(dynamic_nr) '.jpg']), ...
    'BackgroundColor','white');

% Close
close(fig);

% Create GRF and cop figure
fig = figure('Name', [subj.id ' - trial ' subj.move_name{dynamic_nr}], ...
    'WindowState', 'maximized');

for i = 1:n_roi
    subplot(3, n_roi, i);
    p = plot(t(roi(i,1):roi(i,2)), grf.force(roi(i, 1):roi(i, 2), :)./(subj.mass * 9.81), ...
        'LineWidth', 1);
    p(1).Color = cmap(1,:);
    p(2).Color = cmap(2,:);
    p(3).Color = cmap(3,:);
    title(['GRF - rep ' num2str(i)]);
    xlabel('Time (s)', 'Interpreter', 'latex');
    ylabel('Force (body weights)', 'Interpreter', 'latex');

    subplot(3, n_roi, n_roi+i);
    p = plot(t(roi(i,1):roi(i,2)), grf.free_moment(roi(i, 1):roi(i, 2), :)./(subj.mass * 9.81), ...
        'LineWidth', 1);
    p(1).Color = cmap(1,:);
    p(2).Color = cmap(2,:);
    p(3).Color = cmap(3,:);
    title(['Free moment - rep ' num2str(i)]);
    xlabel('Time (s)', 'Interpreter', 'latex');
    ylabel('Free moment ($\frac{N \cdot m}{kg}$)', 'Interpreter', 'latex');

    subplot(3, n_roi, 2*n_roi+i);
    p = plot(t(roi(i,1):roi(i,2)), grf.cop(roi(i, 1):roi(i, 2), :), ...
        'LineWidth', 1);
    p(1).Color = cmap(1,:);
    p(2).Color = cmap(2,:);
    p(3).Color = cmap(3,:);
    title(['COP - rep ' num2str(i)]);
    ylabel('Position (m)', 'Interpreter', 'latex');
end

legend(p, {'X', 'Y', 'Z'}, 'Position', [0.45, 0.02, 0.1, 0.02], ...
    'Orientation', 'horizontal');

% Save NJM figure
exportgraphics(fig, ...
    fullfile(subj.check_path, [subj.id '_GRF_' num2str(dynamic_nr) '.jpg']), ...
    'BackgroundColor','white');

% Close
close(fig);

end