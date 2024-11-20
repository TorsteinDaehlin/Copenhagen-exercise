function  [roi] = IdentifyROI(t, grf, subj, dynamic_nr)

% Truncate force signal
grf_norm = grf ./ (subj.mass * 9.81);

% The vertical ground reaction force is heavily filtered using a moving
% average smoothing approach before identifying when the participant is in
% contact with the force platform. This region is subsequently refined into
% the regions of interest (ROI)
idx = findchangepts(movmean(grf_norm(:, 3), 200), 'MaxNumChanges', 2);
roi = findchangepts(movmean(grf_norm(idx(1):idx(2), 3), 50), 'MaxNumChanges', 6);

% Reshape ideces so the first column contains the start of a region of
% interest and the second constains the end of a region of interst
roi = reshape(idx(1) + roi, 2, [])';

% Prompt user to decide if ROIs are correct
iscorrect = false;
while ~iscorrect
    % Plot current ROis
    fig = PlotROIs(t, grf_norm, roi);

    answer = questdlg('Are regions of interest identified correctly?', ...
        'ROI', 'Yes', 'No', 'No');

    if strcmp(answer, 'Yes')
        iscorrect = true;
    else
        % Close current figure
        close(fig);

        % Plot GRF
        fig2 = figure();
        cmap = parula(9);
        cmap = flipud([cmap(1,:); cmap(4, :); cmap(8,:)]);
        p = plot(grf_norm);
        p(1).Color = cmap(1,:);
        p(2).Color = cmap(2,:);
        p(3).Color = cmap(3,:);
        xlabel('Sample (nr)');
        ylabel('Force (N/kg)');
        legend(p, {'x', 'y', 'z'});
        hold on;
        
        % Get user input
        for i = 1:6
            [x, ~] = ginput(1);
            tmp(i) = round(x);
            if mod(i, 2) == 0
                xline(tmp(i),'k--',['$t_{end, ' num2str(ceil(i/2)) '}$'], ...
                    'Interpreter', 'latex');
            else
                 xline(tmp(i),'k--',['$t_{0, ' num2str(ceil(i/2)) '}$'], ...
                    'Interpreter', 'latex');
            end
        end
        close(fig2);
        
        % Reshape output
        roi = reshape(tmp, 2, [])';
    end
end

% Save to output directory
exportgraphics(fig, ...
    fullfile(subj.check_path, [subj.id '_ROI_' num2str(dynamic_nr) '.jpg']), ...
    'BackgroundColor','white');

% Close figure
close(fig);

end

function fig = PlotROIs(t, grf_norm, roi)

% Create data check figure
cmap = parula(9);
cmap = flipud([cmap(1,:); cmap(4, :); cmap(8,:)]);

fig = figure();
c = size(grf_norm, 2);
for i = 1:c
    p (i) = plot(t, grf_norm(:, i), 'Color', cmap(i, :));
    hold on;
end

for j = 1:size(roi, 1)
    xline(t(roi(j, :)),'k--', {['$t_{0, ' num2str(j) '}$'], ...
        ['$t_{end, ' num2str(j) '}$']}, ...
        'Interpreter', 'latex');
end
ylim([-0.25 1]);
legend([p(1), p(2), p(3)], {'X', 'Y', 'Z'});
xlabel('$Time (s)$', 'Interpreter', 'latex');
ylabel('$GRF (\frac{N}{m \cdot g})$', 'Interpreter', 'latex');
title('$ROI data check$', 'Interpreter', 'latex');

end