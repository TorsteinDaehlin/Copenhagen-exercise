function PlotDynamic(markers, kinematics, grf, roi, subj, dynamic_nr)

% Parse input
lcs = kinematics.dynamic_lcs;
jc = kinematics.jc;

% Plot parameters
axis_scale = 0.1;
line_colors = {'r-','g-','b-'};
global_ax = eye(3);

% Define figure
fig = figure('Name',[subj.id]);
fig.WindowState = 'maximized';

% Loop over regions of interest
roi_mean = round(mean(roi, 2));

for k = 1:numel(roi_mean)
    % Loop over frames and animate motion

    % Plot global coordinate system
    for j = 1:3
        plot3([0 global_ax(1,j)*axis_scale*2],[0 global_ax(2,j)*axis_scale*2],[0 global_ax(3,j)*axis_scale*2],line_colors{j});
        hold on;
    end

    % Label global coordinate axes
    text(0.21,0,0,'X','Color','r','FontWeight','bold');
    text(0,0.21,0,'Y','Color','g','FontWeight','bold');
    text(0,0,0.21,'Z','Color','b','FontWeight','bold');

    % Plot markers
    marker_names = fieldnames(markers);
    for i = 1:length(marker_names)
        % Select marker color
        if contains(marker_names{i}, {'_r_'})
            marker_color = '#77AC30';
        elseif contains(marker_names{i}, {'_l_'})
            marker_color = '#0072BD';
        elseif contains(marker_names{i}, {'external'})
            marker_color = '#ffff3d';
        else
            marker_color = '#D95319';
        end

        % Plot marker
        plot3(markers.(marker_names{i})(roi_mean(k),1),markers.(marker_names{i})(roi_mean(k),2),markers.(marker_names{i})(roi_mean(k),3), ...
            'o','MarkerEdgeColor',marker_color,'MarkerFaceColor',marker_color);
    end

    % Plot local coordinate systems
    lcs_names = fieldnames(lcs);
    axis_names = {'epx','epy','epz'};
    for i = 1:length(lcs_names)
        for j = 1:length(axis_names)
            plot3([lcs.(lcs_names{i}).origin(roi_mean(k),1) lcs.(lcs_names{i}).origin(roi_mean(k),1)+lcs.(lcs_names{i}).(axis_names{j})(roi_mean(k),1)*axis_scale], ...
                [lcs.(lcs_names{i}).origin(roi_mean(k),2) lcs.(lcs_names{i}).origin(roi_mean(k),2)+lcs.(lcs_names{i}).(axis_names{j})(roi_mean(k),2)*axis_scale], ...
                [lcs.(lcs_names{i}).origin(roi_mean(k),3) lcs.(lcs_names{i}).origin(roi_mean(k),3)+lcs.(lcs_names{i}).(axis_names{j})(roi_mean(k),3)*axis_scale], ...
                line_colors{j});
        end
    end

    % Plot joint centres
    jc_names = {'hip_thigh_r','knee_thigh_r','ankle_leg_r'};
    for i = 1:length(jc_names)
        plot3(jc.(jc_names{i})(roi_mean(k),1), jc.(jc_names{i})(roi_mean(k),2), jc.(jc_names{i})(roi_mean(k),3), ...
            'o','MarkerEdgeColor','#A2142F','MarkerFaceColor','#A2142F','MarkerSize',10);
    end

    % Plot force platforms and ground reaction force
    for j = 1:length(grf)
        plot3([grf.cop(roi_mean(k),1) grf.cop(roi_mean(k),1)+grf.force(roi_mean(k),1)/1000], ...
            [grf.cop(roi_mean(k),2) grf.cop(roi_mean(k),2)+grf.force(roi_mean(k),2)/1000], ...
            [grf.cop(roi_mean(k),3) grf.cop(roi_mean(k),3)+grf.force(roi_mean(k),3)/1000],'r-');
    end

    % Format axes
    anchor = [0, 0, 0];
    ax = gca;
    ax.Color = 'k';
    ax.XLim = [anchor(1)-1 anchor(1)+1];
    ax.YLim = [anchor(2)-1 anchor(2)+1];
    ax.ZLim = [-0.1 1.2];
    ax.DataAspectRatio = [1 1 1];
    ax.View = [-105 15];
    ax.XLabel.String = 'Position (m)';
    ax.YLabel.String = 'Position (m)';
    ax.ZLabel.String = 'Position (m)';
    
    % Construct trial name
    trial_name = [subj.id ' - Dynamic ' num2str(dynamic_nr) ', repetition ' ...
        num2str(k)];
    title(trial_name);

    % Save figure as data check
    exportgraphics(fig, fullfile(subj.check_path, [subj.id '_dynamic_' num2str(dynamic_nr) ...
        '_' num2str(k) '.jpg']), ...
        'BackgroundColor','current');

    % Clear figure
    pause(2);
    clf(fig);
end

% Close figure
close(fig);

end