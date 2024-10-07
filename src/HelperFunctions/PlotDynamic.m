function PlotDynamic(markers, kinematics, grf, roi, subj)

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
n_roi = size(roi, 1);

for k = 1:n_roi
    % Loop over frames and animate motion
    for frame = roi(k, 1):roi(k, 2)

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
            plot3(markers.(marker_names{i})(frame,1),markers.(marker_names{i})(frame,2),markers.(marker_names{i})(frame,3), ...
                'o','MarkerEdgeColor',marker_color,'MarkerFaceColor',marker_color);
        end

        % Plot local coordinate systems
        lcs_names = fieldnames(lcs);
        axis_names = {'epx','epy','epz'};
        for i = 1:length(lcs_names)
            for j = 1:length(axis_names)
                plot3([lcs.(lcs_names{i}).origin(frame,1) lcs.(lcs_names{i}).origin(frame,1)+lcs.(lcs_names{i}).(axis_names{j})(frame,1)*axis_scale], ...
                    [lcs.(lcs_names{i}).origin(frame,2) lcs.(lcs_names{i}).origin(frame,2)+lcs.(lcs_names{i}).(axis_names{j})(frame,2)*axis_scale], ...
                    [lcs.(lcs_names{i}).origin(frame,3) lcs.(lcs_names{i}).origin(frame,3)+lcs.(lcs_names{i}).(axis_names{j})(frame,3)*axis_scale], ...
                    line_colors{j});
            end
        end

        % Plot joint centres
        jc_names = {'hip_thigh_r','knee_thigh_r','ankle_leg_r'};
        for i = 1:length(jc_names)
            plot3(jc.(jc_names{i})(frame,1), jc.(jc_names{i})(frame,2), jc.(jc_names{i})(frame,3), ...
                'o','MarkerEdgeColor','#A2142F','MarkerFaceColor','#A2142F','MarkerSize',10);
        end

        % Plot force platforms and ground reaction force
        for j = 1:length(grf)
            plot3([grf.cop(frame,1) grf.cop(frame,1)+grf.force(frame,1)/1000], ...
                [grf.cop(frame,2) grf.cop(frame,2)+grf.force(frame,2)/1000], ...
                [grf.cop(frame,3) grf.cop(frame,3)+grf.force(frame,3)/1000],'r-');
        end

        % Format axes
        anchor = [0, 0, 0];
        ax = gca;
        ax.Color = 'k';
        ax.XLim = [anchor(1)-1 anchor(1)+1];
        ax.YLim = [anchor(2)-1 anchor(2)+1];
        ax.ZLim = [-0.1 1.2];
        ax.DataAspectRatio = [1 1 1];
        ax.View = [-235 20];
        ax.XLabel.String = 'Position (m)';
        ax.YLabel.String = 'Position (m)';
        ax.ZLabel.String = 'Position (m)';

        pause(0.005);
        if ~isequal(frame,roi(k, 2))
            clf(fig);
        end
    end
end
end