function [static_lcs, static_jc, segments] = ProcessStatic(static, meta, subj)

% Get temporal characteristics
nof = meta.nof;
frame_rate = meta.fs;

% Reduce static markers and define anatomical coordinate systems
markers = static.markers;
marker_names = fieldnames(markers);
for i = 1:length(marker_names)
    markers.(marker_names{i}) = mean(markers.(marker_names{i})((nof/2)-(frame_rate/2):(nof/2)+(frame_rate/2),:));
end

% Define local coordinate systems
[static_lcs, static_jc] = DefineLocalSystems(markers);

% Define segment parameters
segment_names = fieldnames(static_lcs);
segments = DefineSegments(markers, static_lcs, static_jc, subj, segment_names);

% Plot static trial
PlotStatic(markers, static_lcs, static_jc, segments, subj);

end

% Private functions
% =================

% Plot static produces a plot of the static trial to assist the user in
% assessing whether the model is created correctly.
function PlotStatic(markers, lcs, jc, segments, subj)

% Plot parameters
axis_scale = 0.1;
line_colors = {'r-','g-','b-'};

% Define figure
fig = figure('Name',['Static trial - ' subj.id]);
fig.WindowState = 'maximized';

% Plot global coordinate system
global_ax = eye(3);
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
    else
        marker_color = '#D95319';
    end

    % Plot marker
    plot3(markers.(marker_names{i})(1),markers.(marker_names{i})(2),markers.(marker_names{i})(3), ...
        'o','MarkerEdgeColor',marker_color,'MarkerFaceColor',marker_color);
end

% Plot local coordinate systems
lcs_names = fieldnames(lcs);
axis_names = {'epx','epy','epz'};
for i = 1:length(lcs_names)
    for j = 1:length(axis_names)
        plot3([lcs.(lcs_names{i}).origin(1) lcs.(lcs_names{i}).origin(1)+lcs.(lcs_names{i}).(axis_names{j})(1)*axis_scale], ...
            [lcs.(lcs_names{i}).origin(2) lcs.(lcs_names{i}).origin(2)+lcs.(lcs_names{i}).(axis_names{j})(2)*axis_scale], ...
            [lcs.(lcs_names{i}).origin(3) lcs.(lcs_names{i}).origin(3)+lcs.(lcs_names{i}).(axis_names{j})(3)*axis_scale], ...
            line_colors{j});
    end

end

% Plot segment centre of mass
segment_names = fieldnames(segments);
for i = 1:length(segment_names)
    R = [lcs.(segment_names{i}).epx' lcs.(segment_names{i}).epy' lcs.(segment_names{i}).epz'];
    com = R*segments.(segment_names{i}).com' + lcs.(segment_names{i}).origin';
    plot3(com(1),com(2), com(3), 'g*');
end

% Plot joint centres
jc_names = {'hip_global','knee_global','ankle_global'};
sides = {'right'};
for i = 1:length(jc_names)
    for j = 1:length(sides)
        plot3(jc.(jc_names{i}).(sides{j})(1),jc.(jc_names{i}).(sides{j})(2),jc.(jc_names{i}).(sides{j})(3), ...
            'o','MarkerEdgeColor','#A2142F','MarkerFaceColor','#A2142F','MarkerSize',10);
    end
end

% Format axes
anchor = lcs.pelvis.origin;
ax = gca;
ax.Color = 'k';
ax.XLim = [anchor(1)-0.75 anchor(1)+0.75];
ax.YLim = [anchor(2)-0.75 anchor(2)+0.75];
ax.ZLim = [-0.1 1.4];
ax.DataAspectRatio = [1 1 1];
ax.View = [-210 10];
ax.XLabel.String = 'Position (m)';
ax.YLabel.String = 'Position (m)';
ax.ZLabel.String = 'Position (m)';

% Annotate participant characteristics
x_pos = ax.XLim(1)+0.1;
y_pos = ax.YLim(end)-0.1;
z_pos = ax.ZLim(end)-0.1;
text(x_pos,y_pos,z_pos,subj.id,'Color','w','HorizontalAlignment','center','FontWeight','bold');
text(x_pos,y_pos,z_pos-0.05,[num2str(subj.height,3) ' m'],'Color','w','HorizontalAlignment','center');
text(x_pos,y_pos,z_pos-0.10,[num2str(subj.mass,3) ' kg'],'Color','w','HorizontalAlignment','center');

% TODO: SAVE FIGURE AS DATA CHECK

% Display figure before closing
pause(5);
close(fig);

end
