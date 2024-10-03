function [dynamic_lcs, jc] = PoseEstimation(dynamic_markers, static_markers, static_lcs, segment_names, static_jc, nof)

% Get the names of the dynamic markers
marker_names = fieldnames(dynamic_markers);

for i = 1:length(segment_names)
    % Get the names of the markers attached to the segments
    segment_markers = marker_names(contains(marker_names, segment_names(i)));
    
    % Perform pose estimation for segment
    dynamic_lcs.(segment_names{i}) = ...
        DynamicLcs(dynamic_markers, static_markers, segment_markers, ...
        static_lcs.(segment_names{i}), nof);
end

% Transform joint centres
joint_names = fieldnames(static_jc);
joint_names = joint_names(~contains(joint_names, {'global'}));

for i = 1:length(joint_names)
    side = fieldnames(static_jc.(joint_names{i}));
    for j = 1:length(side)
        [~, segment_name] = strtok(joint_names{i}, '_');
        segment_name = segment_name(2:end);
        jc.(joint_names{i}) = ...
            TransformJointCentre(dynamic_lcs.(segment_name), static_jc.(joint_names{i}).(side{j}), nof);
    end
end
end


