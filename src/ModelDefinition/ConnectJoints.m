function joints = ConnectJoints(jc, root_name)

% Get joint names
tmp = fieldnames(jc);
joint_segment_tbl = [extractBefore(tmp, '_') extractAfter(tmp, '_')];
rm_idx = strcmp(joint_segment_tbl(:, 2), 'global');
joint_segment_tbl(rm_idx, :) = [];

% Initialize variables
n_joints = length(unique(joint_segment_tbl(:, 1)));
joint_name = [];

% Loop over unique joint labels
for i = 1:n_joints
    % Get the name of the joint connected to the current root segment
    tmp = ...
        joint_segment_tbl(strcmp(joint_segment_tbl(:, 2), root_name), 1);
    if isempty(joint_name)
        joint_name = tmp;
    else
        joint_name(i) = tmp(~strcmp(tmp, joint_name{i-1}));
    end

    % Set the parent frame as the current root frame
    S.parent_frame = root_name;

    % Extract the frames connected by the joint and set the not which is
    % not the current root frame as the child frame
    segment_names = ...
        joint_segment_tbl(strcmp(joint_segment_tbl(:, 1), joint_name{i}), 2);
    S.child_frame = ...
        segment_names{~strcmp(segment_names, root_name)};
    
    % Add to output structure
    joints.(joint_name{i}) = S;
    
    % Set the child frame as current root segment
    root_name = joints.(joint_name{i}).child_frame;
end
end