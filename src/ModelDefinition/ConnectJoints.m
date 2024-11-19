function joints = ConnectJoints(jc, root_name)

% Get joint names
tmp = fieldnames(jc);
joint_segment_tbl = [extractBefore(tmp, '_') extractAfter(tmp, '_')];
rm_idx = strcmp(joint_segment_tbl(:, 2), 'global');
joint_segment_tbl(rm_idx, :) = [];

% Initialize variables
n_joints = length(unique(joint_segment_tbl(:, 1)));
joint_name = [];

% Get links connected to root regment
sides = fieldnames(jc.(tmp{1}));

% Preallocate
joints = struct();

for i = 1:length(sides)
    joints = getChildSegment(joints, joint_segment_tbl, root_name, joint_name, sides{i});
end

% Eliminate any non-existing joints
joint_labels = fieldnames(joints);
for i = 1:length(joint_labels)
    idx = find(contains(joint_segment_tbl(:, 1), strtok(joint_labels{i}, '_')),1);
    side_name = fieldnames(jc.(strjoin(joint_segment_tbl(idx, :), '_')));
    if ~any(contains(side_name, joint_labels{i}(end)))
        joints = rmfield(joints, joint_labels{i});
    end
end
end

function joints = getChildSegment(joints, joint_segment_tbl, root_name, joint_name, side)

% Get the name of the joint connected to the current root segment
tmp = ...
    joint_segment_tbl(strcmp(joint_segment_tbl(:, 2), root_name), 1);

if isempty(joint_name)
    joint_name = tmp;
else
    joint_name = tmp(~strcmp(tmp, joint_name));
end

if ~isempty(joint_name)
    % Set the parent frame as the current root frame
    S.parent_frame = [root_name '_' lower(side(1))];

    % Extract the frames connected by the joint and set the not which is
    % not the current root frame as the child frame
    segment_names = ...
        joint_segment_tbl(strcmp(joint_segment_tbl(:, 1), joint_name), 2);

    % Get child frame
    child_frame = ...
        segment_names{~strcmp(segment_names, root_name)};
    S.child_frame = [child_frame '_' lower(side(1))];

    % Add to output structure
    joints.([joint_name{:} '_' lower(side(1))]) = S;

    % Set the child frame as current root segment
    root_name = child_frame;

    % Recursively call
    joints = getChildSegment(joints, joint_segment_tbl, root_name, joint_name, side);

else
    return;
end
end