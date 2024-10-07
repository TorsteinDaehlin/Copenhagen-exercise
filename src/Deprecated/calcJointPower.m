function joint_power = calcJointPower(njm, angular_velocity, dynamic_lcs, nof)

% Get joint names and side names
joint_names = fieldnames(njm);
side_names = fieldnames(njm.(joint_names{1}));

for joint = 1:length(joint_names)
    for side = 1:length(side_names)
        % Find segment names (1st element of the vector returned is the segment distal to the joint)
        segment_names = FindSegmentNames([joint_names{joint} '_' side_names{side}]);

        for frame = 1:nof
            % Transform joint angular velocity into local coordinate system of the distal segment (njm is already expressed here)
            R = [dynamic_lcs.(segment_names{1}).epx(frame,:)' dynamic_lcs.(segment_names{1}).epy(frame,:)' ...
                dynamic_lcs.(segment_names{1}).epz(frame,:)'];
            omega = (R'*deg2rad(angular_velocity.(joint_names{joint}).(side_names{side})(frame,:))')';

            % Calculate instantaneous power
            joint_power.(joint_names{joint}).(side_names{side})(frame,:) =  ...
                njm.(joint_names{joint}).(side_names{side})(frame,:).*omega;
        end
    end
end
end
