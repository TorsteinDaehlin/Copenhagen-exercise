function njm = calcJointMoments(segments, kinematics, grf, joints, grf_act_on, nof)

% Parse input arguments
dynamic_lcs = kinematics.dynamic_lcs;
jc = kinematics.jc;
position = kinematics.position;
acceleration = kinematics.acceleration;
angular_velocity = kinematics.angular_velocity;
angular_acceleration = kinematics.angular_acceleration;

% Define gravity vector
side = 'r';
g = [0 0 -9.81];

% Extract segment and joint names
joint_names = fieldnames(joints);

% Loop over each joint
for i = 1:length(joint_names)

    % Find names of segments located distal to the joint
    segment_names = getDistalSegments(joints, joint_names{i}, side);
   
    for frame = 1:nof
        % Preallocate
        tau_I = zeros(length(segment_names),3);
        tau = zeros(length(segment_names),3);

        % Loop over segments distal to the joint
        for j = 1:length(segment_names)

            % Calculate net force acting on each segment
            F = segments.(segment_names{j}).mass .* (acceleration.(segment_names{j})(frame, :) - g);

            % Transform angular velocity and accelereation vectors to the
            % segment coordinate system
            R = [dynamic_lcs.(segment_names{j}).epx(frame,:)' dynamic_lcs.(segment_names{j}).epy(frame,:)' ...
                dynamic_lcs.(segment_names{j}).epz(frame,:)'];
            omega = (R'*deg2rad(angular_velocity.(segment_names{j})(frame,:))')';
            alpha = (R'*deg2rad(angular_acceleration.(segment_names{j})(frame,:))')';


            % Calculate inretial moment of the segment and transform it
            % back to the global coordiante system
            tau_I(j,:) = (R * (segments.(segment_names{j}).tensor*alpha' + ...
                cross(omega',(segments.(segment_names{j}).tensor*omega'))))';

            % Find moment arm of centre of mass
            r = position.(segment_names{j})(frame,:) - jc.(joint_names{i})(frame,:);

            % Calculate segment moment
            tau(j,:) = cross(r,F);
        end

        % Find moment arm of ground reaction force
        
        r_grf = grf.cop(frame,:) - jc.(joint_names{i})(frame,:);

        % Calculate net joint moment
        njm_global = sum((tau_I + tau),1) - ...
            grf.free_moment(frame,:) - cross(r_grf,grf.force(frame,:));

        % Transform net joint moment into the coordiante system of the
        % distal segment
        R = [dynamic_lcs.(segment_names{1}).epx(frame,:)' dynamic_lcs.(segment_names{1}).epy(frame,:)' ...
            dynamic_lcs.(segment_names{1}).epz(frame,:)'];
        njm.(joint_names{i})(frame,:) = (R'*njm_global')';
    end
end
end

function segment_names = getDistalSegments(joints, joint_name)

% Get the names of every joint
lbls = fieldnames(joints);

% Get child frame of current joint
segment_names{1} = joints.(joint_name).child_frame;

% Search tree to for segment whose parent frame is previous segment name
for i = 1:length(lbls)
    
    current_segment = segment_names{end};

    if strcmp(joints.(lbls{i}).parent_frame, current_segment)
        segment_names{end+1} = joints.(lbls{i}).child_frame;
    end
end
end