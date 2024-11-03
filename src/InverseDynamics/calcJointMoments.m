function njm = calcJointMoments(segments, kinematics, grf)

% Parse input arguments
dynamic_lcs = kinematics.dynamic_lcs;
jc = kinematics.jc;
position = kinematics.position;
acceleration = kinematics.acceleration;
angular_velocity = kinematics.angular_velocity;
angular_acceleration = kinematics.angular_acceleration;

% Get number of frames
nof = size(grf.force, 1);

% Define gravity vector
g = [0 0 -9.81];

% Extract segment and joint names
joint_names = fieldnames(jc);
joint_names(contains(joint_names, {'hip_pelvis'})) = [];

% Loop over each frame
for frame = 1:nof
    for i = 1:length(joint_names)
        % Find names of segments located distal to the joint
        segment_names = getDistalSegments(joint_names{i});

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
        % TODO:
        % ENSURE THAT THIS IS APPLIED TO THE CORRECT SEGMENT
        r_grf = grf.cop(frame,:) - jc.(joint_names{i})(frame,:);

        % Calculate net joint moment
        njm_global = sum((tau_I + tau),1) - ...
            grf.free_moment(frame,:) - cross(r_grf,grf.force(frame,:));

        % Transform net joint moment into the coordiante system of the
        % distal segment
        R = [dynamic_lcs.(segment_names{1}).epx(frame,:)' dynamic_lcs.(segment_names{1}).epy(frame,:)' ...
            dynamic_lcs.(segment_names{1}).epz(frame,:)'];
        njm.(joint_names{i})(frame,:) = (R*njm_global')';
    end
end
end

function segment_names = getDistalSegments(joint_name)

% Get segments distal to the input joint
% TODO: REFACTOR THIS TO AVOID HARDCODING SEGMENT AND JOINT NAMES. A MORE
% GENERIC SCHEME USING PARENT/CHILD OR ROOT/LINK_I MIGHT BE VIABLE OPTIONS
switch joint_name
    case 'hip_thigh_r'
        segment_names = {'thigh_r','leg_r','foot_r'};
    case 'knee_thigh_r'
        segment_names = {'leg_r','foot_r'};
    case 'ankle_leg_r'
        segment_names = {'foot_r'};
    otherwise
        error(['Invalid joint name:' joint_name]);
end

end