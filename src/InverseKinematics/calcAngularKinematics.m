function [segment_angles, joint_angles, angular_velocity, angular_acceleration] = calcAngularKinematics(dynamic_lcs, static_lcs, nof, time)
%

% Calculate segment angles
% ------------------------
% Extract segment names
segment_names = fieldnames(dynamic_lcs);

% Loop over segments
for i = 1:length(segment_names)
    % Loop over frames
    for frame = 1:nof
        % Define rotation matrix
        R = [dynamic_lcs.(segment_names{i}).epx(frame,:)' dynamic_lcs.(segment_names{i}).epy(frame,:)' ...
            dynamic_lcs.(segment_names{i}).epz(frame,:)'];
        [x_angle, y_angle, z_angle] = EulerAngles(R,'xyz','deg');

        % Assign output
        segment_angles.(segment_names{i})(frame,:) = [x_angle y_angle z_angle];
    end
end

% Calculate joint angles
% ----------------------
% Define joint names
joint_names = {'hip','knee','ankle'};
side = {'right'};

for i = 1:length(joint_names)
    for j = 1:length(side)
        if isequal(side{j},'right')
            post_script = '_r';
        elseif isequal(side{j},'left')
            post_script = '_l';
        else
            error(['Invalid side: ' side{j}]);
        end
        for frame = 1:nof
            % Define proximal and distal coordinate systems
            switch joint_names{i}
                case 'hip'
                    R_prox = [dynamic_lcs.pelvis.epx(frame,:)' dynamic_lcs.pelvis.epy(frame,:)' ...
                        dynamic_lcs.pelvis.epz(frame,:)'];
                    R_prox_static = [static_lcs.pelvis.epx' static_lcs.pelvis.epy' ...
                        static_lcs.pelvis.epz'];
                    R_dist = [dynamic_lcs.(['thigh' post_script]).epx(frame,:)' dynamic_lcs.(['thigh' post_script]).epy(frame,:)' ...
                        dynamic_lcs.(['thigh' post_script]).epz(frame,:)'];
                    R_dist_static = [static_lcs.(['thigh' post_script]).epx' static_lcs.(['thigh' post_script]).epy' ...
                        static_lcs.(['thigh' post_script]).epz'];
                case 'knee'
                    R_prox = [dynamic_lcs.(['thigh' post_script]).epx(frame,:)' dynamic_lcs.(['thigh' post_script]).epy(frame,:)' ...
                        dynamic_lcs.(['thigh' post_script]).epz(frame,:)'];
                    R_prox_static = [static_lcs.(['thigh' post_script]).epx' static_lcs.(['thigh' post_script]).epy' ...
                        static_lcs.(['thigh' post_script]).epz'];
                    R_dist = [dynamic_lcs.(['leg' post_script]).epx(frame,:)' dynamic_lcs.(['leg' post_script]).epy(frame,:)' ...
                        dynamic_lcs.(['leg' post_script]).epz(frame,:)'];
                    R_dist_static = [static_lcs.(['leg' post_script]).epx' static_lcs.(['leg' post_script]).epy' ...
                        static_lcs.(['leg' post_script]).epz'];
                case 'ankle'
                    R_prox = [dynamic_lcs.(['leg' post_script]).epx(frame,:)' dynamic_lcs.(['leg' post_script]).epy(frame,:)' ...
                        dynamic_lcs.(['leg' post_script]).epz(frame,:)'];
                    R_prox_static = [static_lcs.(['leg' post_script]).epx' static_lcs.(['leg' post_script]).epy' ...
                        static_lcs.(['leg' post_script]).epz'];
                    R_dist = [dynamic_lcs.(['foot' post_script]).epx(frame,:)' dynamic_lcs.(['foot' post_script]).epy(frame,:)' ...
                        dynamic_lcs.(['foot' post_script]).epz(frame,:)'];
                    R_dist_static = [static_lcs.(['foot' post_script]).epx' static_lcs.(['foot' post_script]).epy' ...
                        static_lcs.(['foot' post_script]).epz'];
                case 'midfoot'
                    R_prox = [dynamic_lcs.(['rearfoot' post_script]).epx(frame,:)' dynamic_lcs.(['rearfoot' post_script]).epy(frame,:)' ...
                        dynamic_lcs.(['rearfoot' post_script]).epz(frame,:)'];
                    R_prox_static = [static_lcs.(['rearfoot' post_script]).epx' static_lcs.(['rearfoot' post_script]).epy' ...
                        static_lcs.(['rearfoot' post_script]).epz'];
                    R_dist = [dynamic_lcs.(['forefoot' post_script]).epx(frame,:)' dynamic_lcs.(['forefoot' post_script]).epy(frame,:)' ...
                        dynamic_lcs.(['forefoot' post_script]).epz(frame,:)'];
                    R_dist_static = [static_lcs.(['forefoot' post_script]).epx' static_lcs.(['forefoot' post_script]).epy' ...
                        static_lcs.(['forefoot' post_script]).epz'];
                otherwise
                    error(['Invalid joint: ' joint_names{i}]);
            end

            % Calculate joint angles
            R = (R_prox_static'*R_prox)'*(R_dist_static'*R_dist);
            [x_angle(frame,:), y_angle(frame,:), z_angle(frame,:)] = EulerAngles(R,'xyz','deg');
        end

        % Define output
        joint_angles.(joint_names{i}).(side{j}) = [x_angle y_angle z_angle];
    end
end

% Calculate angular velocities and accelerations
% ----------------------------------------------
% Segment angular velocities
for i = 1:length(segment_names)
    % Find Euler rates
    euler_rate = FiniteDiff(segment_angles.(segment_names{i}), (time(2)-time(1)), 1);

    % Calculate angular velocities
    for frame = 1:nof
        % Define trasformation matrix
        E = [cosd(segment_angles.(segment_names{i})(frame,2))*cosd(segment_angles.(segment_names{i})(frame,3)) -sind(segment_angles.(segment_names{i})(frame,3)) 0;...
            cosd(segment_angles.(segment_names{i})(frame,2))*sind(segment_angles.(segment_names{i})(frame,3)) cosd(segment_angles.(segment_names{i})(frame,3)) 0;...
            sind(segment_angles.(segment_names{i})(frame,2)) 0 1];
        angular_velocity.(segment_names{i})(frame,:) = (E*euler_rate(frame,:)')';
    end

    % Calculate angular acceleration
    angular_acceleration.(segment_names{i}) = FiniteDiff(angular_velocity.(segment_names{i}), time(2)-time(1), 1);
end

% Joint angular velocities
for i = 1:length(joint_names)
    for j = 1:length(side)
        % Find Euler rates
        euler_rate = FiniteDiff(joint_angles.(joint_names{i}).(side{j}), time(2)-time(1), 1);

        % Calculate angular velocities
        for frame = 1:nof
            E = [1 0 -sind(joint_angles.(joint_names{i}).(side{j})(frame,2)); ...
                0 cosd(joint_angles.(joint_names{i}).(side{j})(frame,1)) -sind(joint_angles.(joint_names{i}).(side{j})(frame,1))*cosd(joint_angles.(joint_names{i}).(side{j})(frame,2)); ...
                0 sind(joint_angles.(joint_names{i}).(side{j})(frame,1)) cosd(joint_angles.(joint_names{i}).(side{j})(frame,1))*cosd(joint_angles.(joint_names{i}).(side{j})(frame,2))];
            angular_velocity.(joint_names{i}).(side{j})(frame,:) = (E*euler_rate(frame,:)')';
        end

        % Calculate angular acceleration
        angular_acceleration.(joint_names{i}).(side{j}) = FiniteDiff(angular_velocity.(joint_names{i}).(side{j}), time(2)-time(1), 1);
    end
end
end
