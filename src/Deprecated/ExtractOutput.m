function [output, jump_height] = ExtractOutput(static_lcs, segments, kinematics, kinetics, events, old_jump_height, old_output, time)

% Find jump height
R = [static_lcs.pelvis.epx' static_lcs.pelvis.epy' static_lcs.pelvis.epz'];
height_0 = R*segments.pelvis.com' + static_lcs.pelvis.origin';
height_max = max(kinematics.position.pelvis(:,3));
jump_height = height_max - height_0(3,1);

if jump_height > old_jump_height

    % calculate time step
    h = time(2)-time(1);

    % Find net joint work and average net joint moment
    joint_names = fieldnames(kinetics.joint_power);
    side_names = fieldnames(kinetics.joint_power.(joint_names{1}));
    for side = 1:length(side_names)
        for joint = 1:length(joint_names)
            % Extract propulsion phase variables of interest
            work.(joint_names{joint}).(side_names{side}).propulsion = ...
                trapz(h,kinetics.joint_power.(joint_names{joint}).(side_names{side})(events.(side_names{side}).start:events.(side_names{side}).foot_off,1));
            moment.(joint_names{joint}).(side_names{side}).propulsion = ...
                mean(kinetics.njm.(joint_names{joint}).(side_names{side})(events.(side_names{side}).start:events.(side_names{side}).foot_off,1));
            excursion.(joint_names{joint}).(side_names{side}).propulsion = ...
                kinematics.joint_angles.(joint_names{joint}).(side_names{side})(events.(side_names{side}).foot_off,1) - ...
                kinematics.joint_angles.(joint_names{joint}).(side_names{side})(events.(side_names{side}).start,1);
            arch_deformation.(side_names{side}).propulsion = max(kinematics.joint_angles.midfoot.(side_names{side})(events.(side_names{side}).start:events.(side_names{side}).foot_off,1));

            % Extract landing phase variables of interest
            work.(joint_names{joint}).(side_names{side}).landing = ...
                trapz(h,kinetics.joint_power.(joint_names{joint}).(side_names{side})(events.(side_names{side}).foot_ic:events.(side_names{side}).max_knee_flex,1));
            moment.(joint_names{joint}).(side_names{side}).landing = ...
                mean(kinetics.njm.(joint_names{joint}).(side_names{side})(events.(side_names{side}).foot_ic:events.(side_names{side}).max_knee_flex,1));
            excursion.(joint_names{joint}).(side_names{side}).landing = ...
                kinematics.joint_angles.(joint_names{joint}).(side_names{side})(events.(side_names{side}).max_knee_flex,1) - ...
                kinematics.joint_angles.(joint_names{joint}).(side_names{side})(events.(side_names{side}).foot_ic,1);
            arch_deformation.(side_names{side}).landing = max(kinematics.joint_angles.midfoot.(side_names{side})(events.(side_names{side}).foot_ic:events.(side_names{side}).max_knee_flex,1));


            % Define start time
            if events.right.start < events.left.start
                start_frame = events.right.start;
            else
                start_frame = events.left.start;
            end

            % Define end time
            if events.right.max_knee_flex > events.left.max_knee_flex
                end_frame =  events.right.max_knee_flex;
            else
                end_frame =  events.left.max_knee_flex;
            end

            % Extract time series
            time_series.moment.(joint_names{joint}).(side_names{side}) = kinetics.njm.(joint_names{joint}).(side_names{side})(start_frame:end_frame,:);
            time_series.power.(joint_names{joint}).(side_names{side}) = kinetics.joint_power.(joint_names{joint}).(side_names{side})(start_frame:end_frame,:);
            time_series.angles.(joint_names{joint}).(side_names{side}) = kinematics.joint_angles.(joint_names{joint}).(side_names{side})(start_frame:end_frame,:);

        end
        % Extract centre of pressure location
        pressure.(side_names{side}).propulsion = mean(kinetics.relative_cop.(side_names{side})(events.(side_names{side}).start:events.(side_names{side}).foot_off,:));
        pressure.(side_names{side}).landing = mean(kinetics.relative_cop.(side_names{side})(events.(side_names{side}).foot_ic:events.(side_names{side}).max_knee_flex,:));

        % Extract ground reaction force time series
        time_series.grf.(side_names{side}) = kinetics.grf.(side_names{side}).force;
    end

    % Plot events for verification
    PlotEvents(kinetics, kinematics, events, time);

    % Assign output
    output = struct('jump_height',jump_height,'work',work,'moment',moment,'excursion',excursion, ...
        'arch_deformation',arch_deformation,'pressure',pressure,'time_series',time_series);
else
    output = old_output;
end
end
