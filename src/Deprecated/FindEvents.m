function events = FindEvents(kinematics, kinetics, jump_type)

if ismember(jump_type,{'AJL','AJR'})

    % Define parameters
    side_names = {'right','left'};
    threshold_grf = 20; % N
    threshold_knee_vel = 30;
    dly = 20;

    % Loop over sides and extract events
    for side = 1:length(side_names)
        idx(1) = find(kinetics.grf.(side_names{side}).force(1:end,3) > threshold_grf, 1);
        idx(2) = find(kinetics.grf.(side_names{side}).force(idx(1)+dly:end,3) < threshold_grf, 1);
        idx(3) = find(kinetics.grf.(side_names{side}).force(idx(1)+dly+idx(2):end,3) > threshold_grf, 1);
        idx(4) = find(kinematics.angular_velocity.knee.(side_names{side})(idx(1)+dly+idx(2)+idx(3):end,1) <= threshold_knee_vel & ...
            kinematics.angular_velocity.knee.(side_names{side})(idx(1)+dly+idx(2)+idx(3):end,1) >= -threshold_knee_vel, 1);
        events.(side_names{side}).start = idx(1);
        events.(side_names{side}).foot_off = idx(1)+dly+idx(2);
        events.(side_names{side}).foot_ic = idx(1)+dly+idx(2)+idx(3);
        events.(side_names{side}).max_knee_flex = idx(1)+dly+idx(2)+idx(3)+idx(4);
    end

elseif isequal(jump_type,'VJ')

    % Define parameters
    side_names = {'right','left'};
    threshold_grf = 10; % N
    threshold_knee_vel = 30;
    threshold_pelvis_vel = -0.1;
    dly = 60;

    % Loop over sides and extract events
    for side = 1:length(side_names)
        idx(1) = find(kinematics.velocity.pelvis(:,3) < threshold_pelvis_vel, 1);
        idx(2) = find(kinetics.grf.(side_names{side}).force(idx(1)+dly:end,3) < threshold_grf, 1);
        idx(3) = find(kinetics.grf.(side_names{side}).force(idx(1)+dly+idx(2):end,3) > threshold_grf, 1);
        idx(4) = find(kinematics.angular_velocity.knee.(side_names{side})(idx(1)+dly+idx(2)+idx(3):end,1) < threshold_knee_vel & ...
            kinematics.angular_velocity.knee.(side_names{side})(idx(1)+dly+idx(2)+idx(3):end,1) > -threshold_knee_vel, 1);
        events.(side_names{side}).start = idx(1);
        events.(side_names{side}).foot_off = idx(1)+dly+idx(2);
        events.(side_names{side}).foot_ic = idx(1)+dly+idx(2)+idx(3);
        events.(side_names{side}).max_knee_flex = idx(1)+dly+idx(2)+idx(3)+idx(4);

    end
else
    error(['Invalid jump type: ' jump_type]);
end
end
