function [kinematics, kinetics, time] = ProcessDynamic(static_markers, static_lcs, static_jc, segments, data_struct, filter_parameters, trial_name, visit_name)

% Get temporal characteristics
nof = data_struct.Frames;
time = (1:nof)'./data_struct.FrameRate;

% Extract marker and force data
dynamic_markers = DefineMarkers(data_struct);
grf = ForceProcess(data_struct, filter_parameters);

% Filter marker data
filter_parameters.fs = data_struct.FrameRate;
dynamic_markers = FilterData(dynamic_markers, filter_parameters, 'markers');

% Compute dynamic poses
[dynamic_lcs, dynamic_jc] = PoseEstimation(dynamic_markers, static_markers, static_lcs, static_jc, nof);

% Calculate linear kinematics
[position, velocity, acceleration] = FindLinearKinematics(segments, dynamic_lcs, time, nof);

% Calculate angular kinematics
[segment_angles, joint_angles, angular_velocity, angular_acceleration] = ...
    FindAngularKinematics(dynamic_lcs, static_lcs, nof, time);

% Calculate angular kinetics
njm = CalculatejointMoments(segments, dynamic_lcs, dynamic_jc, grf, angular_velocity, angular_acceleration, acceleration, position, nof);

% Calculate joint power
joint_power = CalculateJointPower(njm, angular_velocity, dynamic_lcs, nof);

% Plot dynamic trial
% PlotDynamic(dynamic_markers, dynamic_lcs, dynamic_jc, grf, trial_name, visit_name, nof);

% Find normalized centre of pressure
relative_cop = FindRelativeCop(dynamic_markers, dynamic_lcs, segments, grf, nof);

% Assign output structures
kinematics = struct('position',position,'velocity',velocity,'acceleration',acceleration,'segment_angles',segment_angles, ...
    'joint_angles',joint_angles,'angular_velocity',angular_velocity,'angular_acceleration',angular_acceleration);
kinetics = struct('njm',njm,'joint_power',joint_power,'grf',grf,'relative_cop',relative_cop);

end