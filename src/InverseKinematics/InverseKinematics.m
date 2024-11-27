function [kinematics, time] = InverseKinematics(dynamic_markers, static_markers, static_lcs, static_jc, segments, meta)

% Get temporal characteristics
nof = meta.nof;
time = (1:nof)'./ meta.fs;

% Perform pose estimation
segment_names = fieldnames(static_lcs);
[dynamic_lcs, jc] = PoseEstimation(dynamic_markers, static_markers, ...
    static_lcs, segment_names, static_jc, nof);

% Compute linear kinematics
[position, velocity, acceleration] = ...
    CalcLinearKinematics(segments, dynamic_lcs, time, nof);

% Compute angular kinematics
[segment_angles, joint_angles, angular_velocity, angular_acceleration] = ...
    calcAngularKinematics(dynamic_lcs, static_lcs, nof, time);

% Define output structure
kinematics.dynamic_lcs = dynamic_lcs;
kinematics.jc = jc;
kinematics.position = position;
kinematics.velocity = velocity;
kinematics.acceleration = acceleration;
kinematics.segment_angles = segment_angles;
kinematics.joint_angles = joint_angles;
kinematics.angular_velocity = angular_velocity;
kinematics.angular_acceleration = angular_acceleration;
end