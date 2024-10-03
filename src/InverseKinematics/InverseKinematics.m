function [kinematics, time] = InverseKinematics(dynamic_markers, static_markers, static_lcs, static_jc, segments, meta)

% Get temporal characteristics
nof = meta.nof;
time = (1:nof)'./data_struct.FrameRate;

% Perform pose estimation
segment_names = fieldnames(segments);
[dynamic_lcs, jc] = PoseEstimation(dynamic_markers, static_markers, ...
    static_lcs, segment_names, static_jc, nof);

% Compute linear kinematics
[position, velocity, acceleration] = ...
    CalcLinearKinematics(segments, dynamic_lcs, time, nof);

% Compute angular kinematics


end