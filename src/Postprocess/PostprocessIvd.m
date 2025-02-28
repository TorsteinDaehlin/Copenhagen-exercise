function R = PostprocessIvd(kinematics, grf, njm, subj, roi, meta)

% Define span to extract data from
tspan = 1.0;
h = tspan * meta.fs;

% Define dimensions
dims = {'x', 'y', 'z'};

% Crate id varaible
R.id = subj.id;

% Loop over ROIs and extract variables of interest
n_roi = size(roi, 1);
for i = 1:n_roi
    % Get region to extract data from
    begin = roi(i,1) + round(((roi(i,2) - roi(i,1))/2) - h/2);
    idx = begin:begin+h;

    % Resultant moment
    R.(['hip_njm_resultant_local_' num2str(i)]) = mean(vecnorm(njm.hip_r(idx, :) ./ subj.mass, 2, 2));
    R.(['knee_njm_resultant_local_' num2str(i)]) = mean(vecnorm(njm.knee_r(idx, :) ./ subj.mass, 2, 2));
    R.(['ankle_njm_resultant_local_' num2str(i)]) = mean(vecnorm(njm.ankle_r(idx, :) ./ subj.mass, 2, 2));

    for j = 1:length(dims)
        % Extract NJM
        R.(['hip_njm_' dims{j} '_' num2str(i)]) = mean(njm.hip_r(idx, j) ./ subj.mass);
        R.(['knee_njm_' dims{j} '_' num2str(i)]) = mean(njm.knee_r(idx, j) ./ subj.mass);
        R.(['ankle_njm_' dims{j} '_' num2str(i)]) = mean(njm.ankle_r(idx, j) ./ subj.mass);

        % Extract kinematics
        R.(['pelvis_ang_' dims{j} '_' num2str(i)]) = mean(kinematics.segment_angles.pelvis(idx, j));
        R.(['thigh_ang_' dims{j} '_' num2str(i)]) = mean(kinematics.segment_angles.thigh_r(idx, j));
        R.(['leg_ang_' dims{j} '_' num2str(i)]) = mean(kinematics.segment_angles.leg_r(idx, j));
        R.(['foot_ang_' dims{j} '_' num2str(i)]) = mean(kinematics.segment_angles.v_foot_r(idx, j));

        R.(['hip_ang_' dims{j} '_' num2str(i)]) = mean(kinematics.joint_angles.hip.right(idx, j));
        R.(['knee_ang_' dims{j} '_' num2str(i)]) = mean(kinematics.joint_angles.knee.right(idx, j));
        R.(['ankle_ang_' dims{j} '_' num2str(i)]) = mean(kinematics.joint_angles.ankle.right(idx, j));

        % Extract GRF
        R.(['GRF_' dims{j} '_' num2str(i)]) = mean(grf.force(idx, j) ./ (subj.mass * 9.81));
    end
end

% Rearrange fields in alphabetical order
R = orderfields(R);

end