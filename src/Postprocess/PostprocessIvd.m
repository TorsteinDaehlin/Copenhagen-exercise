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

    for j = 1:length(dims)
        % Extract NJM
        R.(['hip_njm_' dims{j} '_' num2str(i)]) = mean(njm.hip_r(idx, j) ./ subj.mass);
        R.(['knee_njm_' dims{j} '_' num2str(i)]) = mean(njm.knee_r(idx, j) ./ subj.mass);
        R.(['ankle_njm_' dims{j} '_' num2str(i)]) = mean(njm.ankle_r(idx, j) ./ subj.mass);

        % Extract kinematics
        R.(['hip_ang_' dims{j} '_' num2str(i)]) = mean(kinematics.joint_angles.hip.right(idx, j));
        R.(['knee_ang_' dims{j} '_' num2str(i)]) = mean(kinematics.joint_angles.knee.right(idx, j));
        R.(['ankle_ang_' dims{j} '_' num2str(i)]) = mean(kinematics.joint_angles.ankle.right(idx, j));

        % Extract GRF
        R.(['GRF_' dims{j} '_' num2str(i)]) = mean(grf.force(idx, j) ./ (subj.mass * 9.81));
    end
end

% Average over repetitions
% R = structfun(@(x) mean(x), R, 'UniformOutput', false);

% Rearrange fields in alphabetical order
R = orderfields(R);

end