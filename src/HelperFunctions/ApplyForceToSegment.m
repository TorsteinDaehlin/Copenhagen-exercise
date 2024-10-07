function segment_name = ApplyForceToSegment(jc, cop, roi)

% Determine number of ROIs
n_roi = size(roi, 1);

% Determine cop position relative to joint centres in each ROI
for i = 1:n_roi

    for frame = roi(i, 1):roi(i, 2)
        % Calculate distances
        l_hip_cop = norm(cop(frame, :) - jc.hip_thigh_r(frame, :));
        l_hip_knee = norm(jc.knee_thigh_r(frame, :) - jc.hip_thigh_r(frame, :));

        % Check first condition
        if l_hip_cop < l_hip_knee
            segment_name{i} = 'thigh_r';
        else
            segment_name{i} = 'leg_r';
        end
    end
end
end