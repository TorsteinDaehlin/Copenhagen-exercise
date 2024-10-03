function [lcs_dynamic, jc] = PoseEstimation(dynamic_markers, static_markers, static_lcs, segment_names, static_jc, nof)

% Get the names of the dynamic markers
marker_names = fieldnames(dynamic_markers);

for i = 1:length(segment_names)
    % Get the names of the markers attached to the segments
    segment_markers = marker_names(contains(marker_names, segment_names(i)));
    
    % Perform pose estimation for segment
    lcs_dynamic.(segment_names{i}) = ...
        DynamicLcs(dynamic_markers, static_markers, segment_markers, ...
        static_lcs.(segment_names{i}), nof);
end

% Transform joint centres
joint_names = fieldnames(static_jc);
joint_names = joint_names(~contains(joint_names, {'global'}));

for i = 1:length(joint_names)
    side = fieldnames(static_jc.(joint_names{i}));
    for j = 1:length(side)
        [~, segment_name] = strtok(joint_names{i}, '_');
        segment_name = segment_name(2:end);
        jc.(joint_names{i}) = ...
            TransformJointCentre(lcs_dynamic.(segment_name), static_jc.(joint_names{i}).(side{j}), nof);
    end
end

% % Define dynamic pelvis system
% lcs.pelvis = DynamicPelvis(dynamic_markers, static_markers, static_lcs, nof);
% jc.hip.right = TransformJointCentre(lcs.pelvis, static_jc.hip.right, nof);
% jc.hip.left = TransformJointCentre(lcs.pelvis, static_jc.hip.left, nof);
% 
% % Define dynamic thigh systems
% lcs.thigh_r = DynamicThigh(dynamic_markers, static_markers, static_lcs, 'right', nof);
% lcs.thigh_l = DynamicThigh(dynamic_markers, static_markers, static_lcs, 'left', nof);
% jc.knee.right = TransformJointCentre(lcs.thigh_r, static_jc.knee.right, nof);
% jc.knee.left = TransformJointCentre(lcs.thigh_l, static_jc.knee.left, nof);
% 
% % Define dynamic leg systems
% lcs.leg_r = DynamicLeg(dynamic_markers, static_markers, static_lcs, 'right', nof);
% lcs.leg_l = DynamicLeg(dynamic_markers, static_markers, static_lcs, 'left', nof);
% jc.ankle.right = TransformJointCentre(lcs.leg_r, static_jc.ankle.right, nof);
% jc.ankle.left = TransformJointCentre(lcs.leg_l, static_jc.ankle.left, nof);
% 
% % Define dynamic foot systems
% lcs.foot_r = DynamicFoot(dynamic_markers, static_markers, static_lcs, 'right', nof);
% lcs.foot_l = DynamicFoot(dynamic_markers, static_markers, static_lcs, 'left', nof);
% 
% % Define dynamic rearfoot systems
% lcs.rearfoot_r = DynamicRearfoot(dynamic_markers, static_markers, static_lcs, 'right', nof);
% lcs.rearfoot_l = DynamicRearfoot(dynamic_markers, static_markers, static_lcs, 'left', nof);
% 
% % Define dynamic forefoot systems
% lcs.forefoot_r = DynamicForefoot(dynamic_markers, static_markers, static_lcs, 'right', nof);
% lcs.forefoot_l = DynamicForefoot(dynamic_markers, static_markers, static_lcs, 'left', nof);

end


