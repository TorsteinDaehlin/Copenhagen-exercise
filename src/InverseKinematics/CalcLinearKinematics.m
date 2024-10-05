function [position, velocity, acceleration] = CalcLinearKinematics(segments, dynamic_lcs, time, nof)

% Find segment centre of mass positions in the global coordiante system
segment_names = fieldnames(segments);

% Loop over segments
for i = 1:length(segment_names)
    % Loop over frames and determine position
    for j = 1:nof
        % Construct local coordinate system
        R = [dynamic_lcs.(segment_names{i}).epx(j,:)' dynamic_lcs.(segment_names{i}).epy(j,:)' ...
            dynamic_lcs.(segment_names{i}).epz(j,:)'];

        % Transform centre of mass position to global system
        position.(segment_names{i})(j,:) = (R*segments.(segment_names{i}).com' + dynamic_lcs.(segment_names{i}).origin(j,:)')';
    end

    % Find velocity and acceleration
    h = time(2) - time(1);
    velocity.(segment_names{i}) = FiniteDiff(position.(segment_names{i}),h,1);
    acceleration.(segment_names{i}) = FiniteDiff(position.(segment_names{i}),h,2);
end

end