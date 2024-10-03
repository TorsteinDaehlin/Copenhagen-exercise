function lcs_dynamic = DynamicLcs(dynamic_markers, static_markers, marker_names, static_lcs, nof)

% Preallocate
epx = zeros(nof,3);
epy = zeros(nof,3);
epz = zeros(nof,3);
origin = zeros(nof,3);

% Extract static markers and lcs
neutral = [];
for i = 1:length(marker_names)
    neutral = [neutral mean(static_markers.(marker_names{i}))'];
end
R_static = [static_lcs.epx' static_lcs.epy' static_lcs.epz']; % Define static rotation matrix

% Define static positions in local coordinate system
for i = 1:size(neutral,2)
    neutral(:,i) = R_static'*(neutral(:,i) - static_lcs.origin');
end

% Loop over each frame
for i = 1:nof
    move = [];
    for j = 1:length(marker_names)
        move = [move dynamic_markers.(marker_names{j})(i, :)'];
    end
    
    % Perform pose estimation
    T = LeastSquarePose(neutral,move);
    epx(i,:) = T(1:3,1)';
    epy(i,:) = T(1:3,2)';
    epz(i,:) = T(1:3,3)';
    origin(i,:) = T(1:3,4)';
end

% Assign output
lcs_dynamic.epx = epx;
lcs_dynamic.epy = epy;
lcs_dynamic.epz = epz;
lcs_dynamic.origin = origin;

end
