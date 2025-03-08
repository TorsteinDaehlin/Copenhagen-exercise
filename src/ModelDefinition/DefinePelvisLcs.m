function [pelvis_lcs, jc] = DefinePelvisLcs(markers, nof)

% Preallocate
origin = zeros(nof,3);
epx = zeros(nof,3);
epy = zeros(nof,3);
epz = zeros(nof,3);

% Define temporary vectors
temp1 = markers.pelvis_RASIS - markers.pelvis_LASIS;
temp2 = 0.5*(markers.pelvis_RASIS + markers.pelvis_LASIS) - ...
    0.5*(markers.pelvis_RPSIS + markers.pelvis_LPSIS);

% Define origin and local x-axis
for i = 1:nof
    origin(i,:) = 0.5*(markers.pelvis_RASIS(i,:) + markers.pelvis_LASIS(i,:));
    epx(i,:) = temp1(i,:)/norm(temp1(i,:));
end

% Find the proportion of temp2 projected onto the local x-axis
proj = dot(temp2',epx');

% Define local y- and z-axes
for i = 1:nof
    temp3 = temp2(i,:) - proj(i)*epx(i,:);
    epy(i,:) = temp3/norm(temp3);
    temp4 = cross(epx(i,:),epy(i,:));
    epz(i,:) = temp4/norm(temp4);
end

% Find hip joint centre
jc.hip_pelvis = FindHipCentre(markers, nof);

% Transform hip centres to globla system
R = [epx', epy', epz'];
jc.hip_global.right = jc.hip_pelvis.right*R' + origin;
jc.hip_global.left = jc.hip_pelvis.left*R' + origin;

% Assing output
pelvis_lcs.origin = origin;
pelvis_lcs.epx = epx;
pelvis_lcs.epy = epy;
pelvis_lcs.epz = epz;
end