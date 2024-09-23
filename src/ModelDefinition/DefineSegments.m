function segments = DefineSegments(markers, lcs, jc, subj, segment_names)

% Loop over segments and define segment parameters
for i = 1:length(segment_names)
    segments.(segment_names{i}) = CreateSegment(markers, lcs, jc, subj, segment_names{i});
end
end

% Private functions
% =================
function segment = CreateSegment(markers, lcs, jc, subj, segment_name)

% Define length, proximal radius, and distal radius for the given segment
switch segment_name
    case 'pelvis'
        segment.len = dot(lcs.pelvis.origin - 0.5*(jc.hip_global.right + jc.hip_global.left), lcs.pelvis.epz);
        segment.prox_rad = norm(markers.pelvis_RASIS - markers.pelvis_LASIS) / 2;
        segment.dist_rad = norm(jc.hip_global.right - jc.hip_global.left) / 2;
        segment.mass = subj.mass*0.142;
        offset = mean([markers.pelvis_RASIS; markers.pelvis_LASIS; ...
                       markers.pelvis_RPSIS; markers.pelvis_LPSIS]) ...
                       - lcs.pelvis.origin;
    case 'thigh_r'
        segment.len = norm(jc.hip_global.right - 0.5*(markers.thigh_r_dist_med + markers.thigh_r_dist_lat));
        segment.prox_rad = norm(jc.hip_global.right - markers.thigh_r_cluster1);
        segment.dist_rad = norm(markers.thigh_r_dist_med - markers.thigh_r_dist_lat) / 2;
        segment.mass = subj.mass*0.100;
    case 'thigh_l'
        segment.len = norm(jc.hip_global.left - 0.5*(markers.thigh_l_dist_med + markers.thigh_l_dist_lat));
        segment.prox_rad = norm(jc.hip_global.left - markers.thigh_l_cluster1);
        segment.dist_rad = norm(markers.thigh_l_dist_med - markers.thigh_l_dist_lat) / 2;
        segment.mass = subj.mass*0.100;
    case 'leg_r'
        segment.len = norm(jc.knee_global.right - jc.ankle_global.right);
        segment.prox_rad = norm(markers.thigh_r_dist_med - markers.thigh_r_dist_lat) / 2;
        segment.dist_rad = norm(markers.leg_r_dist_med - markers.leg_r_dist_lat) / 2;
        segment.mass = subj.mass*0.0465;
    case 'leg_l'
        segment.len = norm(jc.knee_global.left - jc.ankle_global.left);
        segment.prox_rad = norm(markers.thigh_l_dist_med - markers.thigh_l_dist_lat) / 2;
        segment.dist_rad = norm(markers.leg_l_dist_med - markers.leg_l_dist_lat) / 2;
        segment.mass = subj.mass*0.0465;
    case 'foot_r'
        segment.len = norm(jc.ankle_global.right - 0.5*(markers.foot_r_dist_med + markers.foot_r_dist_lat));
        segment.prox_rad = norm(markers.leg_r_dist_med - markers.leg_r_dist_lat) / 2;
        segment.dist_rad = norm(markers.foot_r_dist_med - markers.foot_r_dist_lat) / 2;
        segment.mass = subj.mass*0.0145;
    case 'foot_l'
        segment.len = norm(jc.ankle_global.left - 0.5*(markers.foot_l_dist_med + markers.foot_l_dist_lat));
        segment.prox_rad = norm(markers.leg_l_dist_med - markers.leg_l_dist_lat) / 2;
        segment.dist_rad = norm(markers.foot_l_dist_med - markers.foot_l_dist_lat) / 2;
        segment.mass = subj.mass*0.0145;
    otherwise
        error(['Invalid segment name: ' segment_name]);
end

% Find segment centre of mass
segment.com = FindCom(segment.len, segment.prox_rad, segment.dist_rad);
segment.tensor = FindTensor(segment.len, segment.prox_rad, segment.dist_rad, segment.mass);

if strcmp(segment_name, 'pelvis')
    segment.com = segment.com + offset;
end
end

function com = FindCom(len, prox_rad, dist_rad)
% FindCom.m
% -------------------------------------------------------------------------
% Finds the centre of mass of a conical frustum with given length and
% radii.
% -------------------------------------------------------------------------
% Syntax and description:
% Centre of mass = FindCom(length, proximal radius, distal radius).
%
% The function takes the length, proximal radius, and distal radius of a
% conical frustum as input and returns a vector that gives the position of
% the frustum's centre of mass from its proximal end.
% -------------------------------------------------------------------------
% Written by Torstein E. Daehlin, August 2021.
% -------------------------------------------------------------------------

% Find centre of mass of conical frustum
if dist_rad < prox_rad
    x = dist_rad/prox_rad;
    sigma = 1 + x + x^2;
    z = ((1 + 2*x + 3*x^2)/(4*sigma)) * len;
else
    x = prox_rad/dist_rad;
    sigma = 1 + x + x^2;
    z = (1 - (1 + 2*x + 3*x^2)/(4*sigma)) * len;
end

% Assign output
com = [0, 0, -z];
end

function tensor = FindTensor(len, prox_rad, dist_rad, mass)
% FindTensor.m
% -------------------------------------------------------------------------
% Finds the inertia tensor about the centre of mass of a segment with the
% shape of a conical frustum.
% -------------------------------------------------------------------------
% Syntax and description:
% Inertia tensor = FindTensor(length, proximal radius, distal radius, mass)
%
% The function takes the length, proximal radius, distal radius, and mass
% of a conical frustum as input and returns the frustum's inertia tensor
% computed about its centre of mass.
% -------------------------------------------------------------------------
% Written by Torstein E. Daehlin, August, 2021.
% -------------------------------------------------------------------------

% Find inertia tensor at centre of mass of conical frustum with given mass
if dist_rad < prox_rad
    x = dist_rad/prox_rad;
else
    x = prox_rad/dist_rad;
end
sigma = 1 + x + x^2;
delta = 3 * mass/(pi*len*(prox_rad^2 + prox_rad*dist_rad + dist_rad^2));
a1 = 9/(20*pi);
a2 = (1 + x + x^2 + x^3 + x^4)/sigma^2;
b1 = 3/80;
b2 = (1 + 4*x + 10*x^2 + 4*x^3 + x^4)/sigma^2;
Ixx = a1+a2+mass^2/(delta*len) + b1*b2*mass*len^2;
Iyy = Ixx;
Izz = 2*a1*a2*mass^2/(delta*len);
tensor = [Ixx 0 0; 0 Iyy 0; 0 0 Izz];
end