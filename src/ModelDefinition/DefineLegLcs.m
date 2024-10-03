function [leg_lcs, jc] = DefineLegLcs(markers, jc, nof, side)

% Preallocate
epx = zeros(nof,3);
epy = zeros(nof,3);
epz = zeros(nof,3);

% Define temporary vectors
if isequal(lower(side),'right')
    temp1 = 0.5*(markers.thigh_r_dist_med + markers.thigh_r_dist_lat) - 0.5*(markers.leg_r_dist_med + markers.leg_r_dist_lat);
    temp2 = markers.leg_r_dist_lat - markers.leg_r_dist_med;
    origin = 0.5*(markers.thigh_r_dist_med + markers.thigh_r_dist_lat);
elseif isequal(lower(side),'left')
    temp1 = 0.5*(markers.thigh_l_dist_med + markers.thigh_l_dist_lat) - 0.5*(markers.leg_l_dist_med + markers.leg_l_dist_lat);
    temp2 = markers.leg_l_dist_med - markers.leg_l_dist_lat;
    origin = 0.5*(markers.thigh_l_dist_med + markers.thigh_l_dist_lat);
else
    error(['Invalid side: ' side]);
end

% Define local coordiante system
for i = 1:nof
    % The local x-axis is along the line connecting the medial and lateral
    % malleoli, pointing to the right
    epz(i,:) = temp1(i,:)/norm(temp1(i,:));

    % The local y-axis is perpendicular to the torsional plane, formed by
    % a vector connecting the knee and ankle joint centres and the local
    % x-axis, pointing anteriorly
    temp3 = cross(epz(i,:),temp2(i,:));
    epy(i,:) = temp3/norm(temp3);

    % The local z-axis is mutually perpendicular to the local x- and y-
    % axes
    epx(i,:) = cross(epy(i,:),epz(i,:));
end

% Find ankle joint centre
R = [epx', epy', epz'];

if isequal(lower(side),'right')
    jc.ankle_global.(side) = 0.5*(markers.leg_r_dist_med + markers.leg_r_dist_lat);
    jc.ankle_leg_r.(side) = (jc.ankle_global.(side) - origin)*R;
elseif isequal(lower(side),'left')
    jc.ankle_global.(side) = 0.5*(markers.leg_l_dist_med + markers.leg_l_dist_lat);
    jc.ankle_leg_l.(side) = (jc.ankle_global.(side) - origin)*R;
end

% Assign ouput
leg_lcs.origin = origin;
leg_lcs.epx = epx;
leg_lcs.epy = epy;
leg_lcs.epz = epz;
end