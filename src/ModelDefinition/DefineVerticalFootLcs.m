function [v_foot_lcs] = DefineVerticalFootLcs(markers, nof, side)

% Preallocate
epx = zeros(nof,3);
epy = zeros(nof,3);
epz = zeros(nof,3);

% Define temporary vectors
if isequal(lower(side),'right')
    temp1 = markers.foot_r_dist_lat - 0.5 * ([markers.leg_r_dist_med(1:2) 0] + [markers.leg_r_dist_lat(1:2) 0]);
    temp2 = markers.foot_r_dist_med - 0.5 * ([markers.leg_r_dist_med(1:2) 0] + [markers.leg_r_dist_lat(1:2) 0]);
    temp3 = 0.5 * (markers.foot_r_dist_med + markers.foot_r_dist_lat) - ...
        0.5 * ([markers.leg_r_dist_med(1:2) 0] + [markers.leg_r_dist_lat(1:2) 0]);
    origin = 0.5 * ([markers.leg_r_dist_med(1:2) 0] + [markers.leg_r_dist_lat(1:2) 0]);
elseif isequal(lower(side),'left')
    temp1 = markers.foot_l_dist_lat - 0.5 * ([markers.leg_l_dist_med(1:2) 0] + [markers.leg_l_dist_lat(1:2) 0]);
    temp2 = markers.foot_l_dist_med - 0.5 * ([markers.leg_l_dist_med(1:2) 0] + [markers.leg_l_dist_lat(1:2) 0]);
    temp3 = 0.5 * (markers.foot_l_dist_med + markers.foot_l_dist_lat) - ...
        0.5 * ([markers.leg_l_dist_med(1:2) 0] + [markers.leg_l_dist_lat(1:2) 0.0065]);
    origin = 0.5 * ([markers.leg_l_dist_med(1:2) 0] + [markers.leg_l_dist_lat(1:2) 0]);
else
    error(['Invalid side: ' side]);
end

% Define vector perpendicular to quasi-transverse plane formed by calcaneus
% markers and 1st and 5th metatarsal markers
for i = 1:nof
    temp_n = cross(temp1(i,:),temp2(i,:));
    epz(i,:) = temp_n/norm(temp_n);
end

% Find the proportion of temp3 projected along temp_n
proj = dot(temp3',epz');

for i = 1:nof
    temp4 = temp3(i,:) - proj(i)*epz(i,:);
    epy(i,:) = temp4/norm(temp4);
    epx(i,:) = cross(epy(i,:),epz(i,:));
end

% Assign ouput
v_foot_lcs.origin = origin;
v_foot_lcs.epx = epx;
v_foot_lcs.epy = epy;
v_foot_lcs.epz = epz;
end
