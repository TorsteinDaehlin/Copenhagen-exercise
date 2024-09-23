function foot_lcs = DefineFootLcs(markers, nof, side)

% Preallocate
epx = zeros(nof,3);
epy = zeros(nof,3);
epz = zeros(nof,3);

% Define temporary vectors
if isequal(lower(side),'right')
    temp1 = markers.foot_r_dist_lat - 0.5 * (markers.leg_r_dist_med + markers.leg_r_dist_lat);
    temp2 = markers.foot_r_dist_med - 0.5 * (markers.leg_r_dist_med + markers.leg_r_dist_lat);
    temp3 = 0.5 * (markers.leg_r_dist_med + markers.leg_r_dist_lat) - ...
        0.5 * (markers.foot_r_dist_med + markers.foot_r_dist_lat);
    origin = 0.5 * (markers.leg_r_dist_med + markers.leg_r_dist_lat);
elseif isequal(lower(side),'left')
    temp1 = markers.foot_l_dist_lat - 0.5 * (markers.leg_l_dist_med + markers.leg_l_dist_lat);
    temp2 = markers.foot_l_dist_med - 0.5 * (markers.leg_l_dist_med + markers.leg_l_dist_lat);
    temp3 = 0.5 * (markers.leg_l_dist_med + markers.leg_l_dist_lat) - ...
        0.5 * (markers.foot_l_dist_med + markers.foot_l_dist_lat);
    origin = 0.5 * (markers.leg_l_dist_med + markers.leg_l_dist_lat);
else
    error(['Invalid side: ' side]);
end

% Define vector perpendicular to quasi-transverse plane formed by calcaneus
% markers and 1st and 5th metatarsal markers
for i = 1:nof
    temp_n = cross(temp1(i,:),temp2(i,:));
    epy(i,:) = temp_n/norm(temp_n);
end

% Find the proportion of temp3 projected along n_hat
proj = dot(temp3',epy');

for i = 1:nof
    temp4 = temp3(i,:) - proj(i)*epy(i,:);
    epz(i,:) = temp4/norm(temp4);
    epx(i,:) = cross(epy(i,:),epz(i,:));
end

% Assign ouput
foot_lcs.origin = origin;
foot_lcs.epx = epx;
foot_lcs.epy = epy;
foot_lcs.epz = epz;
end
