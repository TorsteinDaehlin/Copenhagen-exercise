function forefoot_lcs = DefineForefootLcs(markers, nof, side)

% Preallocate
epx = zeros(nof,3);
epy = zeros(nof,3);
epz = zeros(nof,3);

% Define temporary vectors
if isequal(lower(side),'right')
    temp1 = markers.RMTH5 - 0.5*(markers.RNT + markers.RCU);
    temp2 = markers.RMTH1 - 0.5*(markers.RNT + markers.RCU);
    temp3 = 0.5*(markers.RNT + markers.RCU) - markers.RMTH2;
    origin = 0.5*(markers.RNT + markers.RCU);
elseif isequal(lower(side),'left')
    temp1 = markers.LMTH1 - 0.5*(markers.LNT + markers.LCU);
    temp2 = markers.LMTH5 - 0.5*(markers.LNT + markers.LCU);
    temp3 = 0.5*(markers.LNT + markers.LCU) - markers.LMTH2;
    origin = 0.5*(markers.LNT + markers.LCU);
else
    error(['Invalid side: ' side]);
end

% Define vector perpendicular to quasi-transverse plane formed by the
% midfoot centre and 1st and 5th metatarsal markers
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
forefoot_lcs.origin = origin;
forefoot_lcs.epx = epx;
forefoot_lcs.epy = epy;
forefoot_lcs.epz = epz;
end