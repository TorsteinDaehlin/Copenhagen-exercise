function foot_lcs = DefineFootLcs_calc(markers, nof, side)

% Preallocate
epx = zeros(nof,3);
epy = zeros(nof,3);
epz = zeros(nof,3);

% Define temporary vectors
if isequal(lower(side),'right')
    temp1 = markers.RMTH5 - markers.RPCAL;
    temp2 = markers.RMTH1 - markers.RPCAL;
    temp3 = markers.RPCAL - markers.RMTH2;
    origin = markers.RPCAL;
elseif isequal(lower(side),'left')
    temp1 = markers.LMTH1 - markers.LPCAL;
    temp2 = markers.LMTH5 - markers.LPCAL;
    temp3 = markers.LPCAL - markers.LMTH2;
    origin = markers.LPCAL;
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
