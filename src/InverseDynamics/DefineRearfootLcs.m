function rearfoot_lcs = DefineRearfootLcs(markers, nof, side)

% Preallocate
epx = zeros(nof,3);
epy = zeros(nof,3);
epz = zeros(nof,3);

% Define temporary vectors
if isequal(lower(side),'right')
    temp1 = markers.RPCAL - 0.5*(markers.RNT + markers.RCU);
    temp2 = markers.RCU - markers.RNT;
    origin = markers.RPCAL;
elseif isequal(lower(side),'left')
    temp1 = markers.LPCAL - 0.5*(markers.LNT + markers.LCU);
    temp2 = markers.LNT - markers.LCU;
    origin = markers.LPCAL;
else
    error(['Invalid side: ' side]);
end

% Define z-axis as the long axis of the segment
for i = 1:nof
    epz(i,:) = temp1/norm(temp1);
end

% Find the proportion of temp2 projected along the local z-axis
proj = dot(temp2',epz');

% Find the portion of temp2 perpendicular to the z-axis
for i = 1:nof
    temp3 = temp2(i,:) - proj(i)*epz(i,:);
    epx(i,:) = temp3/norm(temp3);
    epy(i,:) = cross(epz(i,:),epx(i,:));
end

% Assign ouput
rearfoot_lcs.origin = origin;
rearfoot_lcs.epx = epx;
rearfoot_lcs.epy = epy;
rearfoot_lcs.epz = epz;
end