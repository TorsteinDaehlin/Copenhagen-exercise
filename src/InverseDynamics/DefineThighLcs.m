function [thigh_lcs, jc] = DefineThighLcs(markers, jc, nof, side)

% Preallocate
origin = zeros(nof,3);
epx = zeros(nof,3);
epy = zeros(nof,3);
epz = zeros(nof,3);

% Define temporary vectors
if isequal(lower(side),'right')
    temp1 = jc.hip_global.(lower(side)) - 0.5*(markers.RMEP + markers.RLEP);
    temp2 = markers.RLEP - markers.RMEP;
elseif isequal(lower(side),'left')
    temp1 = jc.hip_global.(lower(side)) - 0.5*(markers.LMEP + markers.LLEP);
    temp2 = markers.LMEP - markers.LLEP;
else
    error(['Invalid side: ' side]);
end

% Define local coordinate system
for i = 1:nof
    % Origin coincides with hip centre
    origin(i,:) = jc.hip_global.(lower(side))(i,:);

    % The local z-axis is pointing from the distal to the proximal joint
    % centre
    epz(i,:) = temp1(i,:)/norm(temp1(i,:));
end

% Find the portion of temp2 projected onto the local z-axis
proj = dot(temp2',epz');

for i = 1:nof
    temp3 = temp2(i,:) - proj(i)*epz(i,:);
    epx(i,:) = temp3/norm(temp3);
    epy(i,:) = cross(epz(i,:),epx(i,:));
end

% Find knee joint centre
if isequal(lower(side),'right')
    jc.knee_global.(side) = 0.5*(markers.RMEP + markers.RLEP);
elseif isequal(lower(side),'left')
    jc.knee_global.(side) = 0.5*(markers.LMEP + markers.LLEP);
end

% Transform knee centres to local system
R = [epx', epy', epz'];
jc.knee.(side) = (jc.knee_global.(side) - origin)*R;

% Assign ouput
thigh_lcs.origin = origin;
thigh_lcs.epx = epx;
thigh_lcs.epy = epy;
thigh_lcs.epz = epz;
end
