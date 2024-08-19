function [leg_lcs, jc] = DefineLegLcs(markers, jc, nof, side)

% Preallocate
epx = zeros(nof,3);
epy = zeros(nof,3);
epz = zeros(nof,3);

% Define temporary vectors
if isequal(lower(side),'right')
    temp1 = 0.5*(markers.RMEP + markers.RLEP) - 0.5*(markers.RMMAL + markers.RLMAL);
    temp2 = markers.RLMAL - markers.RMMAL;
    origin = 0.5*(markers.RMEP + markers.RLEP);
elseif isequal(lower(side),'left')
    temp1 = 0.5*(markers.LMEP + markers.LLEP) - 0.5*(markers.LLMAL + markers.LMMAL);
    temp2 = markers.LMMAL - markers.LLMAL;
    origin = 0.5*(markers.LMEP + markers.LLEP);
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

if isequal(lower(side),'right')
    jc.ankle_global.(side) = 0.5*(markers.RMMAL + markers.RLMAL);
elseif isequal(lower(side),'left')
    jc.ankle_global.(side) = 0.5*(markers.LMMAL + markers.LLMAL);
end

% Transform knee centres to local system
R = [epx', epy', epz'];
jc.ankle.(side) = (jc.ankle_global.(side) - origin)*R;

% Assign ouput
leg_lcs.origin = origin;
leg_lcs.epx = epx;
leg_lcs.epy = epy;
leg_lcs.epz = epz;
end