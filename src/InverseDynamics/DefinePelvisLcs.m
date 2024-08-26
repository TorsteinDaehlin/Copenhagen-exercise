function [pelvis_lcs, jc] = DefinePelvisLcs(markers, nof)
% DefinePelvisLcs
% -------------------------------------------------------------------------
% Defines the local coordinate system of the pelvis segement in accordance
% with Cappozzo et al. 1995, Della Croce et al. 1999, and Wu et al. 2002.
% -------------------------------------------------------------------------
% Syntax and description: plevis_local_coordinate_system =
% DefinePelvisLcs(markers, number_of_frames) returns a structure containing
% time series for origin and local coordiante system axes of the pelvis.
% The function takes the structure 'markers' containing marker trajectories
% and the integer 'nof' giving the number of frames contained in the marker
% data set.
% -------------------------------------------------------------------------
% Written by Torstein E. Daehlin, August 2021
% -------------------------------------------------------------------------

% References:
%{
Cappozzo A, Catani F, Della Croce U, Leardini A. Position and orientation
    in space of bones during movement: anatomical frame definition and
    determination. Clin Biomech. 1995, 10(17), pp 1–8.

Della Croce U, Cappozzo A, Kerrigan DC. Pelvis and lower limb anatomical
    landmark calibration precision and its propaga- tion to bone geometry
    and joint angles. Med Biol Eng Comp. 1999, 37, pp 155–161.

Wu G, Siegler S, Allard P, Kirtley C, Leardini A, Rosenbaum D, et al. ISB
    recommendation on definitions of joint coordinate system of various
    joints for the reporting of human joint motion. Part 1: ankle, hip, and
    spine. J Biomech. 2002, 35, pp 543–548.
%}

% Preallocate
origin = zeros(nof,3);
epx = zeros(nof,3);
epy = zeros(nof,3);
epz = zeros(nof,3);

% Define temporary vectors
temp1 = markers.RASIS - markers.LASIS;
temp2 = 0.5*(markers.RASIS + markers.LASIS) - 0.5*(markers.RPSIS + markers.LPSIS);

% Define origin and local x-axis
for i = 1:nof
    origin(i,:) = 0.5*(markers.RASIS(i,:) + markers.LASIS(i,:));
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
jc.hip = FindHipCentre(markers, nof);

% Transform hip centres to globla system
R = [epx', epy', epz'];
jc.hip_global.right = jc.hip.right*R' + origin;
jc.hip_global.left = jc.hip.left*R' + origin;

% Assing output
pelvis_lcs.origin = origin;
pelvis_lcs.epx = epx;
pelvis_lcs.epy = epy;
pelvis_lcs.epz = epz;
end