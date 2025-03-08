function [lcs, jc] = DefineLocalSystems(markers)
% Local coordinate systems are defined as descibed by Wu et al 2002, while
% hip joint centres are calculated using the regression equation from
% Harrington et al. 2007. See main_Copenhagen() for details.

% Pelvis
[lcs.pelvis, jc] = DefinePelvisLcs(markers, 1);

% Thigh
[lcs.thigh_r, jc] = DefineThighLcs(markers, jc, 1, 'right');

% Leg
[lcs.leg_r, jc] = DefineLegLcs(markers, jc, 1, 'right');

% Foot
[lcs.foot_r, jc] = DefineFootLcs(markers, jc, 1, 'right');

% We use a "vertical foot" segment to compute foot kinematics
[lcs.v_foot_r] = DefineVerticalFootLcs(markers, 1, 'right');
end