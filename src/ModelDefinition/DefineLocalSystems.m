function [lcs, jc] = DefineLocalSystems(markers)

% Define pelvis system
% Uses ISB recommendation + Harrington et al. 2007 for joint centres
[lcs.pelvis, jc] = DefinePelvisLcs(markers, 1);

% Define thigh systems
[lcs.thigh_r, jc] = DefineThighLcs(markers, jc, 1, 'right');
% [lcs.thigh_l, jc] = DefineThighLcs(markers, jc, 1, 'left');

% Define leg systems
[lcs.leg_r, jc] = DefineLegLcs(markers, jc, 1, 'right');
% [lcs.leg_l, jc] = DefineLegLcs(markers, jc, 1, 'left');

% Define foot coordinate system
[lcs.foot_r, jc] = DefineFootLcs(markers, jc, 1, 'right');
% lcs.foot_l = DefineFootLcs(markers, 1, 'left');

% Define vertical foot coordinate system
[lcs.v_foot_r] = DefineVerticalFootLcs(markers, 1, 'right');

% Define rearfoot coordinate systems
% lcs.rearfoot_r = DefineRearfootLcs(markers, 1, 'right');
% lcs.rearfoot_l = DefineRearfootLcs(markers, 1, 'left');

% Define forefoot coordinate system
% lcs.forefoot_r = DefineForefootLcs(markers, 1, 'right');
% lcs.forefoot_l = DefineForefootLcs(markers, 1, 'left');
end