function markers = DefineMarkers(data_struct)
% DefineMarkers
% -------------------------------------------------------------------------
% Extracts markers from .mat files exported from Qualisys Track Manager
% (QTM) and organizes the markers in a structure with fields correcsponding
% to each marker in the QTM-file. The fields are named using the marker
% labels given in QTM.
% -------------------------------------------------------------------------
% Syntax and description: markers = defineMarkers(data) returns
% a structure containing trajectories of markers labeled using QTM. The
% function takes the structure 'data' produced by QTM when motion capture
% data is exported in .mat format as input. The fields in the output
% structure are named according to the names of marker labels in the .qtm
% file. Marker trajectories are expressed in the global coordiante system.
% -------------------------------------------------------------------------
% Written by Torstein E. Daehlin, June 2019
% Revised by Torstein E. Daehlin, August 2021
% -------------------------------------------------------------------------

% Extract marker labels
label_names = data_struct.Trajectories.Labeled.Labels;

% Assign output structure
for i = 1:length(label_names)
    markers.(label_names{i}) = permute(data_struct.Trajectories.Labeled.Data(i,1:3,:),[3 2 1])/1000;
end
end
