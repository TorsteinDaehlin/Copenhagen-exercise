function [static, dynamic] = PreprocessMOCAP(subj, marker_reg)

% Parse static file(s)
for i = 1:length(subj.static_name)
    % Load static file
    raw_static(i) = load(fullfile(subj.data_path, subj.static_name{i}));

    % Parse file
    static(i).markers = ParseQTM(raw_static(i), marker_reg); 

    % Filter marker data
    flt.fs = static(i).markers.meta.fs;
    static(i).markers = ...
        FilterData(static(i).markers, flt, 'markers');
end

% Parse dynamic files
for i = length(subj.move_name)
    % Load dynamic file
    raw_dynamic(i) = load(fullfile(subj.data_path, subj.move_name{i}));

    % Parse file
    [dynamic(i).markers, dynamic(i).force] = ParseQTM(raw_dynamic(i), marker_reg);

    % Filter markers
    flt.fs = dynamic(i).markers.meta.fs;
    dynamic(i).markers = ...
        FilterData(dynamic(i).markers, flt, 'markers');

    % Process forces
    dynamic(i).force = ForceProcess(dynamic(i).force, flt);

end

if length(static) > 1
    % Prompt user to assign static file to distinct moving trials
else
    % Assume that all moving trials should be matched to the static trial
end

end

function [markers, varargout] = ParseQTM(qtm_struct, marker_reg)

struct_name = fieldnames(qtm_struct);

% Get metadata
markers.meta.fs = qtm_struct.(struct_name{:}).FrameRate;
markers.meta.nof = qtm_struct.(struct_name{:}).Frames;
markers.meta.n_markers = qtm_struct.(struct_name{:}).Trajectories.Labeled.Count;

% Get marker labels
mrk_labels = (qtm_struct.(struct_name{:}).Trajectories.Labeled.Labels)';

% Assign marker to structure
for i = 1:length(mrk_labels)
    idx = strcmp(marker_reg.label, mrk_labels(i));
    generic_label = [marker_reg.segment{idx} '_' marker_reg.type{idx}];
    markers.(generic_label) = ...
        permute(qtm_struct.(struct_name{:}).Trajectories.Labeled.Data(i,1:3,:),[3 2 1])/1000;
end

% Parse forces
if nargout > 1
    force.meta.SamplinFactor = qtm_struct.(struct_name{:}).Force(1).SamplingFactor;
    force.meta.nof = qtm_struct.(struct_name{:}).Force(1).NrOfFrames;
    force.meta.fs = qtm_struct.(struct_name{:}).Force(1).Frequency;

    for j = 1:length(qtm_struct.(struct_name{:}).Force)
        force(j).fp_location = qtm_struct.(struct_name{:}).Force(j).ForcePlateLocation./1000; % Convert to meters
        force(j).force = qtm_struct.(struct_name{:}).Force(j).Force';
        force(j).moment = qtm_struct.(struct_name{:}).Force(j).Moment';
        force(j).cop = qtm_struct.(struct_name{:}).Force(j).COP'./1000; % Convert to meters
    end
    
    % Assign force data to varargout
    varargout = force;
end

end