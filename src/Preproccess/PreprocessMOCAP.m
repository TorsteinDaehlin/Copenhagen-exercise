function [static, dynamic, meta] = PreprocessMOCAP(subj, marker_reg, flt)

% Parse static file(s)
for i = 1:length(subj.static_name)
    % Load static file
    raw_static = load(fullfile(subj.data_path, subj.static_name{i}));

    % Parse file
    [static(i).markers, meta.static(i)] = ParseQTM(raw_static, marker_reg); 

    % Filter marker data
    flt.fs = meta.static(i).fs;
    static(i).markers = ...
        FilterData(static(i).markers, flt, 'markers');
end

% Parse dynamic files
for i = 1:length(subj.move_name)
    % Load dynamic file
    raw_dynamic = load(fullfile(subj.data_path, subj.move_name{i}));

    % Parse file
    [dynamic(i).markers, meta.dynamic(i), dynamic(i).force] = ParseQTM(raw_dynamic, marker_reg);

    % Filter markers
    flt.fs = meta.dynamic(i).fs * meta.dynamic(i).SamplingFactor;
    dynamic(i).markers = ...
        FilterData(dynamic(i).markers, flt, 'markers');

    % Process forces
    dynamic(i).force = ForceProcess(dynamic(i).force, meta.dynamic(i), flt);

end

if length(static) > 1
    % Prompt user to assign static file to distinct moving trials
    for i = 1:length(static)
        static_name = ['Static ' num2str(i)];
        idx = listdlg('ListString', subj.move_name, 'PromptString', ...
            ['Select trials to be matched with ' static_name]);
        static(i).match_to_move = subj.move_name(idx);
    end
else
    % Assume that all moving trials should be matched to the static trial
    static.match_to_move = subj.move_name;
end

end

function [markers, meta, varargout] = ParseQTM(qtm_struct, marker_reg)

struct_name = fieldnames(qtm_struct);

% Get metadata
meta.fs = qtm_struct.(struct_name{:}).FrameRate;
meta.nof = qtm_struct.(struct_name{:}).Frames;
meta.n_markers = qtm_struct.(struct_name{:}).Trajectories.Labeled.Count;

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
if nargout > 2
    
    % Get sampling factor
    meta.SamplingFactor = qtm_struct.(struct_name{:}).Force(1).SamplingFactor;

    for j = 1:length(qtm_struct.(struct_name{:}).Force)
        force(j).fp_location = qtm_struct.(struct_name{:}).Force(j).ForcePlateLocation./1000; % Convert to meters
        force(j).force = qtm_struct.(struct_name{:}).Force(j).Force';
        force(j).moment = qtm_struct.(struct_name{:}).Force(j).Moment';
        force(j).cop = qtm_struct.(struct_name{:}).Force(j).COP'./1000; % Convert to meters
    end
    
    % Assign force data to varargout
    varargout = {force};
end

end