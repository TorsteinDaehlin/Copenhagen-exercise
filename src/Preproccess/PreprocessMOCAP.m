function [] = PreprocessMOCAP(subj)

% Load static file(s)
for i = 1:length(subj.static_name)
    static(i) = load(fullfile(subj.data_path, subj.static_name{i}));
end

end

