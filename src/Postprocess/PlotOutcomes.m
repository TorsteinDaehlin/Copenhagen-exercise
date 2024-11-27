function PlotOutcomes()

% Load data
[f_name, f_path] = uigetfile('..\*.mat', 'Select results file');
load(fullfile(f_path, f_name), 'tbls');

% Here we select the joints to be plotted
joints = {'hip', 'knee', 'ankle'};

conditions = fieldnames(tbls);
for i = 1:length(joints)
    for j = 1:length(conditions)
        subplot(2, length(joints), i);
        
    end
end

end