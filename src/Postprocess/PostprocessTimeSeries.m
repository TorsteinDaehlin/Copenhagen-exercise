function ts = PostprocessTimeSeries(time, kinematics, grf, njm, roi)

% Loop over ROIs and extract variables of interest
n_roi = size(roi, 1);
for i = 1:n_roi
    % Store time series
    % -----------------
    % Time
    ts.(['repetition_' num2str(i)]).time = time(roi(i, 1):roi(i, 2));
    
    % Kinematics
    vars = fieldnames(kinematics);
    vars(contains(vars, {'dynamic_lcs', 'jc'})) = [];
    for j = 1:length(vars)
        lbls = fieldnames(kinematics.(vars{j}));
        for k = 1:length(lbls)
            if isstruct(kinematics.(vars{j}).(lbls{k}))
                ts.(['repetition_' num2str(i)]).([vars{j} '_' lbls{k}]) = ...
                    kinematics.(vars{j}).(lbls{k}).right(roi(i, 1):roi(i, 2), :);
            else
                ts.(['repetition_' num2str(i)]).([vars{j} '_' lbls{k}]) = ...
                    kinematics.(vars{j}).(lbls{k})(roi(i, 1):roi(i, 2), :);
            end
        end
    end

    % GRF
    vars = {'force', 'cop', 'free_moment'};
    for j = 1:length(vars)
        ts.(['repetition_' num2str(i)]).(vars{j}) = grf.(vars{j})(roi(i, 1):roi(i, 2), :);
    end

    % NJM
    vars = fieldnames(njm);
    for j = 1:length(vars)
        ts.(['repetition_' num2str(i)]).(['njm_' vars{j}]) = njm.(vars{j})(roi(i, 1):roi(i, 2), :);
    end
end
end