function filtered_data = FilterData(raw_data, filter_parameters, type)

if isequal(lower(type),'markers')
    % Construct filter coefficients for butterworth filter
    [B, A] = butter(filter_parameters.order/2, filter_parameters.fc/(filter_parameters.fs/2), filter_parameters.type);

    % Get field names in raw data structure
    label_names = fieldnames(raw_data);

    % Filter raw data
    for i = 1:length(label_names)
        filtered_data.(label_names{i}) = filtfilt(B, A, raw_data.(label_names{i}));
    end

elseif isequal(lower(type),'force')
    % Construct filter coefficients for butterworth filter
    [B, A] = butter(filter_parameters.order/2, filter_parameters.fc/(filter_parameters.fs/2), filter_parameters.type);

    % Get field names in raw data structure and define field names of
    % structures that need filtering
    var_names = {'force','moment','cop','free_moment'};

    % Filter raw data
    for j = 1:length(var_names)
        filtered_data.(var_names{j}) = filtfilt(B, A, raw_data.(var_names{j}));
    end
else
    % Print error message
    error(['Invalid data type name: ' type]);

end
end