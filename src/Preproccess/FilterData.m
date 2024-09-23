function filtered_data = FilterData(raw_data, filter_parameters, data_type)

% Construct filter coefficients for butterworth filter
switch filter_parameters.type
    case 'Butterworth'
        [B, A] = butter(filter_parameters.order/2, ...
            filter_parameters.fc/(filter_parameters.fs/2), ...
            'low');

    case 'Chebyshev I'
        [B, A] = cheby1(filter_parameters.order/2, ...
            filter_parameters.fc/(filter_parameters.fs/2), ...
            'low');

    case 'Chebyshev II'
        [B, A] = cheby2(filter_parameters.order/2, ...
            filter_parameters.fc/(filter_parameters.fs/2), ...
            'low');

    case 'Ellipse'
        [B, A] = ellip(filter_parameters.order/2, ...
            filter_parameters.fc/(filter_parameters.fs/2), ...
            'low');

end    
    

if isequal(lower(data_type),'markers')
    
    % Get field names in raw data structure
    label_names = fieldnames(raw_data);
    
    % Filter raw data
    for i = 1:length(label_names)
        idx = isnan(raw_data.(label_names{i}));
        raw_data.(label_names{i})(idx) = 0;
        filtered_data.(label_names{i}) = filtfilt(B, A, raw_data.(label_names{i}));
    end

elseif isequal(lower(data_type),'force')

    % Get field names in raw data structure and define field names of
    % structures that need filtering
    var_names = {'force','moment','cop','free_moment'};

    % Filter raw data
    for i = 1:length(var_names)
        idx = isnan(raw_data.(var_names{i}));
        raw_data.(var_names{i})(idx) = 0;
        filtered_data.(var_names{i}) = filtfilt(B, A, raw_data.(var_names{i}));
    end
else
    % Print error message
    error(['Invalid data type name: ' data_type]);

end
end