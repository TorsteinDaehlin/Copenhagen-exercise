function grf = ForceProcess(data_struct, filter_parameters)

%{
Forces acting on the force platform (i.e. the action forces, not the ground
reaction forces) are expressed in the local coordiante system of the force
platform and must therefore be transformed to the global coordinate system.
                       _ _ _ _ _ _ _ _ _ _ _ _
                      /                       /
   Z                 /       x               /
   |                /       /               /
   |     Connector /    FP /_ _ _ y        /
LAB|_ _ _ _ Y     /       |               /
  /              /        |              /
 /              /         z             /
X              /_ _ _ _ _ _ _ _ _ _ _ _/
 
%}

% Set parameters
side = {'left','right'};
nof = data_struct.Force(1).NrOfSamples;

% Loop over force platforms
for i = 1:length(side)
    % Extract grf data from data structure
    grf.(side{i}).force = data_struct.Force(i).Force';
    grf.(side{i}).moment = data_struct.Force(i).Moment';
    grf.(side{i}).cop = data_struct.Force(i).COP'./1000;

    % Define force platform origin
    grf.(side{i}).corners = data_struct.Force(i).ForcePlateLocation./1000;
    grf.(side{i}).origin = mean(grf.(side{i}).corners);

    % Define force platform coordinate system
    temp_x = 0.5*(grf.(side{i}).corners(1,:) + grf.(side{i}).corners(4,:)) - ...
        0.5*(grf.(side{i}).corners(2,:) + grf.(side{i}).corners(3,:));
    epx = temp_x/norm(temp_x);
    temp_y = 0.5*(grf.(side{i}).corners(1,:) + grf.(side{i}).corners(2,:)) - ...
        0.5*(grf.(side{i}).corners(3,:) + grf.(side{i}).corners(4,:));
    temp_z = cross(epx,temp_y);
    epz = temp_z/norm(temp_z);
    epy = cross(epz,epx);
    R.(side{i}) = [epx' epy' epz'];

    % Calculate free moment
    grf.(side{i}).free_moment = zeros(nof,3);
    grf.(side{i}).free_moment(:,3) = grf.(side{i}).moment(:,3) - grf.(side{i}).cop(:,2).*grf.(side{i}).force(:,1) ...
        + grf.(side{i}).cop(:,1).*grf.(side{i}).force(:,2);

    % Transform quantities to the global coordinate system
    grf.(side{i}).force = ((R.(side{i})*grf.(side{i}).force')*-1)';
    grf.(side{i}).moment = ((R.(side{i})*grf.(side{i}).moment')*-1)';
    grf.(side{i}).free_moment = ((R.(side{i})*grf.(side{i}).free_moment')*-1)';
    grf.(side{i}).cop = ((R.(side{i})*grf.(side{i}).cop') + grf.(side{i}).origin')';

    % Filter force data
    filter_parameters.fs = data_struct.Force(i).Frequency;
    grf_filt = FilterData(grf.(side{i}), filter_parameters, 'force');

    % Downsample force
    grf.(side{i}).force = downsample(grf_filt.force,5);
    grf.(side{i}).moment = downsample(grf_filt.moment,5);
    grf.(side{i}).free_moment = downsample(grf_filt.free_moment,5);
    grf.(side{i}).cop = downsample(grf_filt.cop,5);
end
end