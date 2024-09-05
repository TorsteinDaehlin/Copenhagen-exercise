function grf = ForceProcess(force, filter_parameters)

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
nof = force(1).meta.nof;

% Loop over force platforms
for i = 1:length(force)
    % Define force platform origin
    force(i).origin = mean(force(i).fp_location);

    % Define force platform coordinate system
    temp_x = 0.5 * (force(i).fp_location(1,:) + force(i).fp_location(4,:)) - ...
        0.5 * (force(i).fp_location(2,:) + force(i).fp_location(3,:));
    epx = temp_x / norm(temp_x);
    temp_y = 0.5 * (force(i).fp_location(1,:) + force(i).fp_location(2,:)) - ...
        0.5 * (force(i).fp_location(3,:) + force(i).fp_location(4,:));
    temp_z = cross(epx,temp_y);
    epz = temp_z / norm(temp_z);
    epy = cross(epz,epx);
    R = [epx' epy' epz'];

    % Calculate free moment
    force(i).free_moment = zeros(nof * force(i).meta.SamplingFactor,3);
    force(i).free_moment(:,3) = force(i).moment(:,3) - force(i).cop(:,2) .* force(i).force(:,1) ...
        + force(i).cop(:,1) .* force(i).force(:,2);

    % Transform quantities to the global coordinate system
    force(i).force = ((R * force(i).force') * -1)';
    force(i).moment = ((R * force(i).moment') * -1)';
    force(i).free_moment = ((R * force(i).free_moment') * -1)';
    force(i).cop = ((R * force(i).cop') + force(i).origin')';

    % Filter force data
    filter_parameters.fs = force(i).meta.fs;
    grf_filt = FilterData(force(i), filter_parameters, 'force');

    % Downsample force
    ds_factor = force(i).meta.SamplingFactor;
    grf(i).force = downsample(grf_filt.force, ds_factor);
    grf(i).moment = downsample(grf_filt.moment, ds_factor);
    grf(i).free_moment = downsample(grf_filt.free_moment, ds_factor);
    grf(i).cop = downsample(grf_filt.cop, ds_factor);
    grf(i).origin = force(i).origin;
end
end