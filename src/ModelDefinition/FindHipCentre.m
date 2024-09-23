function hip_centre = FindHipCentre(markers, nof)

% Preallocate
hip_centre.right = zeros(nof,3);
hip_centre.left = zeros(nof,3);

% Constants from Harrington et al. 2007 (y = ax + b)
a = [0.33, -0.24, -0.30];
b = [7.3 -9.9 -10.9]./1000; % divided by 1000 to convert to meters

% Variables
for i = 1:nof
    pelvis_width = norm(markers.pelvis_RASIS(i,:) - markers.pelvis_LASIS(i,:));
    pelvis_depth = norm(0.5*(markers.pelvis_RASIS(i,:) + markers.pelvis_LASIS(i,:)) ...
        - 0.5*(markers.pelvis_RPSIS(i,:) + markers.pelvis_LPSIS(i,:)));

    % Solve regression equations
    x_hat = a(1)*pelvis_width + b(1);
    y_hat = a(2)*pelvis_depth + b(2);
    z_hat = a(3)*pelvis_width + b(3);

    % Define hip centres
    hip_centre.right(i,:) = [x_hat, y_hat, z_hat];
    hip_centre.left(i,:) = [-x_hat, y_hat, z_hat];
end
end