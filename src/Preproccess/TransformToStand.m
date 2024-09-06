function dynamic = TransformToStand(dynamic)

% Define box weight;
z = -0.038; % m
m_box = 2.6; % kg
g = -9.81; % m/s^2
w_box = [0; 0; m_box * g];

for i = 1:length(dynamic)
    % Preallocate loop variables for speed
    nof = size(dynamic(i).force(2).force, 1);
    GRF = zeros(nof, 3);
    moment = zeros(nof, 3);
    cop = zeros(nof, 3);

    for j = 1:nof

        % Define stand coordinate system
        origin = 0.5 * (dynamic(i).markers.external_cluster1(j,:) ...
            + dynamic(i).markers.external_cluster2(j,:));
        epx = dynamic(i).markers.external_cluster1(j,:) ...
            - dynamic(i).markers.external_cluster2(j,:);
        epx = epx / norm(epx);

        epz = dynamic(i).markers.external_cluster3(j,:) ...
            - dynamic(i).markers.external_cluster2(j,:);
        epz = epz / norm (epz);

        epy = cross(epz, epx);
        epx = cross(epy, epz);

        R = [epx' epy' epz'];

        % Calculate force and moment acting at the top of the stand
        F = (dynamic(i).force(2).force(j,:) + w_box') .* -1;

        r_com = [origin(1:2) origin(3) / 2];
        r_cop = dynamic(i).force(2).cop(j,:) - origin;
        tau = (cross(r_com, w_box) + ...
            cross(r_cop, dynamic(i).force(2).force(j,:))) .* -1;

        % Transform force and moment to stand coordinate system      
        F = (R' * F')';
        tau = (R' * tau')';

        % Resolve into equivalent wrench
        p_x = (tau(2) - z * F(1)) / F(3);
        p_y = (tau(1) + z * F(2)) / F(3);
        p = [p_x, p_y, z];

        % Transform forces, moment, and point of application back to global
        % coordinate system
        GRF(j,:) = (R * F')' .* -1;
        moment(j,:) = (R * tau')' .* -1;
        cop(j,:) = (R * p')' + origin;
        
    end

    % Replace force in dynamic
    dynamic(i).force(2).force = GRF;
    dynamic(i).force(2).moment = moment;
    dynamic(i).force(2).cop = cop;

    clear GRF Moment cop;
end
end