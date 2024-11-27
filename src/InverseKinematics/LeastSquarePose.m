function [T] = LeastSquarePose(x_data,y_data)
% LeastSquarePose.m
% -------------------------------------------------------------------------
% Estimates the pose of a segment based on an arbitraty cluster of markers
% with known positions in the bone-embedded frame of the segment.
% -------------------------------------------------------------------------
% Syntax and description:
% [Rotation matrix, postition vector] = LeastSquarePose(bone-embedded data,
% global data)
%
% The function takes a 3 x m matrix with time-invariant positions of m
% markers expressed in the bone embedded frame and a 3 x m matrix of the
% position of the same markers in an instant in time measured during
% movement and expressed in the global coordinate system. The function
% returns the optimal estimates of the rotation matrix and postition vector
% of the segment for the current instant in time, computed using the least
% squares algorithm described by Söderkvist & Wedin (1993) and extended by
% Cappozzo et al. (1997).
% -------------------------------------------------------------------------
% Written by Torstein E. Daehlin, August, 2021.
% -------------------------------------------------------------------------

% References
%{
Cappozzo, A., Cappello, A., Croce, U. D., & Pensalfini, F. (1997).
    Surface-marker cluster design criteria for 3-D bone movement
    reconstruction. IEEE Transactions on Biomedical Engineering, 44(12),
    1165-1174.
Söderkvist, I., & Wedin, P. Å. (1993). Determining the movements
    of the skeleton using well-configured markers. Journal of biomechanics,
    26(12), 1473-1477.
%}

% Preallocate
X = zeros(size(x_data));
Y = zeros(size(y_data));

% Compute cluster model and marker centroid positions
x_mean = mean(x_data,2);
y_mean = mean(y_data,2);

% Compute cluster model position matrix
for i = 1:size(y_data,2)
    X(:,i) = x_data(:,i) - x_mean;
    Y(:,i) = y_data(:,i) - y_mean;
end

% Compute cluster cross-dispersion
Z = Y*X';

% Decompose cross-dispersion matrix using SVD
[U,~,V] = svd(Z);

% Compute least square estimates of rotation matrix and postition vector
R = U*diag([1 1 det(U*V')])*V';
p = y_mean - R*x_mean;

% Construct transformation matrix
T = [R p; 0 0 0 1];
end
