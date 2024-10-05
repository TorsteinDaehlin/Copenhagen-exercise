function [varargout] = EulerAngles(R, varargin)
% EulerAngles.m
% -------------------------------------------------------------------------
% Calculates angular rotations about the principals axes of a given
% rotation matrix.
% -------------------------------------------------------------------------
% Syntax and description:
% [alpha, beta, gamma] = EulerAngles(R) takes the rotation matrix R as
% input and returns the angles alpha, beta, and gamma describing the
% rotations about the 1st, 2nd, and 3rd axis in the rotation sequence,
% respectively. The rotation sequence 'xyz' is used when only R is provided
% as input, and the resulting angles are expressed in radians.
%
% [alpha, beta, gamma] = EulerAngles(R,sequence) uses the input argument
% 'sequence' to speficy the desired rotation sequence. The input argument
% sequence must be one of the following:
%    'xyz'
%    'xzy'
%    'yxz'
%    'yzx'
%    'zxy'
%    'zyx'
%    'xyx'
%    'xzx'
%    'yxy'
%    'yzy'
%    'zxz'
%    'zyz'
% If an empty string is provided as input for sequence, the function
% returns the 'xyz' rotation sequence angles by default.
%
% [alpha, beta, gamma] = EulerAngles(R,'','deg') returns the output angles
% alpha, beta, and gamma in degrees rather than radians which is the
% default output format. If an empty input is provided, the angles will be
% expressed in radians.
% -------------------------------------------------------------------------
% Written by Torstein E. Daehlin, January 2019.
% -------------------------------------------------------------------------

% Error check input arguments
% Determine number of inputs
n = nargin;

% Error check number of inputs and required input
if n > 3
    error('Error using EulerAngles! Too many input arguments');
elseif n < 1
    error('Error using EulerAngles! Not enough input arguments');
else
    if ~isequal(size(R),[3 3]) % Checks the first input argument
        error('Error using EulerAngles! Rotation matrix must be of dimension 3x3');
    elseif ~isnumeric(R)
        error('Error using EulerAngles! Rotation matrix must contain only numeric values');
    end
end

% Error check optional input arguments
if n >= 2 % Checks the second input argument
    if ~ischar(varargin{1})
        error('Error using EulerAngles! Input ''seq'' must be of type char');
    elseif ~isequal(length(varargin{1}),3)
        error('Error using EulerAngles! Invalid input argument: %s',varargin{1});
    end
end

if n == 3 % Checks the third input argument
    if ~ischar(varargin{2})
        error('Error using EulerAngles! Input ''deg'' must be of type char');
    elseif ~isequal(length(varargin{2}),length('deg'))
        error('Error using EulerAngles! Invalid input argument: %s',varargin{2});
    elseif ~isequal(varargin{2},'deg')
        error('Error using EulerAngles! Invalid input argument: %s',varargin{2});
    end
end

% Calculate joint angles alpha, beta, and gamma
% Select appropriate sequence for specified case
if n == 1
    seq = 'xyz';
elseif isempty(varargin{1})
    seq = 'xyz';
else
    seq = varargin{1};
end

switch seq
    case 'xyz'
        x_angle = atan2(-R(2,3),R(3,3));
        y_angle = atan2(R(1,3),sqrt(R(2,3)^2 + R(3,3)^2));
        z_angle = atan2(-R(1,2),R(1,1));
        varargout{1} = x_angle;
        varargout{2} = y_angle;
        varargout{3} = z_angle;

    case 'xzy'
        x_angle = atan2(R(3,2),R(2,2));
        z_angle = atan2(-R(1,2),sqrt(R(3,2)^2 + R(2,2)^2));
        y_angle = atan2(R(1,3),R(1,1));
        varargout{1} = x_angle;
        varargout{2} = y_angle;
        varargout{3} = z_angle;

    case 'yxz'
        y_angle = atan2(R(1,3),R(3,3));
        x_angle = atan2(-R(2,3),sqrt(R(1,3)^2 + R(3,3)^2));
        z_angle = atan2(R(2,1),R(2,2));
        varargout{1} = x_angle;
        varargout{2} = y_angle;
        varargout{3} = z_angle;

    case 'yzx'
        y_angle = atan2(-R(3,1),R(1,1));
        z_angle = atan2(R(2,1),sqrt(R(1,1)^2 + R(3,1)^2));
        x_angle = atan2(-R(2,3),R(2,2));
        varargout{1} = x_angle;
        varargout{2} = y_angle;
        varargout{3} = z_angle;

    case 'zxy'
        z_angle = atan2(-R(1,2),R(2,2));
        x_angle = atan2(R(3,2),sqrt(R(1,2)^2 + R(2,2)^2));
        y_angle = atan2(-R(3,1),R(3,3));
        varargout{1} = x_angle;
        varargout{2} = y_angle;
        varargout{3} = z_angle;

    case 'zyx'
        z_angle = atan2(R(2,1),R(1,1));
        y_angle = atan2(-R(3,1),sqrt(R(1,1)^2 + R(2,1)^2));
        x_angle = atan2(R(3,2),R(3,3));
        varargout{1} = x_angle;
        varargout{2} = y_angle;
        varargout{3} = z_angle;

    case 'xyx'
        alpha = atan2(R(2,1),-R(3,1));
        beta = atan2(sqrt(R(1,2)^2 + R(1,3)^2),R(1,1));
        gamma = atan2(R(1,2),R(1,3));
        varargout{1} = alpha;
        varargout{2} = beta;
        varargout{3} = gamma;

    case 'xzx'
        alpha = atan2(R(3,1),R(2,1));
        beta = atan2(sqrt(R(2,1)^2 + R(3,1)^2),R(1,1));
        gamma = atan2(R(1,3),-R(1,2));
        varargout{1} = alpha;
        varargout{2} = beta;
        varargout{3} = gamma;

    case 'yxy'
        alpha = atan2(R(1,2),R(3,2));
        beta = atan2(sqrt(R(2,1)^2 + R(2,3)^2),R(2,2));
        gamma = atan2(R(2,1),-R(2,3));
        varargout{1} = alpha;
        varargout{2} = beta;
        varargout{3} = gamma;

    case 'yzy'
        alpha = atan2(R(3,2),-R(1,2));
        beta = atan2(sqrt(R(2,1)^2 + R(2,3)^2),R(2,2));
        gamma = atan2(R(2,3),R(2,1));
        varargout{1} = alpha;
        varargout{2} = beta;
        varargout{3} = gamma;

    case 'zxz'
        alpha = atan2(R(1,3),-R(2,3));
        beta = atan2(sqrt(R(3,1)^2 + R(3,2)^2),R(3,3));
        gamma = atan2(R(3,1),R(3,2));
        varargout{1} = alpha;
        varargout{2} = beta;
        varargout{3} = gamma;

    case 'zyz'
        alpha = atan2(R(2,3),R(1,3));
        beta = atan2(sqrt(R(1,3)^2 + R(2,3)^2),R(3,3));
        gamma = atan2(R(3,2),-R(3,1));
        varargout{1} = alpha;
        varargout{2} = beta;
        varargout{3} = gamma;

    otherwise
        error('Error using EulerAngles! Invalid input argument: %s',varargin{1});
end

% Convert to joint angles from radians to degrees (if specified by input arguments)
if n == 3
    varargout{1} = rad2deg(varargout{1});
    varargout{2} = rad2deg(varargout{2});
    varargout{3} = rad2deg(varargout{3});
end
end
