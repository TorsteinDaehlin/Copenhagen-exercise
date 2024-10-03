% main_Copenhagen.m
% -------------------------------------------------------------------------
% Performes inverse dynamics batch processing of vertical and approach
% jumps collected using a Qualisys motion capture system, 2 AMTI force
% platforms, and a CAST marker set consisting of 28 anatomical markers and
% 27 tracking markers for the incline vs. block heel-raise research project
% conducted in the Sports Biomechanics Laboratory at the University of
% Alberta, AB, Canada 2020-2021.
% -------------------------------------------------------------------------
% Syntax and description: InclineVsBlock2021() - The program prompts the
% user to select a participant directory and subsequently inspects the
% directory, where it expects to find one or more of the following visit
% directories containing .mat-files exported from Qualisys Track Manager:
% 'pretest', 'midtest', and 'posttest' (not case sensitive). If the folder
% does not contain any of the above folders, the program termiantes with an
% error. If only some are present, a warning is printed to the command
% window and to the execution log.
%
% The program then loops the visit directories and inspects them for
% static-files, left and right approach jump files, and vertical jump
% files. Errors are produced if no static file or no motion file of each of
% the jump types are not present. After completion of the error check, the
% program filters the input data and performs inverse dynamics. The program
% prints results files (.mat-format) for each of the trials. These can be
% post processed using the program named "INSERT PROGRAM NAME HERE".
%
% For further details regarding this program and its subroutines, including
% references to research literature, see comments in source code files.
%
% -------------------------------------------------------------------------
% Written by Torstein E. Daehlin, PhD August 2024.
% -------------------------------------------------------------------------

% References:
% -----------
%{


%}
%
%==========================================================================

function main_Copenhagen()

% Initialize workspace
% ====================

% Add paths to dependencies
addpath('.\Initialize\');

% Prompt user to select participant directory
[subj_dir, src_path] = GetSubjDir();

% Look for marker registration file in participant directory. If not there,
% prompt user to open marker registration file
if isfile(fullfile(subj_dir.folder,'model_setup.csv'))
    marker_reg = readtable(fullfile(subj_dir.folder, 'model_setup.csv'));
else
    [reg_file, reg_path] = uigetfile('..\..\*.csv');
    marker_reg = readtable(fullfile(reg_path, reg_file));
end

% Create output directory
dst_path = regexp(src_path, filesep, 'split');
dst_path = fullfile(dst_path{1:end-1},'Output');

if ~isfolder(dst_path)
    mkdir(dst_path);
end

% Requrest processing input (e.g. filter parameters)
[static_id, flt] = GetProcessingParameters();

if isempty(static_id) || isempty(flt)
    fprintf('Please fill provide all processing parameters\n');
    [static_id, flt] = GetProcessingParameters();
end

for s = 1:length(subj_dir)
    % Inspect participant directory
    % =============================
    % Store dirctory as participant name
    subj.id = subj_dir(s).name;

    % Prompt user to provide height and mass
    prompt = {'Enter participant height (m):', ...
        'Enter participant mass (kg):'};
    answer = inputdlg(prompt, [subj.id ': Height and mass']);
    subj.height = str2double(answer{1});
    subj.mass = str2double(answer{2});

    % Get directory contents
    motion_files = dir(fullfile(subj_dir(s).folder, subj_dir(s).name));
    motion_files = motion_files(~(strcmp({motion_files.name}, {'.'}) | ...
        strcmp({motion_files.name}, {'..'})));

    % Error check input
    static_idx = contains({motion_files.name}, static_id);
    if ~any(static_idx)
        error('No static trial whose file name contains "%s".', static_id);
    end

    % Get file name(s) and directories
    subj.data_path = motion_files(1).folder;
    subj.static_name = {motion_files(static_idx).name};
    subj.move_name = {motion_files(~static_idx).name};

    % Preprocess data
    % ===============
    % Add path to dependencies
    addpath('.\Preproccess\');
    addpath('.\ModelDefinition\');
    addpath('.\InverseDynamics\');
    addpath('.\InverseKinematics\');

    % Preprocess input data
    [static, dynamic, meta] = PreprocessMOCAP(subj, marker_reg, flt);

    % Transform force to top of stand (the easiest way to achieve this may be to simply add this as a rigid body
    % to the model. Since there is not angular velocity of the stand, point of force application should be possible to find)
    dynamic = TransformToStand(dynamic);

    % Run inverse dynamics procedure
    % ==============================
    for i = 1:length(static)
        % Generate participant model
        % --------------------------
        [static_lcs, static_jc, segments] = ...
            ProcessStatic(static(i), meta.static(i), subj);

        % Loop over dynamic trials matched to current static
        for j = 1:length(static.match_to_move)

            % Perform inverse kinematics
            [kinematics, time] = ...
                InverseKinematics(dynamic(j).markers, static.markers, ...
                static_lcs, static_jc, segments, meta.dynamic(j));            

            % Determine which external forces are applied to which segments

            % Calculate NJMs using inverse dynamics
        end
    end
    % 
    
    for j = 1:length(static_files)



        % Define subject model

        % Loop over dynamic files
        % -----------------------
        % Initialize variables
        if isequal(static_name,'Static')
            jump_height = struct('AJL', 0, 'AJR', 0, 'VJ', 0);
            output = struct('AJL', [], 'AJR', [], 'VJ', []);
        end

        for i = 1:length(trial_files)
            % Find trial name
            trial_name = strtok(trial_files(i).name,'.');
            fprintf('\nNow processing: %s\n--------------------\n', trial_name);

            % Load dynamic file
            load(fullfile(trial_files(i).folder,trial_files(i).name),trial_name);

            % Process dynamic file
            [kinematics, kinetics, time] = ProcessDynamic(static_markers, static_lcs, static_jc, segments, eval(trial_name), filter_parameters, ...
                trial_name, visit_name);

            % Label events and extract values for export
            events = FindEvents(kinematics, kinetics, trial_name(1:end-1));

            % Extract values for export %MAKE SURE THAT ALL VARIABLES OF
            % INTEREST ARE BEING EXPORTED
            switch trial_name(1:end-1)
                case 'AJL'
                    [output.AJL, jump_height.AJL] = ExtractOutput(static_lcs,segments,kinematics,kinetics,events,jump_height.AJL,output.AJL,time);
                case 'AJR'
                    [output.AJR, jump_height.AJR] = ExtractOutput(static_lcs,segments,kinematics,kinetics,events,jump_height.AJR,output.AJR,time);
                case 'VJ'
                    [output.VJ, jump_height.VJ] = ExtractOutput(static_lcs,segments,kinematics,kinetics,events,jump_height.VJ,output.VJ,time);
                otherwise
                    error(['Invalid trial type: ' trial_name(1:end-1)]);
            end
        end
    end
    % Save output structure
    save(fullfile(participant_path,'Results',[participant.name '_' visit_name '.mat']),'output');
end
end







function dydx = FiniteDiff(y,h,varagin)
% finiteDiff.m
% -------------------------------------------------------------------------
% Differentiates the equally spaced dependent variables y = f(x) using a
% finite difference scheme. The function can return the first or second
% derivative of y.
% -------------------------------------------------------------------------
% Syntax and description:
% dydx = finiteDiff(y,h) returns the first derivative of the input array y
% = f(x) with equally spaced steps specified by the step size h.
%
% dydx = finiteDiff(y,h,order) returns the derivative of the specified
% order. The input argument 'order' can take values 1 or 2, returning the
% first or second order derivative, respectively.
%
% The first-order method utilizes a two-sided two-point scheme to calculate
% intermediate points, while one-sided forward and backward two-point
% schemes is used to calculate the first and last data point, respectively.
%
% The second-order method utilizes a two-sided three-point scheme to
% calculate intermediate points, while one-sided forward and backward
% three-point schemes is used to calculate the first and last data point,
% respectively.
% -------------------------------------------------------------------------
% Written by Torstein E. Daehlin, August 2019
% -------------------------------------------------------------------------

% Detect and assigne variable input arguments
nargs = nargin;
if nargs == 3
    order = varagin(1);
else
    order = 1;
end

% Determine size of input
[row, col] = size(y);

% Transpose input if there are more columns than rows
transpose = false;
if col > row
    y = y';
    [row, col] = size(y);
    transpose = true;
end

% Preallocate output array
dydx = zeros(row, col);

% Select first or second order case
if order == 1
    for cdx = 1:col
        % Differentiate first point using one-sided forward two-point scheme
        dydx(1,cdx) = (y(2,cdx) - y(1,cdx))/h;

        % Differentiate 2:n-1 points using two-sided two-point scheme
        for idx = 2:(row-1)
            dydx(idx,cdx) = (y(idx+1,cdx) - y(idx-1,cdx))/(2*h);
        end

        % Differentiate point n using one-sided backward two-point scheme
        dydx(idx+1,cdx) = (y(idx+1,cdx) - y(idx,cdx))/h;
    end
elseif order == 2
    for cdx = 1:col
        % Differentiate first point using one-sided forward three-point scheme
        dydx(1,cdx) = (y(3,cdx) - 2*y(2,cdx) + y(1,cdx))/h^2;

        % Differentiate 2:n-1 points using two-sided three-point scheme
        for idx = 2:(row-1)
            dydx(idx,cdx) = (y(idx+1,cdx) - 2*y(idx,cdx) + y(idx-1,cdx))/h^2;
        end

        % Differentiate point n using one-sided backward two-point scheme
        dydx(idx+1,cdx) = (y(idx+1,cdx) - 2*y(idx,cdx) + y(idx-1,cdx))/h^2;
    end
else
    error('The input argument ''order'' must take values 1 or 2');
end

% If input data series had time in columns, transpose the output
if transpose
    dydx = dydx';
end

end

function [segment_angles, joint_angles, angular_velocity, angular_acceleration] = FindAngularKinematics(dynamic_lcs, static_lcs, nof, time)
%

% Calculate segment angles
% ------------------------
% Extract segment names
segment_names = fieldnames(dynamic_lcs);

% Loop over segments
for i = 1:length(segment_names)
    % Loop over frames
    for frame = 1:nof
        % Define rotation matrix
        R = [dynamic_lcs.(segment_names{i}).epx(frame,:)' dynamic_lcs.(segment_names{i}).epy(frame,:)' ...
            dynamic_lcs.(segment_names{i}).epz(frame,:)'];
        [x_angle, y_angle, z_angle] = EulerAngles(R,'zyx','deg');

        % Assign output
        segment_angles.(segment_names{i})(frame,:) = [x_angle y_angle z_angle];
    end
end

% Calculate joint angles
% ----------------------
% Define joint names
joint_names = {'hip','knee','ankle','midfoot'};
side = {'right','left'};

for i = 1:length(joint_names)
    for j = 1:length(side)
        if isequal(side{j},'right')
            post_script = '_r';
        elseif isequal(side{j},'left')
            post_script = '_l';
        else
            error(['Invalid side: ' side{j}]);
        end
        for frame = 1:nof
            % Define proximal and distal coordinate systems
            switch joint_names{i}
                case 'hip'
                    R_prox = [dynamic_lcs.pelvis.epx(frame,:)' dynamic_lcs.pelvis.epy(frame,:)' ...
                        dynamic_lcs.pelvis.epz(frame,:)'];
                    R_prox_static = [static_lcs.pelvis.epx' static_lcs.pelvis.epy' ...
                        static_lcs.pelvis.epz'];
                    R_dist = [dynamic_lcs.(['thigh' post_script]).epx(frame,:)' dynamic_lcs.(['thigh' post_script]).epy(frame,:)' ...
                        dynamic_lcs.(['thigh' post_script]).epz(frame,:)'];
                    R_dist_static = [static_lcs.(['thigh' post_script]).epx' static_lcs.(['thigh' post_script]).epy' ...
                        static_lcs.(['thigh' post_script]).epz'];
                case 'knee'
                    R_prox = [dynamic_lcs.(['thigh' post_script]).epx(frame,:)' dynamic_lcs.(['thigh' post_script]).epy(frame,:)' ...
                        dynamic_lcs.(['thigh' post_script]).epz(frame,:)'];
                    R_prox_static = [static_lcs.(['thigh' post_script]).epx' static_lcs.(['thigh' post_script]).epy' ...
                        static_lcs.(['thigh' post_script]).epz'];
                    R_dist = [dynamic_lcs.(['leg' post_script]).epx(frame,:)' dynamic_lcs.(['leg' post_script]).epy(frame,:)' ...
                        dynamic_lcs.(['leg' post_script]).epz(frame,:)'];
                    R_dist_static = [static_lcs.(['leg' post_script]).epx' static_lcs.(['leg' post_script]).epy' ...
                        static_lcs.(['leg' post_script]).epz'];
                case 'ankle'
                    R_prox = [dynamic_lcs.(['leg' post_script]).epx(frame,:)' dynamic_lcs.(['leg' post_script]).epy(frame,:)' ...
                        dynamic_lcs.(['leg' post_script]).epz(frame,:)'];
                    R_prox_static = [static_lcs.(['leg' post_script]).epx' static_lcs.(['leg' post_script]).epy' ...
                        static_lcs.(['leg' post_script]).epz'];
                    R_dist = [dynamic_lcs.(['foot' post_script]).epx(frame,:)' dynamic_lcs.(['foot' post_script]).epy(frame,:)' ...
                        dynamic_lcs.(['foot' post_script]).epz(frame,:)'];
                    R_dist_static = [static_lcs.(['foot' post_script]).epx' static_lcs.(['foot' post_script]).epy' ...
                        static_lcs.(['foot' post_script]).epz'];
                case 'midfoot'
                    R_prox = [dynamic_lcs.(['rearfoot' post_script]).epx(frame,:)' dynamic_lcs.(['rearfoot' post_script]).epy(frame,:)' ...
                        dynamic_lcs.(['rearfoot' post_script]).epz(frame,:)'];
                    R_prox_static = [static_lcs.(['rearfoot' post_script]).epx' static_lcs.(['rearfoot' post_script]).epy' ...
                        static_lcs.(['rearfoot' post_script]).epz'];
                    R_dist = [dynamic_lcs.(['forefoot' post_script]).epx(frame,:)' dynamic_lcs.(['forefoot' post_script]).epy(frame,:)' ...
                        dynamic_lcs.(['forefoot' post_script]).epz(frame,:)'];
                    R_dist_static = [static_lcs.(['forefoot' post_script]).epx' static_lcs.(['forefoot' post_script]).epy' ...
                        static_lcs.(['forefoot' post_script]).epz'];
                otherwise
                    error(['Invalid joint: ' joint_names{i}]);
            end

            % Calculate joint angles
            R = (R_prox_static'*R_prox)'*(R_dist_static'*R_dist);
            [x_angle(frame,:), y_angle(frame,:), z_angle(frame,:)] = EulerAngles(R,'xyz','deg');
        end

        % Define output
        joint_angles.(joint_names{i}).(side{j}) = [x_angle y_angle z_angle];
    end
end

% Calculate angular velocities and accelerations
% ----------------------------------------------
% Segment angular velocities
for i = 1:length(segment_names)
    % Find Euler rates
    euler_rate = FiniteDiff(segment_angles.(segment_names{i}), (time(2)-time(1)), 1);

    % Calculate angular velocities
    for frame = 1:nof
        % Define trasformation matrix
        E = [cosd(segment_angles.(segment_names{i})(frame,2))*cosd(segment_angles.(segment_names{i})(frame,3)) -sind(segment_angles.(segment_names{i})(frame,3)) 0;...
            cosd(segment_angles.(segment_names{i})(frame,2))*sind(segment_angles.(segment_names{i})(frame,3)) cosd(segment_angles.(segment_names{i})(frame,3)) 0;...
            sind(segment_angles.(segment_names{i})(frame,2)) 0 1];
        angular_velocity.(segment_names{i})(frame,:) = (E*euler_rate(frame,:)')';
    end

    % Calculate angular acceleration
    angular_acceleration.(segment_names{i}) = FiniteDiff(angular_velocity.(segment_names{i}), time(2)-time(1), 1);
end

% Joint angular velocities
for i = 1:length(joint_names)
    for j = 1:length(side)
        % Find Euler rates
        euler_rate = FiniteDiff(joint_angles.(joint_names{i}).(side{j}), time(2)-time(1), 1);

        % Calculate angular velocities
        for frame = 1:nof
            E = [1 0 -sind(joint_angles.(joint_names{i}).(side{j})(frame,2)); ...
                0 cosd(joint_angles.(joint_names{i}).(side{j})(frame,1)) -sind(joint_angles.(joint_names{i}).(side{j})(frame,1))*cosd(joint_angles.(joint_names{i}).(side{j})(frame,2)); ...
                0 sind(joint_angles.(joint_names{i}).(side{j})(frame,1)) cosd(joint_angles.(joint_names{i}).(side{j})(frame,1))*cosd(joint_angles.(joint_names{i}).(side{j})(frame,2))];
            angular_velocity.(joint_names{i}).(side{j})(frame,:) = (E*euler_rate(frame,:)')';
        end

        % Calculate angular acceleration
        angular_acceleration.(joint_names{i}).(side{j}) = FiniteDiff(angular_velocity.(joint_names{i}).(side{j}), time(2)-time(1), 1);
    end
end

end

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

function njm = CalculatejointMoments(segments, dynamic_lcs, jc, grf, angular_velocity, angular_acceleration, acceleration, position, nof)

% Define gravity vector
g = [0 0 -9.81];

% Extract segment and joint names
sides = {'right','left'};
joint_names = fieldnames(jc);

% Loop over frames
for frame = 1:nof
    for i = 1:length(joint_names)
        for side = 1:length(sides)
            % Find names of segments located distal to the joint
            segment_names = FindSegmentNames([joint_names{i} '_' sides{side}]);

            % Preallocate
            tau_I = zeros(length(segment_names),3);
            tau = zeros(length(segment_names),3);

            % Loop over segments distal to the joint
            for j = 1:length(segment_names)

                % Calculate net force acting on each segment
                F = segments.(segment_names{j}).mass * (acceleration.(segment_names{j})(frame,:) - g);

                % Transform angular velocity and accelereation vectors to the
                % segment coordinate system
                R = [dynamic_lcs.(segment_names{j}).epx(frame,:)' dynamic_lcs.(segment_names{j}).epy(frame,:)' ...
                    dynamic_lcs.(segment_names{j}).epz(frame,:)'];
                omega = (R'*deg2rad(angular_velocity.(segment_names{j})(frame,:))')';
                alpha = (R'*deg2rad(angular_acceleration.(segment_names{j})(frame,:))')';

                % Calculate inretial moment of the segment and transform it
                % back to the global coordiante system
                tau_I(j,:) = (R * (segments.(segment_names{j}).tensor*alpha' + ...
                    cross(omega',(segments.(segment_names{j}).tensor*omega'))))';

                % Find moment arm of centre of mass
                r = position.(segment_names{j})(frame,:) - jc.(joint_names{i}).(sides{side})(frame,:);

                % Calculate segment moment
                tau(j,:) = cross(r,F);
            end

            % Find moment arm of ground reaction force
            r_grf = grf.(sides{side}).cop(frame,:) - jc.(joint_names{i}).(sides{side})(frame,:);

            % Calculate net joint moment
            njm_global = sum((tau_I + tau),1) - ...
                grf.(sides{side}).free_moment(frame,:) - cross(r_grf,grf.(sides{side}).force(frame,:));

            % Transform net joint moment into the coordiante system of the
            % distal segment
            R = [dynamic_lcs.(segment_names{1}).epx(frame,:)' dynamic_lcs.(segment_names{1}).epy(frame,:)' ...
                dynamic_lcs.(segment_names{1}).epz(frame,:)'];
            njm.(joint_names{i}).(sides{side})(frame,:) = (R*njm_global')';
        end
    end
end
end

function joint_power = CalculateJointPower(njm, angular_velocity, dynamic_lcs, nof)

% Get joint names and side names
joint_names = fieldnames(njm);
side_names = fieldnames(njm.(joint_names{1}));

for joint = 1:length(joint_names)
    for side = 1:length(side_names)
        % Find segment names (1st element of the vector returned is the segment distal to the joint)
        segment_names = FindSegmentNames([joint_names{joint} '_' side_names{side}]);

        for frame = 1:nof
            % Transform joint angular velocity into local coordinate system of the distal segment (njm is already expressed here)
            R = [dynamic_lcs.(segment_names{1}).epx(frame,:)' dynamic_lcs.(segment_names{1}).epy(frame,:)' ...
                dynamic_lcs.(segment_names{1}).epz(frame,:)'];
            omega = (R'*deg2rad(angular_velocity.(joint_names{joint}).(side_names{side})(frame,:))')';

            % Calculate instantaneous power
            joint_power.(joint_names{joint}).(side_names{side})(frame,:) =  ...
                njm.(joint_names{joint}).(side_names{side})(frame,:).*omega;
        end
    end
end
end

function PlotDynamic(markers, lcs, jc, grf, trial_name, visit_name, nof)

% Plot parameters
axis_scale = 0.1;
line_colors = {'r-','g-','b-'};

% Define figure
fig = figure('Name',[trial_name ' - ' visit_name]);
fig.WindowState = 'maximized';
pause(2);

% Loop over frames nad animate motion
for frame = 1:nof

    % Plot global coordinate system
    global_ax = eye(3);
    for j = 1:3
        plot3([0 global_ax(1,j)*axis_scale*2],[0 global_ax(2,j)*axis_scale*2],[0 global_ax(3,j)*axis_scale*2],line_colors{j});
        hold on;
    end

    % Label global coordinate axes
    text(0.21,0,0,'X','Color','r','FontWeight','bold');
    text(0,0.21,0,'Y','Color','g','FontWeight','bold');
    text(0,0,0.21,'Z','Color','b','FontWeight','bold');

    % Plot markers
    marker_names = fieldnames(markers);
    for i = 1:length(marker_names)
        % Select marker color
        if marker_names{i}(1) == 'R'
            marker_color = '#77AC30';
        elseif marker_names{i}(1) == 'L'
            marker_color = '#0072BD';
        else
            marker_color = '#D95319';
        end

        % Plot marker
        plot3(markers.(marker_names{i})(frame,1),markers.(marker_names{i})(frame,2),markers.(marker_names{i})(frame,3), ...
            'o','MarkerEdgeColor',marker_color,'MarkerFaceColor',marker_color);
    end

    % Plot local coordinate systems
    lcs_names = fieldnames(lcs);
    axis_names = {'epx','epy','epz'};
    for i = 1:length(lcs_names)
        for j = 1:length(axis_names)
            plot3([lcs.(lcs_names{i}).origin(frame,1) lcs.(lcs_names{i}).origin(frame,1)+lcs.(lcs_names{i}).(axis_names{j})(frame,1)*axis_scale], ...
                [lcs.(lcs_names{i}).origin(frame,2) lcs.(lcs_names{i}).origin(frame,2)+lcs.(lcs_names{i}).(axis_names{j})(frame,2)*axis_scale], ...
                [lcs.(lcs_names{i}).origin(frame,3) lcs.(lcs_names{i}).origin(frame,3)+lcs.(lcs_names{i}).(axis_names{j})(frame,3)*axis_scale], ...
                line_colors{j});
        end
    end

    % Plot joint centres
    jc_names = {'hip','knee','ankle'};
    sides = {'right','left'};
    for i = 1:length(jc_names)
        for j = 1:length(sides)
            plot3(jc.(jc_names{i}).(sides{j})(frame,1),jc.(jc_names{i}).(sides{j})(frame,2),jc.(jc_names{i}).(sides{j})(frame,3), ...
                'o','MarkerEdgeColor','#A2142F','MarkerFaceColor','#A2142F','MarkerSize',10);
        end
    end

    % Plot force platforms and ground reaction force
    for j = 1:length(sides)
        for i = 1:3
            plot3([grf.(sides{j}).corners(i,1) grf.(sides{j}).corners(i+1,1)], ...
                [grf.(sides{j}).corners(i,2) grf.(sides{j}).corners(i+1,2)], ...
                [grf.(sides{j}).corners(i,3) grf.(sides{j}).corners(i+1,3)],'w-');
        end
        plot3([grf.(sides{j}).corners(1,1) grf.(sides{j}).corners(4,1)], ...
            [grf.(sides{j}).corners(1,2) grf.(sides{j}).corners(4,2)], ...
            [grf.(sides{j}).corners(1,3) grf.(sides{j}).corners(4,3)],'w-');
        plot3([grf.(sides{j}).cop(frame,1) grf.(sides{j}).cop(frame,1)+grf.(sides{j}).force(frame,1)/1000], ...
            [grf.(sides{j}).cop(frame,2) grf.(sides{j}).cop(frame,2)+grf.(sides{j}).force(frame,2)/1000], ...
            [grf.(sides{j}).cop(frame,3) grf.(sides{j}).cop(frame,3)+grf.(sides{j}).force(frame,3)/1000],'r-');
    end

    % Format axes
    anchor = 0.5*(markers.RCREST(end,:) + markers.LCREST(end,:));
    ax = gca;
    ax.Color = 'k';
    ax.XLim = [anchor(1)-1 anchor(1)+1];
    ax.YLim = [anchor(2)-1 anchor(2)+1];
    ax.ZLim = [-0.1 1.6];
    ax.DataAspectRatio = [1 1 1];
    ax.View = [-235 20];
    ax.XLabel.String = 'Position (m)';
    ax.YLabel.String = 'Position (m)';
    ax.ZLabel.String = 'Position (m)';

    pause(0.005);
    if ~isequal(frame,nof)
        clf(fig);
    end
end
end

function segment_names = FindSegmentNames(joint_name)

switch joint_name
    case 'hip_right'
        segment_names = {'thigh_r','leg_r','foot_r'};
    case 'hip_left'
        segment_names = {'thigh_l','leg_l','foot_l'};
    case 'knee_right'
        segment_names = {'leg_r','foot_r'};
    case 'knee_left'
        segment_names = {'leg_l','foot_l'};
    case 'ankle_right'
        segment_names = {'foot_r'};
    case 'ankle_left'
        segment_names = {'foot_l'};
    otherwise
        error(['Invalid joint name:' joint_name]);
end
end

function events = FindEvents(kinematics, kinetics, jump_type)


if ismember(jump_type,{'AJL','AJR'})

    % Define parameters
    side_names = {'right','left'};
    threshold_grf = 20; % N
    threshold_knee_vel = 30;
    dly = 20;

    % Loop over sides and extract events
    for side = 1:length(side_names)
        idx(1) = find(kinetics.grf.(side_names{side}).force(1:end,3) > threshold_grf, 1);
        idx(2) = find(kinetics.grf.(side_names{side}).force(idx(1)+dly:end,3) < threshold_grf, 1);
        idx(3) = find(kinetics.grf.(side_names{side}).force(idx(1)+dly+idx(2):end,3) > threshold_grf, 1);
        idx(4) = find(kinematics.angular_velocity.knee.(side_names{side})(idx(1)+dly+idx(2)+idx(3):end,1) <= threshold_knee_vel & ...
            kinematics.angular_velocity.knee.(side_names{side})(idx(1)+dly+idx(2)+idx(3):end,1) >= -threshold_knee_vel, 1);
        events.(side_names{side}).start = idx(1);
        events.(side_names{side}).foot_off = idx(1)+dly+idx(2);
        events.(side_names{side}).foot_ic = idx(1)+dly+idx(2)+idx(3);
        events.(side_names{side}).max_knee_flex = idx(1)+dly+idx(2)+idx(3)+idx(4);
    end

elseif isequal(jump_type,'VJ')

    % Define parameters
    side_names = {'right','left'};
    threshold_grf = 10; % N
    threshold_knee_vel = 30;
    threshold_pelvis_vel = -0.1;
    dly = 60;

    % Loop over sides and extract events
    for side = 1:length(side_names)
        idx(1) = find(kinematics.velocity.pelvis(:,3) < threshold_pelvis_vel, 1);
        idx(2) = find(kinetics.grf.(side_names{side}).force(idx(1)+dly:end,3) < threshold_grf, 1);
        idx(3) = find(kinetics.grf.(side_names{side}).force(idx(1)+dly+idx(2):end,3) > threshold_grf, 1);
        idx(4) = find(kinematics.angular_velocity.knee.(side_names{side})(idx(1)+dly+idx(2)+idx(3):end,1) < threshold_knee_vel & ...
            kinematics.angular_velocity.knee.(side_names{side})(idx(1)+dly+idx(2)+idx(3):end,1) > -threshold_knee_vel, 1);
        events.(side_names{side}).start = idx(1);
        events.(side_names{side}).foot_off = idx(1)+dly+idx(2);
        events.(side_names{side}).foot_ic = idx(1)+dly+idx(2)+idx(3);
        events.(side_names{side}).max_knee_flex = idx(1)+dly+idx(2)+idx(3)+idx(4);

    end
else
    error(['Invalid jump type: ' jump_type]);
end

end

function [output, jump_height] = ExtractOutput(static_lcs, segments, kinematics, kinetics, events, old_jump_height, old_output, time)

% Find jump height
R = [static_lcs.pelvis.epx' static_lcs.pelvis.epy' static_lcs.pelvis.epz'];
height_0 = R*segments.pelvis.com' + static_lcs.pelvis.origin';
height_max = max(kinematics.position.pelvis(:,3));
jump_height = height_max - height_0(3,1);

if jump_height > old_jump_height

    % calculate time step
    h = time(2)-time(1);

    % Find net joint work and average net joint moment
    joint_names = fieldnames(kinetics.joint_power);
    side_names = fieldnames(kinetics.joint_power.(joint_names{1}));
    for side = 1:length(side_names)
        for joint = 1:length(joint_names)
            % Extract propulsion phase variables of interest
            work.(joint_names{joint}).(side_names{side}).propulsion = ...
                trapz(h,kinetics.joint_power.(joint_names{joint}).(side_names{side})(events.(side_names{side}).start:events.(side_names{side}).foot_off,1));
            moment.(joint_names{joint}).(side_names{side}).propulsion = ...
                mean(kinetics.njm.(joint_names{joint}).(side_names{side})(events.(side_names{side}).start:events.(side_names{side}).foot_off,1));
            excursion.(joint_names{joint}).(side_names{side}).propulsion = ...
                kinematics.joint_angles.(joint_names{joint}).(side_names{side})(events.(side_names{side}).foot_off,1) - ...
                kinematics.joint_angles.(joint_names{joint}).(side_names{side})(events.(side_names{side}).start,1);
            arch_deformation.(side_names{side}).propulsion = max(kinematics.joint_angles.midfoot.(side_names{side})(events.(side_names{side}).start:events.(side_names{side}).foot_off,1));

            % Extract landing phase variables of interest
            work.(joint_names{joint}).(side_names{side}).landing = ...
                trapz(h,kinetics.joint_power.(joint_names{joint}).(side_names{side})(events.(side_names{side}).foot_ic:events.(side_names{side}).max_knee_flex,1));
            moment.(joint_names{joint}).(side_names{side}).landing = ...
                mean(kinetics.njm.(joint_names{joint}).(side_names{side})(events.(side_names{side}).foot_ic:events.(side_names{side}).max_knee_flex,1));
            excursion.(joint_names{joint}).(side_names{side}).landing = ...
                kinematics.joint_angles.(joint_names{joint}).(side_names{side})(events.(side_names{side}).max_knee_flex,1) - ...
                kinematics.joint_angles.(joint_names{joint}).(side_names{side})(events.(side_names{side}).foot_ic,1);
            arch_deformation.(side_names{side}).landing = max(kinematics.joint_angles.midfoot.(side_names{side})(events.(side_names{side}).foot_ic:events.(side_names{side}).max_knee_flex,1));


            % Define start time
            if events.right.start < events.left.start
                start_frame = events.right.start;
            else
                start_frame = events.left.start;
            end

            % Define end time
            if events.right.max_knee_flex > events.left.max_knee_flex
                end_frame =  events.right.max_knee_flex;
            else
                end_frame =  events.left.max_knee_flex;
            end

            % Extract time series
            time_series.moment.(joint_names{joint}).(side_names{side}) = kinetics.njm.(joint_names{joint}).(side_names{side})(start_frame:end_frame,:);
            time_series.power.(joint_names{joint}).(side_names{side}) = kinetics.joint_power.(joint_names{joint}).(side_names{side})(start_frame:end_frame,:);
            time_series.angles.(joint_names{joint}).(side_names{side}) = kinematics.joint_angles.(joint_names{joint}).(side_names{side})(start_frame:end_frame,:);

        end
        % Extract centre of pressure location
        pressure.(side_names{side}).propulsion = mean(kinetics.relative_cop.(side_names{side})(events.(side_names{side}).start:events.(side_names{side}).foot_off,:));
        pressure.(side_names{side}).landing = mean(kinetics.relative_cop.(side_names{side})(events.(side_names{side}).foot_ic:events.(side_names{side}).max_knee_flex,:));

        % Extract ground reaction force time series
        time_series.grf.(side_names{side}) = kinetics.grf.(side_names{side}).force;
    end

    % Plot events for verification
    PlotEvents(kinetics, kinematics, events, time);

    % Assign output
    output = struct('jump_height',jump_height,'work',work,'moment',moment,'excursion',excursion, ...
        'arch_deformation',arch_deformation,'pressure',pressure,'time_series',time_series);
else
    output = old_output;
end
end

function PlotEvents(kinetics, kinematics, events, time)

% Define time vector
sides = sort(fieldnames(events));
event_names = fieldnames(events.right);
line_specs = {'b-','r-'};
labels = {'on','off','ic','pkf'};

% Create figure
figure();

for s = 1:length(sides)
    % plot ground reaction force in subplot 1
    subplot(2,2,0+s);
    plot(time,kinetics.grf.(sides{s}).force(:,3),line_specs{s},'Linewidth',1.5);
    hold on;
    ax = gca;

    for e = 1:length(event_names)
        plot([time(events.(sides{s}).(event_names{e})) time(events.(sides{s}).(event_names{e}))], ...
            ax.YLim,'k--');
        text(time(events.(sides{s}).(event_names{e})),ax.YLim(2)-50,labels{e},'HorizontalAlignment','center');
    end

    ylabel('Force (N)');
    xlabel('Time (s)');
    xlim([0 time(end)]);
    title(['Vertical ground reaction force ' sides{s}]);

    % plot knee angle in subplot 2
    subplot(2,2,2+s);
    plot(time,kinematics.joint_angles.knee.(sides{s})(:,1),line_specs{s},'Linewidth',1.5);
    hold on;
    ax = gca;
    for e = 1:length(event_names)
        plot([time(events.(sides{s}).(event_names{e})) time(events.(sides{s}).(event_names{e}))], ...
            ax.YLim,'k--');
        text(time(events.(sides{s}).(event_names{e})),ax.YLim(2)-5,labels{e},'HorizontalAlignment','center');
    end
    ylabel('Angle (deg)');
    xlabel('Time (s)');
    xlim([0 time(end)]);
    title(['Knee flexion/extension angle ' sides{s}]);
end


end

function relative_cop = FindRelativeCop(dynamic_markers, dynamic_lcs, segments, grf, nof)

sides = fieldnames(grf);

for side = 1:length(sides)
    origo = 0.5*(dynamic_markers.([upper(sides{side}(1)) 'MTH1']) + dynamic_markers.([upper(sides{side}(1)) 'MTH5']));
    vec = grf.(sides{side}).cop - origo;

    for frame = 1:nof
        if isequal(lower(sides{side}),'left')
            relative_cop.(sides{side})(frame,1) = (dot(vec(frame,:),[dynamic_lcs.(['foot_' sides{side}(1)]).epx(frame,1:2) 0])/ ...
                segments.(['foot_' sides{side}(1)]).dist_rad) * -100;
        else
            relative_cop.(sides{side})(frame,1) = (dot(vec(frame,:),[dynamic_lcs.(['foot_' sides{side}(1)]).epx(frame,1:2) 0])/ ...
                segments.(['foot_' sides{side}(1)]).dist_rad) * 100;
        end
        relative_cop.(sides{side})(frame,2) = (dot(vec(frame,:),[dynamic_lcs.(['foot_' sides{side}(1)]).epz(frame,1:2) 0])/ ...
            segments.(['foot_' sides{side}(1)]).len) * 100;
    end
end

end
