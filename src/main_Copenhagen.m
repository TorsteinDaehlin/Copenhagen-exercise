% main_Copenhagen.m
% -------------------------------------------------------------------------
% Performes inverse dynamics batch processing of Copenhagen exercise data
% collected using a Qualisys motion capture system, 2 AMTI force platforms,
% and a CAST marker set consisting of 28 anatomical markers and 27 tracking
% markers for the incline vs. block heel-raise research project conducted
% in the Sports Biomechanics Laboratory at the University of Alberta, AB,
% Canada 2020-2021.
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
addpath(['.' filesep 'Initialize' filesep]);
addpath(['.' filesep 'Preproccess' filesep]);
addpath(['.' filesep 'ModelDefinition' filesep]);
addpath(['.' filesep 'InverseDynamics' filesep]);
addpath(['.' filesep 'InverseKinematics' filesep]);
addpath(['.' filesep 'HelperFunctions' filesep]);
addpath(['.' filesep 'Postprocess' filesep]);

% Prompt user to select participant directory
[subj_dir, src_path, start_idx] = GetSubjDir();

% Look for marker registration file in participant directory. If not there,
% prompt user to open marker registration file
if isfile(fullfile(subj_dir(1).folder,'model_setup.csv'))
    marker_reg = readtable(fullfile(subj_dir(1).folder, 'model_setup.csv'));
else
    [reg_file, reg_path] = uigetfile('..\..\*.csv');
    marker_reg = readtable(fullfile(reg_path, reg_file));
end

% Look for participant characteristics file in participant directory
if isfile(fullfile(subj_dir(1).folder, 'Copenhagen I - Participant Characteristics.xlsx'))
    subj_char = readtable(fullfile(subj_dir(1).folder, 'Copenhagen I - Participant Characteristics.xlsx'), 'NumHeaderLines', 1);
    subj_list = subj_char.Code;
else
    subj_char = [];
    subj_list = {};
end

% Create output directory
dst_path = regexp(src_path, filesep, 'split');
dst_path = fullfile(dst_path{1:end-1},'Output');

if ~isfolder(dst_path)
    mkdir(dst_path);
end

% Load previous results if they already exist
if isfile(fullfile(dst_path, 'CPH_results.mat'))
    load(fullfile(dst_path, 'CPH_results.mat'), 'tbls');
else
    tbls = struct('C_A', [], 'C_B', [], 'C_C', []);
end

% Requrest processing input (e.g. filter parameters)
[static_id, flt] = GetProcessingParameters();

if isempty(static_id) || isempty(flt)
    fprintf('Please fill provide all processing parameters\n');
    [static_id, flt] = GetProcessingParameters();
end

try
    for s = 1:length(subj_dir)
        % Inspect participant directory
        % =============================
        % Get subject characteristics
        if ~isempty(subj_char)
            subj.id = subj_char.Code{(start_idx + s) - 1};
            subj.height = subj_char.Height((start_idx + s) - 1);
            subj.mass = subj_char.BodyMass((start_idx + s) - 1);
        else
            subj.id = subj_dir(s).name;

            % Prompt user to provide height and mass
            prompt = {'Enter participant height (m):', ...
                'Enter participant mass (kg):'};
            answer = inputdlg(prompt, [subj.id ': Height and mass']);
            subj.height = str2double(answer{1});
            subj.mass = str2double(answer{2});

            subj_list{end} = subj.id;
        end

        % Create paths for subject outputs
        subj.out_path = fullfile(dst_path, subj.id);
        if ~isfolder(subj.out_path)
            mkdir(subj.out_path);
        end

        subj.check_path = fullfile(subj.out_path, 'DataChecks');
        if ~isfolder(subj.check_path)
            mkdir(subj.check_path);
        end

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
            [static_lcs, static_jc, segments, joints] = ...
                ProcessStatic(static(i), meta.static(i), subj, i);

            % Loop over dynamic trials matched to current static
            for j = 1:length(static.match_to_move)

                % Perform inverse kinematics
                [kinematics, time] = ...
                    InverseKinematics(dynamic(j).markers, static(i).markers, ...
                    static_lcs, static_jc, segments, meta.dynamic(j));

                % Determine which external forces are applied to which segments
                roi = IdentifyROI(time, dynamic(j).force(2).force, kinematics.position.thigh_r(:,3), subj, static(i).match_to_move(j));
                grf_act_on = ApplyForceToSegment(kinematics.jc, dynamic(j).force(2).cop, roi);

                % Visualize dynamic trial
                PlotDynamic(dynamic(j).markers, kinematics, dynamic(j).force(2), ...
                    roi, subj, static(i).match_to_move(j));

                %  Error check grf_act_on
                if ~isequal(grf_act_on{:})
                    warning('COP segment mismatch');

                    % Open figure
                    chck_figs = dir(fullfile(subj.check_path, ['*_dynamic_' num2str(static(i).match_to_move(j)) '_*.jpg']));
                    img = figure('WindowState','maximized');
                    for k = 1:length(chck_figs)
                        subplot(1, length(chck_figs), k);
                        I = imread(fullfile(chck_figs(k).folder, chck_figs(k).name));
                        imshow(I);
                    end

                    % Prompt user to select appropriate segment
                    answer = listdlg('SelectionMode','single','PromptString', 'Select segment GRF acts on', ...
                        'ListString', grf_act_on);
                    grf_act_on = grf_act_on{answer};

                    % Close figure
                    close(img);
                else
                    grf_act_on = grf_act_on{1};
                end

                % Calculate NJMs using inverse dynamics
                njm = ...
                    calcJointMoments(segments, kinematics, dynamic(j).force(2), ...
                    joints, grf_act_on, meta.dynamic(j).nof);

                % Plot NJM data checks
                PlotKinematicsChecks(time, kinematics, roi, subj, static(i).match_to_move(j));
                PlotNjmChecks(time, njm, dynamic(j).force(2), roi, subj, static(i).match_to_move(j));

                % Store outputs
                R.(strtok(subj.move_name{static(i).match_to_move(j)}, '.'))(s) = ...
                    PostprocessIvd(kinematics, dynamic(j).force(2), njm, subj, roi, meta.dynamic(j));

                % Store time series
                ts(s).(strtok(subj.move_name{static(i).match_to_move(j)}, '.')) = ...
                    PostprocessTimeSeries(time, kinematics, dynamic(j).force(2), ...
                    njm, roi);
                ts(s).id = subj.id;
            end
        end
    end
    % Export data
    ExportResults(tbls, R, ts, dst_path, subj_list);

catch
    % Print warning message
    warning('Error occured while processing data from participant %s. Exiting.', subj.id);

    % Export results to current point
    ExportResults(tbls, R, ts, dst_path, subj_list);
end
end