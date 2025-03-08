% main_Copenhagen.m
% =========================================================================
% MIT License
% 
% Copyright (c) 2025 Torstein Eriksen Dæhlin, PhD
% 
% Permission is hereby granted, free of charge, to any person obtaining a
% copy of this software and associated documentation files (the
% "Software"), to deal in the Software without restriction, including
% without limitation the rights to use, copy, modify, merge, publish,
% distribute, sublicense, and/or sell copies of the Software, and to permit
% persons to whom the Software is furnished to do so, subject to the
% following conditions:
% 
% The above copyright notice and this permission notice shall be included
% in all copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
% MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
% NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
% DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
% OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
% USE OR OTHER DEALINGS IN THE SOFTWARE.
% 
% -------------------------------------------------------------------------
% This software performs inverse dynamics batch processing of Copenhagen
% exercise motion capture data collected. The software may be modified for
% analysis of different tasks, marker sets, and motion capture setups.
% -------------------------------------------------------------------------
% Written by Torstein E. Daehlin, PhD.
% -------------------------------------------------------------------------

% References:
% -----------
%{
Cappozzo A, Catani F, Della Croce U, Leardini A. (1995) Position and
    orientation in space of bones during movement: anatomical frame definition
    and determination. Clin Biomech. 10(17), pp 1–8.

Cappozzo, A., Cappello, A., Croce, U. D., & Pensalfini, F. (1997)
    Surface-marker cluster design criteria for 3-D bone movement
    reconstruction. IEEE Transactions on Biomedical Engineering, 44(12),
    1165-1174.

Della Croce U, Cappozzo A, Kerrigan DC. (1999) Pelvis and lower limb
    anatomical landmark calibration precision and its propaga- tion to bone
    geometry and joint angles. Med Biol Eng Comp. 37, pp 155–161.

Dempster, W. T. (1955) Space requirements of the seated operator,
    geometrical, kinematic, and mechanical aspects of the body with special
    reference to the limbs. University of Michigan.

Grood, E. S., & Suntay, W. J. (1983) A joint coordinate system for the
    clinical description of three-dimensional motions: Application to the knee.
    Journal of Biomechanical Engineering, 105(2), 136–144.

Söderkvist, I., & Wedin, P. Å. (1993) Determining the movements of the
    skeleton using well-configured markers. Journal of biomechanics, 26(12),
    1473-1477.

Wu G, Siegler S, Allard P, Kirtley C, Leardini A, Rosenbaum D, et al.
    (2002) ISB recommendation on definitions of joint coordinate system of
    various joints for the reporting of human joint motion. Part 1: ankle, hip,
    and spine. J Biomech. 35, pp 543–548.

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

% We use a marker registration file called 'model_setup.csv' to allow more
% flexible marker labelling. The program looks for this marker registration
% file in the participant directory by default. If it is not found, the
% user is prompted to open the marker registration file
if isfile(fullfile(subj_dir(1).folder,'model_setup.csv'))
    marker_reg = readtable(fullfile(subj_dir(1).folder, 'model_setup.csv'));
else
    [reg_file, reg_path] = uigetfile('..\..\*.csv');
    marker_reg = readtable(fullfile(reg_path, reg_file));
end

% Participant characteristics are provided in a participant characteristics
% file. This version uses a .xlsx file containing containing the column
% "Code" containing participant identifiers, "Height" containing
% participant heights, and "BodyMass" containing participant mass. If a
% characteristics file is not provided, the folder name of the folder
% containing participant data is used as participant name, and the user is
% prompted to provided participant height and mass. 
if isfile(fullfile(subj_dir(1).folder, 'Copenhagen I - Participant Characteristics.xlsx'))
    subj_char = readtable(fullfile(subj_dir(1).folder, 'Copenhagen I - Participant Characteristics.xlsx'), 'NumHeaderLines', 1);
    subj_list = subj_char.Code;
else
    subj_char = [];
    subj_list = {};
end

% Files are saved to this output directory
dst_path = regexp(src_path, filesep, 'split');
dst_path = fullfile(dst_path{1:end-1},'Output');

if ~isfolder(dst_path)
    mkdir(dst_path);
end

% If previous results have been saved, these are reloaded before starting
% the processing loop
if isfile(fullfile(dst_path, 'CPH_results.mat'))
    load(fullfile(dst_path, 'CPH_results.mat'), 'tbls');
else
    tbls = struct('C_A', [], 'C_B', [], 'C_C', []);
end

% The user is prompted to provide some filter parameters and a common name
% to identify static trials by
[static_id, flt] = GetProcessingParameters();

if isempty(static_id) || isempty(flt)
    fprintf('Please fill provide all processing parameters\n');
    [static_id, flt] = GetProcessingParameters();
end

try
    for s = 1:length(subj_dir)
        % We store the subject characteristics in the subj structure
        if ~isempty(subj_char)
            subj.id = subj_char.Code{(start_idx + s) - 1};
            subj.height = subj_char.Height((start_idx + s) - 1);
            subj.mass = subj_char.BodyMass((start_idx + s) - 1);
        else
            subj.id = subj_dir(s).name;

            prompt = {'Enter participant height (m):', ...
                'Enter participant mass (kg):'};
            answer = inputdlg(prompt, [subj.id ': Height and mass']);
            subj.height = str2double(answer{1});
            subj.mass = str2double(answer{2});

            subj_list{end} = subj.id;
        end

        % Participant specific outputs are stored here
        subj.out_path = fullfile(dst_path, subj.id);
        if ~isfolder(subj.out_path)
            mkdir(subj.out_path);
        end

        subj.check_path = fullfile(subj.out_path, 'DataChecks');
        if ~isfolder(subj.check_path)
            mkdir(subj.check_path);
        end

        motion_files = dir(fullfile(subj_dir(s).folder, subj_dir(s).name));
        motion_files = motion_files(~(strcmp({motion_files.name}, {'.'}) | ...
            strcmp({motion_files.name}, {'..'})));

        static_idx = contains({motion_files.name}, static_id);
        if ~any(static_idx)
            error('No static trial whose file name contains "%s".', static_id);
        end

        subj.data_path = motion_files(1).folder;
        subj.static_name = {motion_files(static_idx).name};
        subj.move_name = {motion_files(~static_idx).name};

        % Preprocess data
        % ===============
        % Markers are forces are loaded and restructured into the format
        % used for the remainder of the program. See further comments
        % inside this function for supported data formats
        [static, dynamic, meta] = PreprocessMOCAP(subj, marker_reg, flt);

        % We transform the recorded ground reactions force to the top of
        % stand by summing the moments about its top centre. See manuscript
        % for details. PS! this should be removed if no implement is placed
        % on the force platform.
        dynamic = TransformToStand(dynamic);

        % Run inverse dynamics procedure
        % ==============================
        for i = 1:length(static)
            % Generate participant model
            % --------------------------
            % We generate one model per static recording in the participant
            % folder. See PreprocessMOCAP() for details on how multiple
            % static trials are handled.
            [static_lcs, static_jc, segments, joints] = ...
                ProcessStatic(static(i), meta.static(i), subj, i);

            for j = 1:length(static.match_to_move)

                % Perform inverse kinematics
                [kinematics, time] = ...
                    InverseKinematics(dynamic(j).markers, static(i).markers, ...
                    static_lcs, static_jc, segments, meta.dynamic(j));

                % We attempt to automatically detect the regions of
                % interest and apply the stand reaction force to the
                % appropriate segment. Plots are then procduced and the
                % user is asked to verify ROIs and segemnet allocation.
                roi = IdentifyROI(time, dynamic(j).force(2).force, kinematics.position.thigh_r(:,3), subj, static(i).match_to_move(j));
                grf_act_on = ApplyForceToSegment(kinematics.jc, dynamic(j).force(2).cop, roi);

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

                % Newton-Euler iterative inverse dynamic are used to
                % compute net joint moments which are expressed in the
                % coordinate system of the distal ("child") segment of each
                % joint. 
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