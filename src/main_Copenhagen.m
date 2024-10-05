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
                InverseKinematics(dynamic(j).markers, static(i).markers, ...
                static_lcs, static_jc, segments, meta.dynamic(j));            

            % Determine which external forces are applied to which segments
            roi = IdentifyROI(time, dynamic(j).force(2).force, subj.mass);


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
