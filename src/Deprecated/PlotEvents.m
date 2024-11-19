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