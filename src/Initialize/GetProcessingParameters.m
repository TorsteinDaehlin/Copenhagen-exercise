function GetProcessingParameters()

% Create GUI figure
fig = uifigure('Name','Processing parameters','Visible','off', ...
    'HandleVisibility','on');
fig.Position(3:4) = [506 327];

% Get static identifier. The script will ask for a generic identifier that
% is used to identify which trial in the participant directory that
% contains the static trial. The identifier should be a substring that is
% generic to all particiants. E.g., if P[id]_static.mat is used as
% identifier, "static" could be provided as generic identifier
static_instruction = 'Provide generic identifier for static trials';
static_pnl = uipanel(fig, 'Title', 'Static:', 'Position',[20 243 462 64], ...
'FontSize',12,'FontWeight','bold');

static_lbl = uilabel(static_pnl, 'HorizontalAlignment', 'right', 'Position', ...
    [11 9 215 22], 'Text', static_instruction);
static_field = uieditfield(static_pnl, 'text', 'Position',[241 9 203 22]);

% Get desired filter parameters
filter_pnl = uipanel(fig, 'Title', 'Filter', 'Position', [20 64 462 152], ...
    'FontSize',12,'FontWeight','bold');

fc_lbl = uilabel(filter_pnl, 'HorizontalAlignment', 'right', 'Position', ...
    [142 82 192 22], 'Text', 'Desired filter cut-off frequency (Hz)');
fc_field = uieditfield(filter_pnl, 'numeric', 'Position', [349 82 100 22]);

order_lbl = uilabel(filter_pnl, 'HorizontalAlignment', 'right', 'Position', ...
    [230 37 104 22], 'Text', 'Desired filter order');
order_field = uieditfield(filter_pnl, 'numeric', 'Position', ...
    [349 37 100 22]);


for i = 1:length(T)
    uitable(fig,'Data',T(i).tbl,'ColumnEditable',[false false false],'Position',tbl_pos(i,:));
    lbl(i) = uilabel(fig,'Text',nm{i},'Position',lbl_pos(i,:),'FontWeight','bold');
end
uiimage(fig,'ImageSource',icon,'Position',icon_pos);
uilabel(fig,'Text',msg,'Position',msg_pos,'WordWrap',true,'FontSize',16);
pnl = uipanel(fig,'Title','Instructions:','Position',[40 70 220 99], ...
    'FontSize',12,'FontWeight','bold');
uilabel(pnl,'Text',instruction,'WordWrap',true,'FontSize',12,'Position',instruction_pos);
lb = uilistbox(fig, 'Items', exp_lbls, 'Position',[40 201 220 99]);
btn1 = uibutton(fig,'Text','Continue','Position',[230 25 100 22]);
btn1.ButtonPushedFcn = @(src,~)uiresume(ancestor(src,'figure'));
fig.Visible = 'on';
uiwait(fig);
new_lbl = lb.Value;
close(fig);


end



            
           

            % Create TypeButtonGroup
            app.TypeButtonGroup = uibuttongroup(app.FilterPanel);
            app.TypeButtonGroup.Title = 'Type';
            app.TypeButtonGroup.Position = [10 8 123 116];

            % Create ButterworthButton
            app.ButterworthButton = uiradiobutton(app.TypeButtonGroup);
            app.ButterworthButton.Text = 'Butterworth';
            app.ButterworthButton.Position = [11 70 83 22];
            app.ButterworthButton.Value = true;

            % Create ChebyshevIButton
            app.ChebyshevIButton = uiradiobutton(app.TypeButtonGroup);
            app.ChebyshevIButton.Text = 'Chebyshev I';
            app.ChebyshevIButton.Position = [11 48 88 22];

            % Create ChebyshevIIButton
            app.ChebyshevIIButton = uiradiobutton(app.TypeButtonGroup);
            app.ChebyshevIIButton.Text = 'Chebyshev II';
            app.ChebyshevIIButton.Position = [11 26 92 22];

            % Create EllipticButton
            app.EllipticButton = uiradiobutton(app.TypeButtonGroup);
            app.EllipticButton.Text = 'Elliptic';
            app.EllipticButton.Position = [11 3 56 22];

            % Create ContinueButton
            app.ContinueButton = uibutton(app.UIFigure, 'push');
            app.ContinueButton.Position = [123 20 100 22];
            app.ContinueButton.Text = 'Continue';

            % Create CancelButton
            app.CancelButton = uibutton(app.UIFigure, 'push');
            app.CancelButton.Position = [305 20 100 22];
            app.CancelButton.Text = 'Cancel';

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = app1

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end