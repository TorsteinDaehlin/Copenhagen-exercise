function [static_id, flt] = GetProcessingParameters()

% Create GUI figure
fig = uifigure('Name','Processing parameters','Visible','off', ...
    'HandleVisibility','on');
fig.Position(3:4) = [520 327];

% Get static identifier. The script will ask for a generic identifier that
% is used to identify which trial in the participant directory that
% contains the static trial. The identifier should be a substring that is
% generic to all particiants. E.g., if P[id]_static.mat is used as
% identifier, "static" could be provided as generic identifier
static_instruction = 'Provide generic identifier for static trials';
static_pnl = uipanel(fig, 'Title', 'Static:', 'Position',[20 243 462 64], ...
'FontSize',12,'FontWeight','bold');

uilabel(static_pnl, 'HorizontalAlignment', 'right', 'Position', ...
    [11 9 215 22], 'Text', static_instruction);
static_field = uieditfield(static_pnl, 'text', 'Position',[241 9 203 22]);

% Get desired filter parameters
filter_pnl = uipanel(fig, 'Title', 'Filter', 'Position', [20 64 462 152], ...
    'FontSize',12,'FontWeight','bold');

uilabel(filter_pnl, 'HorizontalAlignment', 'right', 'Position', ...
    [142 82 192 22], 'Text', 'Desired filter cut-off frequency (Hz)');
fc_field = uieditfield(filter_pnl, 'numeric', 'Position', [356 82 100 22]);

uilabel(filter_pnl, 'HorizontalAlignment', 'right', 'Position', ...
    [230 37 104 22], 'Text', 'Desired filter order');
order_field = uieditfield(filter_pnl, 'numeric', 'Position', ...
    [356 37 100 22]);

uilabel(filter_pnl, 'HorizontalAlignment', 'right', 'Position', ...
    [10 96 31 22], 'Text', 'Type');
type_list = uilistbox(filter_pnl, 'Items', ...
    {'Butterworth', 'Chebyshev I', 'Chebyshev II', 'Ellipse'}, ...
    'Position', [51 13 92 107], ...
    'Value', 'Butterworth', 'Multiselect','off');

cont_btn = uibutton(fig, 'push', 'Position', [204,23,100,22], ...
    'Text', 'Continue');
cont_btn.ButtonPushedFcn = @(src,~)uiresume(ancestor(src,'figure'));

fig.Visible = 'on';
uiwait(fig);

% Get values after button press
if ~isempty(static_field.Value)
    static_id = static_field.Value;
else
    static_id = [];
end

% Error check inputs
if ~any([fc_field.Value == 0, order_field.Value == 0, ...
        isempty(type_list.Value)])
    flt.fc = fc_field.Value;
    flt.order = order_field.Value;
    flt.type = type_list.Value;
else
    flt = [];
end

close(fig);
end