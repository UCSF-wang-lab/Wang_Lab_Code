function AlignDataGUI()
%% Check to make sure it is in the correct path to run GUI
curr_path = mfilename('fullpath');
slashes = strfind(curr_path,'/');
if ~strcmp(curr_path(1:slashes(end)-1),pwd)
    cd(curr_path(1:slashes(end)));
end

%% Initialize GUI
% Create windowed figure
main_window = figure('Visible','on','units','normalized','outerposition',[0.10,0.15,0.8,0.7],...
    'MenuBar','none',...
    'DockControls','off',...
    'NumberTitle','off',...
    'Name','RC+S Data Alignment',...
    'color',[128,128,128]./256,...
    'UserData',struct(),...
    'Tag','AlignData_main_window');

% Data panel
createDataSourcePanel(main_window);

% Plot panels
createPlotPanel(main_window);

% Plotting option panel
createPlotOptionPanel(main_window);

% Alignment tool panel
createAlignmentToolsPanel(main_window);

% Event panel
event_panel = uicontrol('Parent',main_window,...
    'Style','listbox',...
    'BackgroundColor','white',...
    'Units','normalized',...
    'Value',0,...
    'Position',[0.01,0.01,0.98,0.18],...
    'Enable','inactive');

% Create struct of user data. This is so we can pass around data between
% the GUI.
main_window.UserData.LFP_data = [];
main_window.UserData.Accel_data = [];
main_window.UserData.DeviceSettings = [];
main_window.UserData.LogTable = [];
main_window.UserData.Delsys_data = [];
main_window.UserData.Xsens_data = [];
main_window.UserData.FP_data = [];
main_window.UserData.Rover_data = [];
main_window.UserData.stim_condition = 'OFF';
main_window.UserData.med_condition = 'ON';
main_window.UserData.trial_num = [];
main_window.UserData.file_names = {};
main_window.UserData.alignment_times = nan(1,9); % Left LFP|Right LFP|Left Accel|Right Accel|Delsys|Xsens|Force plate|Rover Left|Rover Right
main_window.UserData.alignment_source = cell(1,9); % Left LFP|Right LFP|Left Accel|Right Accel|Delsys|Xsens|Force plate|Rover Left|Rover Right
main_window.UserData.logger = event_panel;
main_window.UserData.basePath = [];
end