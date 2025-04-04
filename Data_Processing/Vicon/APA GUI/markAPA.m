function markAPA()
screen_size = get(0,'ScreenSize');
if screen_size(3) < 1800
    button_font_size = 14;
else
    button_font_size = 20;
end

% Main figure window
main_window = figure('Visible','on','units','normalized','outerposition',[0.10,0.15,0.8,0.7],...
    'MenuBar','none',...
    'DockControls','off',...
    'NumberTitle','off',...
    'Name','Gait Initiation APA Marking',...
    'color',[128,128,128]./256,...
    'UserData',struct(),...
    'Tag','main_window');

% Raw force and TTL pulse figure
force_ttl_plot = axes('Parent',main_window,'Units','normalized','Position',[0.03,0.34,0.35,0.62]);
title(force_ttl_plot,'Raw force and Delsys TTL Pulse','FontSize',24);
cop_marking_plot = axes('Parent',main_window,'Units','normalized','Position',[0.41,0.34,0.35,0.62]);
title(cop_marking_plot,'CoP APA Marking','FontSize',24);

cop_apa_example_plot = axes('Parent',main_window,'Units','normalized','Position',[0.79,0.32,0.2,0.65],'Visible','off');
image=imread('apa_cop_example.jpg'); % How to read in jpg image
imshow(image,'Border','tight','Parent',cop_apa_example_plot) % How to plot with a specific axes handel
title(cop_apa_example_plot,'Landmark Guide','FontSize',24);

set(cop_apa_example_plot,'xtick',[],'xticklabels',[],'ytick',[],'yticklabels',[],...
    'xcolor','k','ycolor','k','linewidth',1.5);
set(force_ttl_plot,'ytick',[],'yticklabels',[],'xcolor','k','ycolor','k','linewidth',1.5);
set(cop_marking_plot,'xcolor','k','ycolor','k','linewidth',1.5);

try
    box([force_ttl_plot,cop_marking_plot,cop_apa_example_plot],'on');
catch
    box(force_ttl_plot,'On');
    box(cop_marking_plot,'On');
    box(cop_apa_example_plot,'On');
end

% Event panel
event_panel = uicontrol('Parent',main_window,...
    'Style','listbox',...
    'BackgroundColor','white',...
    'Units','normalized',...
    'Value',0,...
    'Position',[0.01,0.01,0.98,0.16],...
    'Enable','inactive');

% Button to load vicon data
 uicontrol('Parent',main_window,...
    'Style','pushbutton',...
    'BackgroundColor','white',...
    'String','Load aligned data',...
    'FontSize',16,...
    'Units','normalized',...
    'Position',[0.03 0.24 0.11 0.05],...
    'Callback',@loadData);

% Drop down menu to select plot line 1 and 2
 uicontrol('Parent',main_window,...
    'Style','text',...
    'String','FP Right:',...
    'BackgroundColor',[128,128,128]./256,...
    'HorizontalAlignment','center',...
    'FontSize',14,...
    'Units','normalized',...
    'Position',[0.15 0.27 0.05 0.02]);
 uicontrol('Parent',main_window,...
    'Callback',@updateForcePlot,...
    'Style','popupmenu',...
    'String',' ',...
    'Fontsize',12,...
    'BackgroundColor','white',...
    'Units','normalized',...
    'Position',[0.15 0.215 0.05 0.05],...
    'Tag','Plot_Line_FP1',...
    'Enable','inactive');

 uicontrol('Parent',main_window,...
    'Style','text',...
    'String','FP Left:',...
    'BackgroundColor',[128,128,128]./256,...
    'HorizontalAlignment','center',...
    'FontSize',14,...
    'Units','normalized',...
    'Position',[0.205 0.27 0.05 0.02]);
 uicontrol('Parent',main_window,...
    'Callback',@updateForcePlot,...
    'Style','popupmenu',...
    'String',' ',...
    'Fontsize',12,...
    'BackgroundColor','white',...
    'Units','normalized',...
    'Position',[0.205 0.215 0.05 0.05],...
    'Tag','Plot_Line_FP2',...
    'Enable','inactive');

% Button to zoom
 uicontrol('Parent',main_window,...
    'Style','togglebutton',...
    'BackgroundColor','white',...
    'String','Zoom',...
    'FontSize',16,...
    'Units','normalized',...
    'Position',[0.27 0.24 0.11 0.05],...
    'Callback',@setZoom);

% Button to set Cue
 uicontrol('Parent',main_window,...
    'Style','togglebutton',...
    'BackgroundColor','white',...
    'String','Mark cue',...
    'FontSize',20,...
    'Units','normalized',...
    'Position',[0.03 0.18 0.11 0.05],...
    'Callback',@markMovement);

% Button to set movement onset
 uicontrol('Parent',main_window,...
    'Style','togglebutton',...
    'BackgroundColor','white',...
    'String','Mark movement onset',...
    'FontSize',button_font_size,...
    'Units','normalized',...
    'Position',[0.15 0.18 0.11 0.05],...
    'Callback',@markMovement);

% Button to set initial stance leg toe-off
 uicontrol('Parent',main_window,...
    'Style','togglebutton',...
    'BackgroundColor','white',...
    'String','Mark stance leg toe-off',...
    'FontSize',button_font_size,...
    'Units','normalized',...
    'Position',[0.27 0.18 0.11 0.05],...
    'Callback',@markMovement);

% Button to plot CoP
 uicontrol('Parent',main_window,...
    'Style','pushbutton',...
    'BackgroundColor','white',...
    'String','Plot CoP',...
    'FontSize',20,...
    'Units','normalized',...
    'Position',[0.41 0.24 0.35 0.05],...
    'Callback',@plotCoP);

% Button to mark landmark 1
 uicontrol('Parent',main_window,...
    'Style','togglebutton',...
    'BackgroundColor','white',...
    'String','Mark Landmark 1',...
    'FontSize',20,...
    'Units','normalized',...
    'Position',[0.41 0.18 0.17 0.05],...
    'Callback',@markLandmark);

% Button to mark landmark 2
 uicontrol('Parent',main_window,...
    'Style','togglebutton',...
    'BackgroundColor','white',...
    'String','Mark Landmark 2',...
    'FontSize',20,...
    'Units','normalized',...
    'Position',[0.59 0.18 0.17 0.05],...
    'Callback',@markLandmark);

% Button to calculate APA metrics
uicontrol('Parent',main_window,...
    'Style','pushbutton',...
    'BackgroundColor','white',...
    'String','Calc APA Metrics',...
    'FontSize',20,...
    'Units','normalized',...
    'Position',[0.815 0.24 0.15 0.05],...
    'Callback',@calcAPAMetrics);

% Button to save data
 uicontrol('Parent',main_window,...
    'Style','pushbutton',...
    'BackgroundColor','white',...
    'String','Save markings',...
    'FontSize',20,...
    'Units','normalized',...
    'Position',[0.815 0.18 0.15 0.05],...
    'Callback',@saveData);

main_window.UserData.logger = event_panel;
main_window.UserData.plot_axes = [force_ttl_plot,cop_marking_plot,cop_apa_example_plot];
main_window.UserData.plot_markers = nan(1,5);
main_window.UserData.Vicon_data = [];
main_window.UserData.Vicon_time = [];
main_window.UserData.Delsys_data = [];
main_window.UserData.Marked_events = nan(1,5);
main_window.UserData.save_path = [];
main_window.UserData.data_file = [];
end