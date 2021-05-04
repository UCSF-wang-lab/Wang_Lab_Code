function createAlignmentToolsPanel(main_window)
%% Alignment tools
alignment_tools_panel = uipanel('Title','Alignment tools:',...
    'FontSize',24,...
    'BackgroundColor','white',...
    'BorderType','none',...
    'Units','normalized',...
    'Position',[0.01 0.20 0.2 0.29]);

uicontrol('Parent',alignment_tools_panel,...
    'Style','togglebutton',...
    'BackgroundColor','white',...
    'String','Mark alignment points',...
    'FontSize',16,...
    'Units','normalized',...
    'Position',[0.01 0.84 0.98 0.15],...
    'Tag','Left_LFP_button',...
    'Callback',@markData);

uicontrol('Parent',alignment_tools_panel,...
    'Style','pushbutton',...
    'BackgroundColor','white',...
    'String','Align data',...
    'FontSize',16,...
    'Units','normalized',...
    'Position',[0.01 0.71 0.485 0.12],...
    'Tag','Left_LFP_button',...
    'Callback',@alignData);

uicontrol('Parent',alignment_tools_panel,...
    'Style','text',...
    'BackgroundColor','white',...
    'String',{'Pre-alignment'; 'time in sec:'},...
    'FontSize',12,...
    'Units','normalized',...
    'Position',[0.505 0.715 0.20 0.12],...
    'Tag','Left_LFP_button');

uicontrol('Parent',alignment_tools_panel,...
    'Style','edit',...
    'BackgroundColor','white',...
    'String','0',...
    'FontSize',16,...
    'Units','normalized',...
    'Position',[0.705 0.71 0.285 0.12],...
    'Tag','Left_LFP_button',...
    'Callback',@changePreAlignmentTime);

uicontrol('Parent',alignment_tools_panel,...
    'Style','pushbutton',...
    'BackgroundColor','white',...
    'String','Show alignment',...
    'FontSize',16,...
    'Units','normalized',...
    'Position',[0.01 0.58 0.98 0.10],...
    'Tag','Left_LFP_button',...
    'Callback',@showAlignment);

uicontrol('Parent',alignment_tools_panel,...
    'Style','pushbutton',...
    'BackgroundColor','white',...
    'String','Save alignment',...
    'FontSize',16,...
    'Units','normalized',...
    'Position',[0.01 0.47 0.98 0.10],...
    'Tag','Left_LFP_button',...
    'Callback',@saveData);

% Alignment time that have been marked
alignment_times_text(1) = uicontrol('Parent',alignment_tools_panel,...
    'Style','text',...
    'HorizontalAlignment','left',...
    'BackgroundColor','white',...
    'String',sprintf('Left RC+S = %i',nan),...
    'FontSize',14,...
    'Units','normalized',...
    'Position',[0.01 0.38 0.48 0.08]);

alignment_times_text(2) = uicontrol('Parent',alignment_tools_panel,...
    'Style','text',...
    'HorizontalAlignment','left',...
    'BackgroundColor','white',...
    'String',sprintf('Right RC+S = %i',nan),...
    'FontSize',14,...
    'Units','normalized',...
    'Position',[0.01 0.29 0.48 0.08]);

alignment_times_text(3) = uicontrol('Parent',alignment_tools_panel,...
    'Style','text',...
    'HorizontalAlignment','left',...
    'BackgroundColor','white',...
    'String',sprintf('Left RC+S Accel = %i',nan),...
    'FontSize',14,...
    'Units','normalized',...
    'Position',[0.01 0.20 0.48 0.08]);

alignment_times_text(4) = uicontrol('Parent',alignment_tools_panel,...
    'Style','text',...
    'HorizontalAlignment','left',...
    'BackgroundColor','white',...
    'String',sprintf('Right RC+S Accel = %i',nan),...
    'FontSize',14,...
    'Units','normalized',...
    'Position',[0.01 0.11 0.48 0.08]);

alignment_times_text(5) = uicontrol('Parent',alignment_tools_panel,...
    'Style','text',...
    'HorizontalAlignment','left',...
    'BackgroundColor','white',...
    'String',sprintf('Delsys = %i',nan),...
    'FontSize',14,...
    'Units','normalized',...
    'Position',[0.51 0.38 0.48 0.08]);

alignment_times_text(6) = uicontrol('Parent',alignment_tools_panel,...
    'Style','text',...
    'HorizontalAlignment','left',...
    'BackgroundColor','white',...
    'String',sprintf('Xsens = %i',nan),...
    'FontSize',14,...
    'Units','normalized',...
    'Position',[0.51 0.29 0.48 0.08]);

alignment_times_text(7) = uicontrol('Parent',alignment_tools_panel,...
    'Style','text',...
    'HorizontalAlignment','left',...
    'BackgroundColor','white',...
    'String',sprintf('Force Plate = %i',nan),...
    'FontSize',14,...
    'Units','normalized',...
    'Position',[0.51 0.20 0.48 0.08]);

alignment_times_text(8) = uicontrol('Parent',alignment_tools_panel,...
    'Style','text',...
    'HorizontalAlignment','left',...
    'BackgroundColor','white',...
    'String',sprintf('Teensey = %i',nan),...
    'FontSize',14,...
    'Units','normalized',...
    'Position',[0.51 0.11 0.48 0.08]);


% Save variables
main_window.UserData.alignment_times_text = alignment_times_text;
main_window.UserData.pre_alignment_time = 0;
    
end