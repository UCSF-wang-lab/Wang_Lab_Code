function txt = setPoint(~,cursor_info)
% Grabs and saves the cursor position depending on the button that was
% pressed.
%
% Author:   Kenneth Louie
% Date:     05/04/21

main_window = findobj('Tag','AlignData_main_window');
rcs_time_menu = findobj('Tag','RCSTimeUnit');

plot_window_name = cursor_info.Target.Parent.Title.String;
colon_ind = strfind(plot_window_name,':');
data_source_name = plot_window_name(colon_ind+2:end);
if contains(plot_window_name,'Left') && contains(plot_window_name,'LFP')
    if strcmp(rcs_time_menu.String{rcs_time_menu.Value},'0 Time Start')  
        main_window.UserData.alignment_times(1) = cursor_info.Position(1);
        main_window.UserData.alignment_times_text(1).String = sprintf('Left RC+S = %7.4f',main_window.UserData.alignment_times(1));
    else
        main_window.UserData.alignment_times(1) = cursor_info.Target.Children.DataIndex;  % Use the data index instead of the data tip x position because it gives it in releative units
        main_window.UserData.alignment_times_text(1).String = sprintf('Left RC+S = %s',string(timeofday(cursor_info.Target.XData(main_window.UserData.alignment_times(1))),'hh:mm:ss.SSS'));
    end
    main_window.UserData.alignment_source{1} = data_source_name;
elseif contains(plot_window_name,'Right') && contains(plot_window_name,'LFP')
    if strcmp(rcs_time_menu.String{rcs_time_menu.Value},'0 Time Start')  
        main_window.UserData.alignment_times(2) = cursor_info.Position(1);
        main_window.UserData.alignment_times_text(2).String = sprintf('Right RC+S = %7.4f',main_window.UserData.alignment_times(2));
    else
        main_window.UserData.alignment_times(2) = cursor_info.Target.Children.DataIndex;  % Use the data index instead of the data tip x position because it gives it in releative units
        main_window.UserData.alignment_times_text(2).String = sprintf('Right RC+S = %s',string(timeofday(cursor_info.Target.XData(main_window.UserData.alignment_times(2))),'hh:mm:ss.SSS'));
    end
    main_window.UserData.alignment_source{2} = data_source_name;
elseif contains(plot_window_name,'Left') && contains(plot_window_name,'Accel') && ~contains(plot_window_name,'Linear')
    if strcmp(rcs_time_menu.String{rcs_time_menu.Value},'0 Time Start')  
        main_window.UserData.alignment_times(3) = cursor_info.Position(1);
        main_window.UserData.alignment_times_text(3).String = sprintf('Left RC+S Accel = %7.4f',main_window.UserData.alignment_times(3));
    else
        main_window.UserData.alignment_times(3) = cursor_info.Target.Children.DataIndex;  % Use the data index instead of the data tip x position because it gives it in releative units
        main_window.UserData.alignment_times_text(3).String = sprintf('Left RC+S Accel = %s',string(timeofday(cursor_info.Target.XData(main_window.UserData.alignment_times(3))),'hh:mm:ss.SSS'));
    end
    main_window.UserData.alignment_source{3} = data_source_name;
elseif contains(plot_window_name,'Right') && contains(plot_window_name,'Accel') && ~contains(plot_window_name,'Linear')
    if strcmp(rcs_time_menu.String{rcs_time_menu.Value},'0 Time Start')  
        main_window.UserData.alignment_times(4) = cursor_info.Position(1);
        main_window.UserData.alignment_times_text(4).String = sprintf('Right RC+S Accel = %7.4f',main_window.UserData.alignment_times(4));
    else
        main_window.UserData.alignment_times(4) = cursor_info.Target.Children.DataIndex;  % Use the data index instead of the data tip x position because it gives it in releative units
        main_window.UserData.alignment_times_text(4).String = sprintf('Right RC+S Accel = %s',string(timeofday(cursor_info.Target.XData(main_window.UserData.alignment_times(4))),'hh:mm:ss.SSS'));
    end
    main_window.UserData.alignment_source{4} = data_source_name;
elseif contains(plot_window_name,'Delsys')
    main_window.UserData.alignment_times(5) = cursor_info.Position(1);
    main_window.UserData.alignment_times_text(5).String = sprintf('Delsys = %7.4f',main_window.UserData.alignment_times(5));
    main_window.UserData.alignment_source{5} = data_source_name;
elseif contains(plot_window_name,'Xsens')
    main_window.UserData.alignment_times(6) = cursor_info.Position(1);
    main_window.UserData.alignment_times_text(6).String = sprintf('Xsens = %7.4f',main_window.UserData.alignment_times(6));
    main_window.UserData.alignment_source{6} = data_source_name;
elseif contains(plot_window_name,'Force')
    main_window.UserData.alignment_times(7) = cursor_info.Position(1);
    main_window.UserData.alignment_times_text(7).String = sprintf('Force plate = %7.4f',main_window.UserData.alignment_times(7));
    main_window.UserData.alignment_source{7} = data_source_name;
elseif contains(plot_window_name,'Left Rover')
    main_window.UserData.alignment_times(8) = cursor_info.Target.Children.DataIndex;  % Use the data index instead of the data tip x position because it gives it in releative units
    main_window.UserData.alignment_times_text(8).String = sprintf('Left Rover = %s',string(cursor_info.Target.XData(main_window.UserData.alignment_times(8)),'hh:mm:ss.SSS'));
    main_window.UserData.alignment_source{8} = data_source_name;
elseif contains(plot_window_name,'Right Rover')
    main_window.UserData.alignment_times(9) = cursor_info.Target.Children.DataIndex;  % Use the data index instead of the data tip x position because it gives it in releative units
    main_window.UserData.alignment_times_text(9).String = sprintf('Right Rover = %s',string(cursor_info.Target.XData(main_window.UserData.alignment_times(9)),'hh:mm:ss.SSS'));
    main_window.UserData.alignment_source{9} = data_source_name;
end

if ~contains(plot_window_name,'Rover')
    txt = {sprintf('t = %7.4f',cursor_info.Position(1));...
        sprintf('y = %7.4f',cursor_info.Position(2))};
else
    txt = {sprintf('t = %s',string(cursor_info.Target.XData(cursor_info.Target.Children.DataIndex),'hh:mm:ss.SSS'));...
        sprintf('y = %7.4f',cursor_info.Position(2))};
end

end