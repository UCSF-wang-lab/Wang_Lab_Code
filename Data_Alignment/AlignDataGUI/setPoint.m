function txt = setPoint(~,cursor_info)
% Grabs and saves the cursor position depending on the button that was
% pressed.
%
% Author:   Kenneth Louie
% Date:     12/28/20

main_window = findobj('Tag','AlignData_main_window');

plot_window_name = cursor_info.Target.Parent.Title.String;
colon_ind = strfind(plot_window_name,':');
data_source_name = plot_window_name(colon_ind+2:end);
if contains(plot_window_name,'Left') && contains(plot_window_name,'LFP')
    main_window.UserData.alignment_times(1) = cursor_info.Position(1);
    main_window.UserData.alignment_times_text(1).String = sprintf('Left RC+S = %4.3f',main_window.UserData.alignment_times(1));
    main_window.UserData.alignment_source{1} = data_source_name;
elseif contains(plot_window_name,'Right') && contains(plot_window_name,'LFP')
    main_window.UserData.alignment_times(2) = cursor_info.Position(1);
    main_window.UserData.alignment_times_text(2).String = sprintf('Right RC+S = %4.3f',main_window.UserData.alignment_times(2));
    main_window.UserData.alignment_source{2} = data_source_name;
elseif contains(plot_window_name,'Left') && contains(plot_window_name,'Accel')
    main_window.UserData.alignment_times(3) = cursor_info.Position(1);
    main_window.UserData.alignment_times_text(3).String = sprintf('Left RC+S Accel = %4.3f',main_window.UserData.alignment_times(3));
    main_window.UserData.alignment_source{3} = data_source_name;
elseif contains(plot_window_name,'Right') && contains(plot_window_name,'Accel')
    main_window.UserData.alignment_times(4) = cursor_info.Position(1);
    main_window.UserData.alignment_times_text(4).String = sprintf('Right RC+S Accel = %4.3f',main_window.UserData.alignment_times(4));
    main_window.UserData.alignment_source{4} = data_source_name;
elseif contains(plot_window_name,'Delsys')
    main_window.UserData.alignment_times(5) = cursor_info.Position(1);
    main_window.UserData.alignment_times_text(5).String = sprintf('Delsys = %4.3f',main_window.UserData.alignment_times(5));
    main_window.UserData.alignment_source{5} = data_source_name;
elseif contains(plot_window_name,'Xsens')
    main_window.UserData.alignment_times(6) = cursor_info.Position(1);
    main_window.UserData.alignment_times_text(6).String = sprintf('Xsens = %4.3f',main_window.UserData.alignment_times(6));
    main_window.UserData.alignment_source{6} = data_source_name;
elseif contains(plot_window_name,'Force')
    main_window.UserData.alignment_times(7) = cursor_info.Position(1);
    main_window.UserData.alignment_times_text(7).String = sprintf('Force plate = %4.3f',main_window.UserData.alignment_times(7));
    main_window.UserData.alignment_source{7} = data_source_name;
end

txt = {sprintf('t = %4.3f',cursor_info.Position(1));...
    sprintf('y = %4.3f',cursor_info.Position(2))};
end