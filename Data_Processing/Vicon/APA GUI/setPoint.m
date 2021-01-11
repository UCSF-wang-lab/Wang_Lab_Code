function txt = setPoint(~,cursor_info,button)
% Grabs and saves the cursor position depending on the button that was
% pressed.
%
% Author:   Kenneth Louie
% Date:     12/28/20

main_window = findobj('Tag','main_window');

if contains(button.String,'cue')
    main_window.UserData.Marked_events(1) = cursor_info.Position(1);
    
    output_string = sprintf('Cue time: t = %3.6f, f = %3.6f',cursor_info.Position(1),cursor_info.Position(2));
    main_window.UserData.logger.String = [main_window.UserData.logger.String;output_string];
    main_window.UserData.logger.Value = length(main_window.UserData.logger.String);
    
elseif contains(button.String,'onset')
    main_window.UserData.Marked_events(2) = cursor_info.Position(1);
    
    output_string = sprintf('Movement onset time: t = %3.6f, f = %3.6f',cursor_info.Position(1),cursor_info.Position(2));
    main_window.UserData.logger.String = [main_window.UserData.logger.String;output_string];
    main_window.UserData.logger.Value = length(main_window.UserData.logger.String);
    
elseif contains(button.String,'toe-off')
    main_window.UserData.Marked_events(3) = cursor_info.Position(1);
    
    output_string = sprintf('Stance leg toe-off time: t = %3.6f, f = %3.6f',cursor_info.Position(1),cursor_info.Position(2));
    main_window.UserData.logger.String = [main_window.UserData.logger.String;output_string];
    main_window.UserData.logger.Value = length(main_window.UserData.logger.String);
    
elseif contains(button.String,'Landmark 1')
    main_window.UserData.Marked_events(4) = getCursorInfo(button.Parent.UserData.dataCursorObj).DataIndex;
    
    main_window.UserData.logger.String = [main_window.UserData.logger.String;'Landmark 1 set.'];
    main_window.UserData.logger.Value = length(main_window.UserData.logger.String);
    
elseif contains(button.String,'Landmark 2')
    main_window.UserData.Marked_events(5) = getCursorInfo(button.Parent.UserData.dataCursorObj).DataIndex;
    
    main_window.UserData.logger.String = [main_window.UserData.logger.String;'Landmark 2 set.'];
    main_window.UserData.logger.Value = length(main_window.UserData.logger.String);
end

if contains(button.String,'Landmark')
    txt = {sprintf('x = %3.6f',cursor_info.Position(1));...
        sprintf('y = %3.6f',cursor_info.Position(2))};
else
    txt = {sprintf('t = %3.6f',cursor_info.Position(1));...
        sprintf('f = %3.6f',cursor_info.Position(2))};
end
end