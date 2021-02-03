function addEvent(msg,append)
% Find main window to holds the logger 
main = findobj('Tag','AlignData_main_window');

% Add msg to the logger and display
if ~exist('append','var')
    append = 0;
end
   
if append
    cur_string = main.UserData.logger.String{main.UserData.logger.Value};
    main.UserData.logger.String{main.UserData.logger.Value} = [cur_string,msg];
else
    main.UserData.logger.Value = main.UserData.logger.Value + 1;
    if isempty(main.UserData.logger.String)
        main.UserData.logger.String = {msg};
    else
        main.UserData.logger.String = [main.UserData.logger.String;msg];
    end
end
end