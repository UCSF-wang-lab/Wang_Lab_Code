function change_screen_color(src,event,trial_type,trial_state)
clear_screen();

% Get main_window object
try
    main_window = src.Parent;
catch
    main_window = findobj('Tag','GI_Task_main_window');
end

% play sound according to trial type
play_sound(trial_type);

% Change screen color
switch trial_state
    case 1
        main_window.Color = main_window.UserData.colors.stop;
    case 2
        main_window.Color = main_window.UserData.colors.prepare;
    case 3
        main_window.Color = main_window.UserData.colors.go;
    case 4
        main_window.Color = [128,128,128]./256;
        pause(2);
        trial_option_screen(main_window);
end

% Write change to save file
write_input(main_window,trial_type,trial_state,toc(main_window.UserData.clock));

% Start time to change color if beginning of trial
if trial_state == 1 && trial_type == 3
    start(main_window.UserData.timer_prepare);
    start(main_window.UserData.timer_cue);
elseif trial_state == 1
    start(main_window.UserData.timer_cue);
elseif trial_state == 3
    main_window.UserData.timer_stop.StartDelay = 5;
    main_window.UserData.timer_stop.TimerFcn = {@change_screen_color,trial_type,4};
    start(main_window.UserData.timer_stop);
end

end