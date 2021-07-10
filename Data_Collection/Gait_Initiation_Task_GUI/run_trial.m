function run_trial(src,event)
%%
% run_trial(...)
%
% Author: Kenneth Louie
% Project: General Gait Protocol
% Date: 05/27/2021
% Version: 1.0
%
% Description:
%   Run the type of trial for the specific amount of time. Changes screen
%   color and logs event times.
%
% Inputs:   src             [=] Used to check what type of trial to run
%           event           [=] Ignored. Not used as a callback function
%
% Outputs:  NONE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write to file button press, play sound, and prepare timer
switch src.String
    case 'Cued'
        trial_type = 1;
        trial_time = toc(src.Parent.UserData.clock);
        
        src.Parent.UserData.timer_cue.StartDelay = 3;
        src.Parent.UserData.timer_cue.TimerFcn = {@change_screen_color,trial_type,3};
        
    case 'Self-Initiated'
        trial_type = 2;
        trial_time = toc(src.Parent.UserData.clock);
        
    case 'Prepare-and-Go'
        trial_type = 3;
        trial_time = toc(src.Parent.UserData.clock);
        
        src.Parent.UserData.timer_prepare.StartDelay = 3;
        src.Parent.UserData.timer_prepare.TimerFcn = {@change_screen_color,trial_type,2};
        
        prep_time = round(src.Parent.UserData.random_timer_settings.std * randn + src.Parent.UserData.random_timer_settings.mean,3);
        src.Parent.UserData.timer_cue.StartDelay = 3+prep_time;
        src.Parent.UserData.timer_cue.TimerFcn = {@change_screen_color,trial_type,3};
end
play_sound(trial_type);
write_input(src.Parent,trial_type,0,trial_time);

% Clear screen of annotation and buttons
% Text and countdown to start of trial
for i = 5:-1:1
    clear_screen();
    
    str1 = sprintf('%s trial starts in:',src.String);
    dim1 = [0.125,0.70,0.75,0.1];
    str2 = num2str(i);
    dim2 = [0.125,0.40,0.75,0.1];
    annotation('textbox',dim1,'String',str1,...
        'FontName','Arial','FontSize',36,...
        'EdgeColor','none',...
        'HorizontalAlignment','center','VerticalAlignment','middle');
    annotation('textbox',dim2,'String',str2,...
        'FontName','Arial','FontSize',80,...
        'EdgeColor','none',...
        'HorizontalAlignment','center','VerticalAlignment','middle');
    
    pause(1);
end

% Change screen to red to indicate start of trial for cued and
% prepare-and-go
if trial_type == 1 || trial_type == 3
    change_screen_color(src,[],trial_type,1);
else
    change_screen_color(src,[],trial_type,3);
end

end