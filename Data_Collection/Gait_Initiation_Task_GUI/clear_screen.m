function clear_screen(src,event)
%%
% clear_screen(...)
%
% Author: Kenneth Louie
% Project: General Gait Protocol
% Date: 05/27/2021
% Version: 1.0
%
% Description:
%   Helper function to clear the screen of the task. 
%
% Inputs:   src             [=] Ignored. Not used as a callback function
%           event           [=] Ignored. Not used as a callback function
%
% Outputs:  NONE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
main_window = findobj('Tag','GI_Task_main_window');
delete(findall(main_window,'type','annotation'));

set(main_window.UserData.trial_buttons,'Visible','Off');
end