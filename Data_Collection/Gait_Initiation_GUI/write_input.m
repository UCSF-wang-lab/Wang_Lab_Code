function write_input(main_window,trial_type,trial_state,event_time)
%%
% write_input(...)
%
% Author: Kenneth Louie
% Project: General Gait Protocol
% Date: 05/27/2021
% Version: 1.0
%
% Description:
%   Helper function to write events and key presses to a binary data file
%   that is used as an auto save feature. In some cases garbage is placed
%   in some columns because it is an event marker.
%
% Inputs:   main_window     [=] Object handle for the figure window used
%                               for the task.
%           trial_type      [=] The type of trial that is current occuring
%                               (1 = cued; 2 = self-initiated; 3 =
%                               prepare-and-go)
%           trial_state     [=] State of the GI trial. (for cued: 1 = trial
%                               start, 3 = go cue, 4 = trial end; for
%                               self-initiated: 1 = start, 3 = go cue, 4 =
%                               trial end; for prepare-and-go: 1 = start, 2
%                               = prepare cue, 3 = go cue, 4 = trial end).
%                               0 = initial click of the trial type.
%           event_time      [=] Time since program started running when
%                               trial_state is updated 
%
% Outputs:  NONE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
fwrite(main_window.UserData.auto_save_file,[trial_type;trial_state;event_time],'double');
end