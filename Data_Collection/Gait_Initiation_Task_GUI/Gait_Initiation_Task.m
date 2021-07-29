function Gait_Initiation_Task()
%% Header
% Author:   Kenneth Louie
% Date:     05/20/2021
% Project:  General Gait Protocol
% Version:  1.0
% 
% Description:
%   Main file to run the Gait Initiation task of the general gait protocol.
%   Reads in configuration settings and changes the behavior of the program
%   appropriately. 
%
% Inputs:   None
% Output:   None

clc; clear all;
delete(timerfindall)

%% Set processor priority
if ispc
    % "High"
    [~,~] = dos('wmic process where name="matlab.exe" CALL setpriority 128');
else
end

%% Load config file
config_settings = load_config();

%% Initialize GUI
% Create fullscreen figure
main_window = figure('Visible','on','units','normalized','outerposition',[0,0,1,1],...
    'MenuBar','none',...
    'DockControls','off',...
    'NumberTitle','off',...
    'Name','Gait Initiation',...
    'color',[128,128,128]./256,...
    'WindowState','fullscreen',...
    'UserData',struct(),...
    'Tag','GI_Task_main_window');

% Create struct of user data. This is so we can pass around data between
% the GUI.
main_window.UserData.subid = config_settings.subid;
main_window.UserData.save_path = config_settings.save_path;
main_window.UserData.colors = struct('go',config_settings.color_go,'prepare',config_settings.color_prepare,'stop',config_settings.color_stop);
main_window.UserData.random_timer_settings = config_settings.random_time;
main_window.UserData.clock = tic;
main_window.UserData.ready = 0;
main_window.UserData.timer_cue = timer('ExecutionMode','singleShot','TimerFcn',@change_screen_color);
main_window.UserData.timer_prepare = timer('ExecutionMode','singleShot','TimerFcn',@change_screen_color);
main_window.UserData.timer_stop = timer('ExecutionMode','singleShot','TimerFcn',@change_screen_color);
main_window.UserData.n_pg_trials = 0;


%% Call play sound to initialize the local variables
play_sound(0);

%% Start program
trial_option_screen(main_window);
end