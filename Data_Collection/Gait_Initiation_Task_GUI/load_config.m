function config_settings = load_config()
%%
% load_config(...)
%
% Author: Kenneth Louie
% Project: General Gait Protocol
% Date: 05/20/2021
% Version: 1.0
%
% Description:
%   Reads in a configuration text file that will set various paramters of
%   the task.
%
% Inputs:   NONE
%
% Outputs:  NONE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Default settings
config_settings.subid = [];
config_settings.save_path = [];
config_settings.color_go = [];
config_settings.color_prepare = [];
config_settings.color_stop = [];
config_settings.random_time_mean = [];
config_settings.random_time_std = [];

% Try to open config file
try
    fid = fopen('test_config.txt');
catch
    error('File cannot be opened.');
end

% Read content of config file
val = fgetl(fid);
while val ~= -1
    if contains(val,'ID')
        config_settings.subid = readString(val);
    elseif contains(val,'save_path')
        config_settings.save_path = readString(val);
    elseif contains(val,'color_go')
        config_settings.color_go = readString(val);
    elseif contains(val,'color_prepare')
        config_settings.color_prepare = readString(val);
    elseif contains(val,'color_stop')
        config_settings.color_stop = readString(val);
    elseif contains(val,'random_time_mean')
        config_settings.random_time_mean = getVal(val);
    elseif contains(val,'random_time_std')
        config_settings.random_time_std = getVal(val);
    end
    val = fgetl(fid);
end

fclose(fid);
end

function name = readString(val)
if contains(val,'"')
    loc = strfind(val,'"');
    name = val(loc(1)+1:loc(2)-1);
else
    name = [];
end
end

function value = getVal(val)
if contains(val,'=')
    loc = strfind(val,'=');
    value = str2num(val(loc+2:end));
else
    value = [];
end
end