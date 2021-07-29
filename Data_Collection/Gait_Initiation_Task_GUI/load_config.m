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
config_settings.random_time = [];

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
    elseif contains(val,'random_time')
        config_settings.random_time = getVal(val);    
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
value = [];
if contains(val,'=')
    loc = strfind(val,'=');
    if contains(val,',')
        loc_comma = strfind(val,',');
        locs = [loc,loc_comma];
        for i = 1:length(locs)-1
            value = [value,str2double(val(locs(i)+1:locs(i+1)-1))];
        end
        value = [value,str2double(val(locs(end):end))];
    else
        value = str2num(val(loc+2:end));
    end
elseif contains(val,',')
        loc_comma = strfind(val,',');
else
    value = [];
end
end