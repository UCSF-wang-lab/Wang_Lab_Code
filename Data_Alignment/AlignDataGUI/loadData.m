function loadData(src,varargin)

if isempty(src.Parent.Parent.UserData.basePath)
    basePath = pwd;
else
    basePath = src.Parent.Parent.UserData.basePath;
end


if contains(src.Tag,'Xsens') || contains(src.Tag,'Force') || contains(src.Tag,'rover')
    [filename,path] = uigetfile([basePath,'/*.csv']);
elseif contains(src.Tag,'Delsys')
    [filename,path] = uigetfile([basePath,'/*.mat']);
else
    path = uigetdir();
    if path == 0
        addEvent('No folder selected.');
        return;
    end
end

if exist('filename','var') && exist('path','var')
    try
        [~,~,ext] = fileparts(filename);
    catch
        addEvent('No file selected.');
        return;
    end
    
    try
        switch ext
            case '.mat'
                data = load(fullfile(path,filename));
            case '.csv'
                data = readtable(fullfile(path,filename));
        end
    catch
        warning('Could not open selected file.');
    end
    
elseif ~exist('filename','var') && exist('path','var')
    file_list = dir(fullfile(path,'*.mat'));

    % Check if there is a filtered time domain file
    if sum(cellfun(@(x)contains(x,'filtered'),{file_list.name}))>0
        filt_data_array = cellfun(@(x)contains(x,'filtered'),{file_list.name});
        td_data_array = cellfun(@(x)contains(x,'RawDataTD'),{file_list.name});
        remove_inds = xor(filt_data_array,td_data_array);
        file_list(remove_inds) = [];
    end

    filenames = [];
    for i = 1:length(file_list)
        if file_list(i).isdir == 0
            filename = fullfile(file_list(i).folder,file_list(i).name);
            filenames(end+1) = filename;
            switch file_list(i).name
                case "DeviceSettings.mat"
                    settings_data = load(filename);
                    str = sprintf('Loaded RC+S DeviceSettings file: %s',filename);
                    addEvent(str);
                case "LogTable.mat"
                    log_data = load(filename);
                    str = sprintf('Loaded RC+S LogTable file: %s',filename);
                    addEvent(str);
                case "RawDataAccel.mat"
                    accel_data = load(filename);
                    str = sprintf('Loaded RC+S Acceleration file: %s',filename);
                    addEvent(str);
                case {"RawDataTD.mat","RawDataTD_filtered.mat"}
                    td_data = load(filename);
                    str = sprintf('Loaded RC+S LFP time domain file: %s',filename);
                    addEvent(str);
            end
        end
    end
    filename = filenames;
end
    

switch src.Tag
    case 'Left_INS_button'
        if exist('td_data','var')
            src.Parent.Parent.UserData.LFP_data.Left = td_data;
            src.Parent.Parent.UserData.indicators.left_LFP.ForegroundColor = [20,136,7]./256;
        end
        
        if exist('accel_data','var')
            src.Parent.Parent.UserData.Accel_data.Left = accel_data;
            src.Parent.Parent.UserData.indicators.left_accel.ForegroundColor = [20,136,7]./256;
        end
        
        if exist('settings_data','var')
            src.Parent.Parent.UserData.DeviceSettings.Left = settings_data.DeviceSettings;
            src.Parent.Parent.UserData.indicators.left_settings.ForegroundColor = [20,136,7]./256;
        end
        
        if exist('log_data','var')
            src.Parent.Parent.UserData.LogTable.Left = log_data;
            src.Parent.Parent.UserData.indicators.left_logs.ForegroundColor = [20,136,7]./256;
        end
        src.Parent.Parent.UserData.indicators.left_INS.Value = true;
    case 'Right_INS_button'
        if exist('td_data','var')
            src.Parent.Parent.UserData.LFP_data.Right = td_data;
            src.Parent.Parent.UserData.indicators.right_LFP.ForegroundColor = [20,136,7]./256;
        end
        
        if exist('accel_data','var')
            src.Parent.Parent.UserData.Accel_data.Right = accel_data;
            src.Parent.Parent.UserData.indicators.right_accel.ForegroundColor = [20,136,7]./256;
        end
        
        if exist('settings_data','var')
            src.Parent.Parent.UserData.DeviceSettings.Right = settings_data.DeviceSettings;
            src.Parent.Parent.UserData.indicators.right_settings.ForegroundColor = [20,136,7]./256;
        end
        
        if exist('log_data','var')
            src.Parent.Parent.UserData.LogTable.Right = log_data;
            src.Parent.Parent.UserData.indicators.right_logs.ForegroundColor = [20,136,7]./256;
        end
        
        src.Parent.Parent.UserData.indicators.right_INS.Value = true;
    case 'Delsys_button'
        src.Parent.Parent.UserData.Delsys_data = data;
        src.Parent.Parent.UserData.indicators.delsys.Value = true;
        str = sprintf('Loaded Delsys file: %s',fullfile(path,filename));
    case 'Xsens_button'
        src.Parent.Parent.UserData.Xsens_data = data;
        src.Parent.Parent.UserData.indicators.xsens.Value = true;
        str = sprintf('Loaded Xsens file: %s',fullfile(path,filename));
    case 'Force_plate_button'
        src.Parent.Parent.UserData.FP_data = data;
        src.Parent.Parent.UserData.indicators.fp.Value = true;
        str = sprintf('Loaded gait initiation force plate file: %s',fullfile(path,filename));
    case 'Left_rover_button'
        src.Parent.Parent.UserData.Rover_data.Left = data;
        src.Parent.Parent.UserData.indicators.left_rover.Value = true;
        str = sprintf('Loaded left Rover file: %s',fullfile(path,filename));
    case 'Right_rover_button'
        src.Parent.Parent.UserData.Rover_data.Right = data;
        src.Parent.Parent.UserData.indicators.right_rover.Value = true;
        str = sprintf('Loaded right Rover file: %s',fullfile(path,filename));
end

for i = 1:length(filename)
    src.Parent.Parent.UserData.file_names{end+1} = fullfile(path,filename(i));
end

% if isempty(src.Parent.Parent.UserData.basePath)
%     curr_path = path(1:end-1);
%     [parent_dir,curr_dir,~] = fileparts(curr_path);
%     while ~strcmp(curr_dir,'Processed Data')
%          curr_path = parent_dir;
%          [parent_dir,curr_dir,~] = fileparts(curr_path);
%     end
%     src.Parent.Parent.UserData.basePath = curr_path;
% end

% Add event to logger
if ~contains(src.Tag,'INS')
    addEvent(str);
end

% Update data selection options for plot windows
updatePlotSelectionOptions([],[]);

end