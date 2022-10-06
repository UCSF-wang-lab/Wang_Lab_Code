function loadData(src,varargin)

if isempty(src.Parent.Parent.UserData.basePath)
    basePath = pwd;
else
    basePath = src.Parent.Parent.UserData.basePath;
end


if contains(src.Tag,'Xsens') || contains(src.Tag,'Force') || contains(src.Tag,'Teensey')
    [filename,path] = uigetfile([basePath,'/*.csv']);
else
    [filename,path] = uigetfile(basePath);
end

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

switch src.Tag
    case 'Left_LFP_button'
        src.Parent.Parent.UserData.LFP_data.Left = data;
        src.Parent.Parent.UserData.indicators.left_LFP.Value = true;
        str = sprintf('Loaded RC+S left brain LFP time domain file: %s',fullfile(path,filename));
    case 'Right_LFP_button'
        src.Parent.Parent.UserData.LFP_data.Right = data;
        src.Parent.Parent.UserData.indicators.right_LFP.Value = true;
        str = sprintf('Loaded RC+S right brain LFP time domain file: %s',fullfile(path,filename));
    case 'Left_Accel_button'
        src.Parent.Parent.UserData.Accel_data.Left = data;
        src.Parent.Parent.UserData.indicators.left_accel.Value = true;
        str = sprintf('Loaded RC+S left brain acceleration file: %s',fullfile(path,filename));
    case 'Right_Accel_button'
        src.Parent.Parent.UserData.Accel_data.Right = data;
        src.Parent.Parent.UserData.indicators.right_accel.Value = true;
        str = sprintf('Loaded RC+S right brain acceleration file: %s',fullfile(path,filename));
    case 'Left_Device_Settings_button'
        src.Parent.Parent.UserData.DeviceSettings.Left = data.DeviceSettings;
        src.Parent.Parent.UserData.indicators.left_device_settings.Value = true;
        str = sprintf('Loaded RC+S left brain recording setting file: %s',fullfile(path,filename));
    case 'Right_Device_Settings_button'
        src.Parent.Parent.UserData.DeviceSettings.Right = data.DeviceSettings;
        src.Parent.Parent.UserData.indicators.right_device_settings.Value = true;
        str = sprintf('Loaded RC+S right brain recording setting file: %s',fullfile(path,filename));
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
    case 'Teensey_button'
        src.Parent.Parent.UserData.Teensey_data = data;
        src.Parent.Parent.UserData.indicators.teensey.Value = true;
        str = sprintf('Loaded Teensey file: %s',fullfile(path,filename));
end

src.Parent.Parent.UserData.file_names{end+1} = fullfile(path,filename);

if isempty(src.Parent.Parent.UserData.basePath)
    curr_path = path(1:end-1);
    [parent_dir,curr_dir,~] = fileparts(curr_path);
    while ~strcmp(curr_dir,'Processed Data')
         curr_path = parent_dir;
         [parent_dir,curr_dir,~] = fileparts(curr_path);
    end
    src.Parent.Parent.UserData.basePath = curr_path;
end

% Add event to logger
addEvent(str);

% Update data selection options for plot windows
updatePlotSelectionOptions([],[]);

end