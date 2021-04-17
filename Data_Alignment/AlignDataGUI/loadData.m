function loadData(src,varargin)

if contains(src.Tag,'Xsens') || contains(src.Tag,'Force')
    [filename,path] = uigetfile('.csv');
else
    [filename,path] = uigetfile();
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
end

src.Parent.Parent.UserData.file_names{end+1} = fullfile(path,filename);

% Add event to logger
addEvent(str);

% Update data selection options for plot windows
updatePlotSelectionOptions([],[]);

end