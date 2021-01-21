function varargout = processRCS(folder_path)
% Processes the raw JSON files from RC+S and converts them to a matlab
% struct.
% Author:   Kenneth Louie
% Date:     12/11/20

% Check if folder path is valid
if ~exist('folder_path','var') || isempty(folder_path)
    folder_path = uigetdir();
end

% these files should exist, if one is missing show error
% DeviceSettings must be processed first!
% (future: AdaptiveLog.json,DiagnosticLogs.json,ErrorLog.json,EventLog.json,StimLog.json,TimeSync.json)
files2Convert = {'DeviceSettings.json','RawDataTD.json',...
    'RawDataAccel.json','RawDataPower.json','RawDataFFT.json'};
filesNotFound = {};
for i = length(files2Convert)
    if ~isfile(fullfile(folder_path,files2Convert{i}))
        filesNotFound{end+1} = files2Convert{i};
    end
end

if ~isempty(filesNotFound)
    for i = 1:length(filesNotFound)
        if i == 1
            error_msg = 'Required files not found:\n';
        end
        error_msg = [error_msg,['\t',filesNotFound{i},'\n']];
    end
    error(sprintf(error_msg));
end

fprintf([repmat('=',1,100),'\n']);
fprintf('Processing directory: %s\n',folder_path)
for i = 1:length(files2Convert)
    filename = files2Convert{i};
    fprintf('Processing file: %s...\n',filename)
    switch filename
        case 'DeviceSettings.json'
            [timeDomainSettings,powerSettings,fftSettings,metaData] = processDeviceSettings(folder_path);
            DeviceSettings.timeDomainSettings = timeDomainSettings;
            DeviceSettings.powerSettings = powerSettings;
            DeviceSettings.fftSettings = fftSettings;
            DeviceSettings.metaData = metaData;
        case 'RawDataTD.json'
            timeDomainDataTable = processTDData(fullfile(folder_path,filename));
        case 'RawDataAccel.json'
            accelDataTable = processAccelData(fullfile(folder_path,filename));
        case 'RawDataPower.json'
            powerDataTable = processPowerData(folder_path,powerSettings);
        case 'RawDataFFT.json'
            fftDataTable = processFFTData(fullfile(folder_path,filename),fftSettings);
    end
    fprintf('Complete.\n\n');
end

% Save data
fprintf('Saving data...');
save(fullfile(folder_path,'DeviceSettings.mat'),'DeviceSettings');
save(fullfile(folder_path,'RawDataTD.mat'),'timeDomainDataTable');
save(fullfile(folder_path,'RawDataAccel.mat'),'accelDataTable');
save(fullfile(folder_path,'RawDataPower.mat'),'powerDataTable');
save(fullfile(folder_path,'RawDataFFT.mat'),'fftDataTable');
fprintf('Complete.\n\n');

fprintf('All files in directory processed.\n');
fprintf([repmat('=',1,100),'\n']);
end

function varargout = processDeviceSettings(folder_path)
[timeDomainSettings,powerSettings,fftSettings,metaData] = createDeviceSettingsTable(folder_path);
varargout{1} = timeDomainSettings;
varargout{2} = powerSettings;
varargout{3} = fftSettings;
varargout{4} = metaData;
end

function timeDomainData = processTDData(filename)
jsonobj_TD = deserializeJSON(filename);
if ~isempty(jsonobj_TD.TimeDomainData)
    fprintf('\tLoading Time Domain Data');
    [outtable_TD, srates_TD] = createTimeDomainTable(jsonobj_TD);
    fprintf('\tCreating derivedTimes for time domain:')
    timeDomainData = assignTime(outtable_TD);
else
    timeDomainData = [];
end
end

function accelData = processAccelData(filename)
jsonobj_Accel = deserializeJSON(filename);
if ~isempty(jsonobj_Accel.AccelData)
    fprintf('\tLoading Accelerometer Data');
    [outtable_Accel, srates_Accel] = createAccelTable(jsonobj_Accel);
    fprintf('\tCreating derivedTimes for accelerometer:');
    accelData = assignTime(outtable_Accel);
else
    accelData = [];
end
end

function powerData = processPowerData(folder_path,powerSettings)
fprintf('\tLoading Power Data');
% Checking if power data is empty happens within createPowerTable
% function
[outtable_Power] = createPowerTable(folder_path);

% Calculate power band cutoffs (in Hz) and add column to powerSettings
if ~isempty(outtable_Power)
    % Translate powerSettings.powerBands into Hz
    numSettings = size(powerSettings,1);
    for iSetting = 1:numSettings
        powerBands_toConvert = powerSettings.powerBands{iSetting};
        currentTDsampleRate = powerSettings.TDsampleRates(iSetting);
        currentFFTconfig = powerSettings.fftConfig(iSetting);
        [currentPowerBands] = getPowerBands(powerBands_toConvert,currentFFTconfig,currentTDsampleRate);
        powerSettings.powerBandsInHz(iSetting) = currentPowerBands;
    end
    
    % Add samplerate and packetsizes column to outtable_Power -- samplerate is inverse
    % of fftConfig.interval; in principle this interval could change
    % over the course of the recording
    
    % Determine if more than one sampling rate across recording
    for iSetting = 1:numSettings
        all_powerFs(iSetting) =  1/((powerSettings.fftConfig(iSetting).interval)/1000);
    end
    
    if length(unique(all_powerFs)) > 1
        error('More than one sampling rate for power channels -- code development needed')
        % Need to loop through each setting in powerSettings; find
        % closest times between powerSettings and outtable_Power to
        % assign corresponding samplerate
    else
        % Same sample rate for power data for the full file
        powerDomain_sampleRate = unique(all_powerFs);
        outtable_Power.samplerate(:) = powerDomain_sampleRate;
        outtable_Power.packetsizes(:) = 1;
    end
    
    powerData = assignTime(outtable_Power);
else
    powerData = [];
end
end

function FFTData = processFFTData(filename,fftSettings)
jsonobj_FFT = deserializeJSON(filename);
if ~isempty(jsonobj_FFT.FftData)
    fprintf('\tLoading FFT Data');
    outtable_FFT = createFFTtable(jsonobj_FFT);
    
    % Add FFT parameter info to fftSettings
    numSettings = size(fftSettings,1);
    for iSetting = 1:numSettings
        currentFFTconfig = fftSettings.fftConfig(iSetting);
        currentTDsampleRate = fftSettings.TDsampleRates(iSetting);
        fftParameters = getFFTparameters(currentFFTconfig,currentTDsampleRate);
        fftSettings.fftParameters(iSetting) = fftParameters;
    end
    % Add samplerate and packetsizes column to outtable_Power -- samplerate is inverse
    % of fftConfig.interval; in principle this interval could change
    % over the course of the recording
    
    % Determine if more than one sampling rate across recording
    for iSetting = 1:numSettings
        all_powerFs(iSetting) =  1/((fftSettings.fftConfig(iSetting).interval)/1000);
    end
    
    if length(unique(all_powerFs)) > 1
        error('More than one sampling rate for FFT channels -- code development needed')
        % Need to loop through each setting in fftSettings; find
        % closest times between fftSettings and outtable_FFT to
        % assign corresponding samplerate
    else
        % Same sample rate for FFT data for the full file
        FFT_sampleRate = unique(all_powerFs);
        outtable_FFT.samplerate(:) = FFT_sampleRate;
        outtable_FFT.packetsizes(:) = 1;
    end
    
    fprintf('\tCreating derivedTimes for FFT:');
    FFTData = assignTime(outtable_FFT);
else
    FFTData = [];
end
end