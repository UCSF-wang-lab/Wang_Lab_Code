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
    'RawDataAccel.json','RawDataPower.json','RawDataFFT.json',...
    'EventLog.json','StimLog.json','AdaptiveLog.json'};
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

shortGaps_systemTick = 0;

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
            timeDomainDataTable = processTDData(fullfile(folder_path,filename),shortGaps_systemTick);
        case 'RawDataAccel.json'
            accelDataTable = processAccelData(fullfile(folder_path,filename),shortGaps_systemTick);
        case 'RawDataPower.json'
            powerDataTable = processPowerData(folder_path,powerSettings,shortGaps_systemTick);
        case 'RawDataFFT.json'
            fftDataTable = processFFTData(fullfile(folder_path,filename),fftSettings,shortGaps_systemTick);
        case 'EventLog.json'
            eventLogDataTable = processEventLog(folder_path);
        case 'StimLog.json'
            stimLogDataTable = processStimLog(folder_path);
        case 'AdaptiveLog.json'
            adaptiveLogDataTable = processAdaptiveLog(fullfile(folder_path,filename),shortGaps_systemTick);
    end
    fprintf('Complete.\n\n');
end

% Save data
save_path = strrep(folder_path,'Raw Data','Processed Data');
slashes = strfind(save_path,'/');
if exist(save_path(1:slashes(end-1)),'dir') ~= 7
    mkdir(save_path(1:slashes(end-1)));
end

if exist(save_path(1:slashes(end)),'dir') ~= 7
    mkdir(save_path(1:slashes(end)));
end

if exist(save_path,'dir') ~= 7
    mkdir(save_path);
end

fprintf('Saving data...');
save(fullfile(save_path,'DeviceSettings.mat'),'DeviceSettings');
save(fullfile(save_path,'RawDataTD.mat'),'timeDomainDataTable');
save(fullfile(save_path,'RawDataAccel.mat'),'accelDataTable');
save(fullfile(save_path,'RawDataPower.mat'),'powerDataTable');
save(fullfile(save_path,'RawDataFFT.mat'),'fftDataTable');
save(fullfile(save_path,'EventLog.mat'),'eventLogDataTable');
save(fullfile(save_path,'StimLog.mat'),'stimLogDataTable');
save(fullfile(save_path,'AdaptiveLog.mat'),'adaptiveLogDataTable');
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

function timeDomainData = processTDData(filename,shortGaps_systemTick)
jsonobj_TD = deserializeJSON(filename);
if isfield(jsonobj_TD,'TimeDomainData') && ~isempty(jsonobj_TD.TimeDomainData)
    disp('Loading Time Domain Data')
    [outtable_TD, srates_TD] = createTimeDomainTable(jsonobj_TD);
    disp('Creating derivedTimes for time domain:')
    timeDomainData = assignTime(outtable_TD,shortGaps_systemTick);
else
    timeDomainData = [];
end
end

function accelData = processAccelData(filename,shortGaps_systemTick)
jsonobj_Accel = deserializeJSON(filename);
if isfield(jsonobj_Accel,'AccelData') && ~isempty(jsonobj_Accel.AccelData)
    disp('Loading Accelerometer Data')
    [outtable_Accel, srates_Accel] = createAccelTable(jsonobj_Accel);
    disp('Creating derivedTimes for accelerometer:')
    accelData = assignTime(outtable_Accel,shortGaps_systemTick);
else
    accelData = [];
end
end

function powerData = processPowerData(folder_path,powerSettings,shortGaps_systemTick)
disp('Loading Power Data')
% Checking if power data is empty happens within createPowerTable
% function
[outtable_Power] = createPowerTable(folder_path);

% Calculate power band cutoffs (in Hz) and add column to powerSettings
if ~isempty(outtable_Power)
    % Add samplerate and packetsizes column to outtable_Power -- samplerate is inverse
    % of fftConfig.interval
    numSettings = size(powerSettings,1);
    % Determine if more than one sampling rate across recording
    for iSetting = 1:numSettings
        all_powerFs(iSetting) =  1/((powerSettings.fftConfig(iSetting).interval)/1000);
    end
    
    if length(unique(all_powerFs)) > 1
        % Multiple sample rates for power data in the full file
        powerData = createDataTableWithMultipleSamplingRates(all_powerFs,powerSettings,outtable_Power,shortGaps_systemTick);
    else
        % Same sample rate for power data for the full file
        powerDomain_sampleRate = unique(all_powerFs);
        outtable_Power.samplerate(:) = powerDomain_sampleRate;
        outtable_Power.packetsizes(:) = 1;
        powerData = assignTime(outtable_Power,shortGaps_systemTick);
    end
end
end

function FFTData = processFFTData(filename,fftSettings,shortGaps_systemTick)
jsonobj_FFT = deserializeJSON(filename);
if isfield(jsonobj_FFT,'FftData') && ~isempty(jsonobj_FFT.FftData)
    disp('Loading FFT Data')
    outtable_FFT = createFFTtable(jsonobj_FFT);
    
    % Add FFT parameter info to fftSettings
    numSettings = size(fftSettings,1);
    for iSetting = 1:numSettings
        currentFFTconfig = fftSettings.fftConfig(iSetting);
        currentTDsampleRate = fftSettings.TDsampleRates(iSetting);
        fftParameters = getFFTparameters(currentFFTconfig,currentTDsampleRate);
        fftSettings.fftParameters(iSetting) = fftParameters;
    end
    % Add samplerate and packetsizes column to outtable_FFT -- samplerate is inverse
    % of fftConfig.interval; in principle this interval could change
    % over the course of the recording
    
    % Determine if more than one sampling rate across recording
    for iSetting = 1:numSettings
        all_fftFs(iSetting) =  1/((fftSettings.fftConfig(iSetting).interval)/1000);
    end
    
    if length(unique(all_fftFs)) > 1
        FFTData = createDataTableWithMultipleSamplingRates(all_fftFs,fftSettings,outtable_FFT,shortGaps_systemTick);
    else
        % Same sample rate for FFT data for the full file
        FFT_sampleRate = unique(all_powerFs);
        outtable_FFT.samplerate(:) = FFT_sampleRate;
        outtable_FFT.packetsizes(:) = 1;
        disp('Creating derivedTimes for FFT:')
        FFTData = assignTime(outtable_FFT,shortGaps_systemTick);
    end
else
    FFTData = [];
end
end

function eventLogData = processEventLog(folderPath)
disp('Collecting Event Information from Event Log file')
eventLogData = createEventLogTable(folderPath);
end

function stimLogData = processStimLog(folderPath)
disp('Collecting Stimulation Settings from Device Settings file')
DeviceSettings_fileToLoad = [folderPath filesep 'DeviceSettings.json'];
if isfile(DeviceSettings_fileToLoad)
    [stimSettingsOut, stimMetaData] = createStimSettingsFromDeviceSettings(folderPath);
else
    warning('No DeviceSettings.json file - could not extract stimulation settings')
end

disp('Collecting Stimulation Settings from Stim Log file')
StimLog_fileToLoad = [folderPath filesep 'StimLog.json'];
if isfile(StimLog_fileToLoad)
    [stimLogData] = createStimSettingsTable(folderPath,stimMetaData);
else
    warning('No StimLog.json file')
end
end

function AdaptiveData = processAdaptiveLog(filename,shortGaps_systemTick)
jsonobj_Adaptive = deserializeJSON(filename);
if isfield(jsonobj_Adaptive,'AdaptiveUpdate') && ~isempty(jsonobj_Adaptive(1).AdaptiveUpdate)
    try
        disp('Loading Adaptive Data')
        outtable_Adaptive = createAdaptiveTable(jsonobj_Adaptive);
        % Note: StateTime must still be converted to sec in
        % outtable_Adaptive
        
        % Calculate adaptive_sampleRate - determine if more than one
        if size(fftSettings,1) == 1
            adaptive_sampleRate =  1/((fftSettings.fftConfig(1).interval)/1000);
            outtable_Adaptive.samplerate(:) = adaptive_sampleRate;
            outtable_Adaptive.packetsizes(:) = 1;
            outtable_Adaptive.StateTime = outtable_Adaptive.StateTime * (fftSettings.fftConfig(1).interval/1000);
            
            disp('Creating derivedTimes for Adaptive:')
            AdaptiveData = assignTime(outtable_Adaptive, shortGaps_systemTick);
        else
            for iSetting = 1:size(fftSettings,1)
                all_adaptiveFs(iSetting) =  1/((fftSettings.fftConfig(iSetting).interval)/1000);
            end
            if length(unique(all_adaptiveFs)) > 1
                AdaptiveData = createDataTableWithMultipleSamplingRates(all_adaptiveFs,fftSettings,outtable_Adaptive,shortGaps_systemTick);
            else
                adaptive_sampleRate = all_adaptiveFs(1);
                outtable_Adaptive.samplerate(:) = adaptive_sampleRate;
                outtable_Adaptive.packetsizes(:) = 1;
                
                disp('Creating derivedTimes for Adaptive:')
                AdaptiveData = assignTime(outtable_Adaptive, shortGaps_systemTick);
            end
        end
    catch
        warning('Adaptive data present but error when trying to process. Will be excluded.')
        AdaptiveData = [];
    end
else
    AdaptiveData = [];
end
end