function artifactRemovalRCSLFP(lfp_data,filt_type,varargin)
% stim_freq,sampling_rate,ekg_template_indexes,save_data
% lfp_data can either be a file path to the raw time domain LFP data or a
% table object
%
% filt_type can be "stim" or "ekg" or "stim_ekg"
%
% stim_freq is a float value
%
% sampling_rate is a float value
%
% ekg_template is a 4x2 matrix, where each row is the template EKG to
% remove
%
% Example run:
%   file_path = '/Volumes/dwang3_shared/Patient Data/RC+S Data/gait_RCS_02/v7_2022-10-31_pre-programming_second_side_implant/Data/Processed Data/RC+S/Left INS/Session1667247503283/DeviceNPC700517H/RawDataTD.mat';
%   ekg_template_indexes(1,:) = [39992,40311]   % Index of the ekg artifact template
%   artifactRemovalRCSLFP(file_path,'stim_gait_ekg','sampling_rate',500,'stim_freq',150.6,'ekg_template_indexes',ekg_template_indexes,'save_data',false,'show_plots',true)


for i = 1:2:nargin-2
    switch varargin{i}
        case 'stim_freq'
            stim_freq = varargin{i+1};      % Stimulation frequency
        case 'sampling_rate'
            sampling_rate = varargin{i+1};  % Sampling rate of the LFP data
        case 'ekg_template_indexes'
            ekg_template_indexes = varargin{i+1};   % 4x2 matrix
        case 'show_plots'
            show_plots = varargin{i+1}; % Boolean, true or false
        case 'save_data'
            save_data = varargin{i+1};  % Boolean, true or false
    end
end

% Set default values
if isempty(filt_type)
    filt_type = 'none';
end

if ~exist('stim_freq','var') || isempty(stim_freq)
    stim_freq = 130.2;
end

if ~exist('sampling_rate','var') || isempty(sampling_rate)
    sampling_rate = 500;
end

if ~exist('ekg_template_indexes','var') || isempty(ekg_template_indexes)
    ekg_template_indexes = nan(4,2);
end

if ~exist('show_plots','var') || isempty(show_plots)
    show_plots = false;
end

if ~exist('save_data','var') || isempty(save_data)
    save_data = false;
end

% check to see if the lfp data passed in was a file path or the table
if ~istable(lfp_data)
    file_path = lfp_data;
    lfp_data = load(file_path);
    lfp_data = lfp_data.timeDomainDataTable;
end

% Grab the recording channels indicated by the keyword "key"
table_variables = lfp_data.Properties.VariableNames;
recording_channels = table_variables(contains(table_variables,'key'));

% Setup output variable
filt_data = lfp_data;

% Show unfiltered data
if show_plots
    figure;
        for i = 1:length(recording_channels)
            subplot(2,2,i);
            plot(lfp_data.(recording_channels{i}));
            ylabel('mV');
            xlabel('Data Sample')
            title(recording_channels{i});
        end
        sgtitle('Unfiltered LFP');
        
        figure;
        for i = 1:length(recording_channels)
            subplot(2,2,i);
            pwelch(lfp_data.(recording_channels{i}),sampling_rate,round(sampling_rate*0.9),[1:200],sampling_rate);
            ylabel('dB/Hz');
            xlabel('Frequency')
            title(recording_channels{i});
        end
        sgtitle('Unfiltered LFP PSD');
end

%% Remove walking artifacts (~1-2Hz) and any high frequency noise (+180Hz)
if contains(filt_type,'gait')
    [b,a] = butter(4,180/(sampling_rate/2),'low');
    [d,c] = butter(6,2/(sampling_rate/2),'high');
    
    for i = 1:length(recording_channels)
        filt_data.(recording_channels{i}) = filtfilt(b,a,filt_data.(recording_channels{i}));
        filt_data.(recording_channels{i}) = filtfilt(d,c,filt_data.(recording_channels{i}));
    end
    
    % Show the filtered LFP data and power spectral density
    if show_plots
        figure;
        for i = 1:length(recording_channels)
            subplot(2,2,i);
            plot(filt_data.(recording_channels{i}));
            ylabel('mV');
            xlabel('Data Sample')
            title(recording_channels{i});
        end
        sgtitle('Gait Filtered LFP');
    
        figure;
        for i = 1:length(recording_channels)
            subplot(2,2,i);
            pwelch(filt_data.(recording_channels{i}),sampling_rate,round(sampling_rate*0.9),[1:200],sampling_rate);
            ylabel('dB/Hz');
            xlabel('Frequency')
            title(recording_channels{i});
        end
        sgtitle('Gait Filtered LFP PSD');
    end
end

%% Filter stimulation artifacts
if contains(filt_type,'stim')
    [b,a] = butter(6,[stim_freq-7.5,stim_freq+7.5]/(sampling_rate/2),'stop');
    
    for i = 1:length(recording_channels)
        filt_data.(recording_channels{i}) = filtfilt(b,a,filt_data.(recording_channels{i}));
    end
    
    % Show the raw LFP data and power spectral density
    if show_plots
        % Show the filtered LFP data and power spectral density
        figure;
        for i = 1:length(recording_channels)
            subplot(2,2,i);
            plot(filt_data.(recording_channels{i}));
            ylabel('mV');
            xlabel('Data Sample')
            title(recording_channels{i});
        end
        if contains(filt_type,'gait')
            sgtitle('Gait and Stim Filtered LFP');
        else
            sgtitle('Stim Filtered LFP');
        end
        
        figure;
        for i = 1:length(recording_channels)
            subplot(2,2,i);
            pwelch(filt_data.(recording_channels{i}),sampling_rate,round(sampling_rate*0.9),[1:200],sampling_rate);
            ylabel('dB/Hz');
            xlabel('Frequency')
            title(recording_channels{i});
        end
        if contains(filt_type,'gait')
            sgtitle('Gait and Stim Filtered LFP');
        else
            sgtitle('Stim Filtered LFP');
        end
    end
end

%% Remove EKG artifacts
if contains(filt_type,'ekg')
    for i = 1:size(ekg_template_indexes,1)
        if ~isnan(ekg_template_indexes(i,1))
           ekg_template = filt_data.(recording_channels{i})(ekg_template_indexes(i,1):ekg_template_indexes(i,2));
           [sigClean, artRem, artLog, template] = RemoveNoiseTempMatch(filt_data.(recording_channels{i}), sampling_rate, ekg_template, [1 5], [0,0], [], 0, 1, 1);
           filt_data.(recording_channels{i}) = sigClean;
        end
    end
    
    % Show the EKG and stim filtered LFP data and power spectral density
    if show_plots
        figure;
        for i = 1:length(recording_channels)
            subplot(2,2,i);
            plot(filt_data.(recording_channels{i}));
            ylabel('mV');
            xlabel('Data Sample')
            title(recording_channels{i});
        end
        if contains(filt_type,'stim') && contains(filt_type,'gait')
            sgtitle('Gait, Stim, and EKG Filtered LFP');
        elseif contains(filt_type,'stim') && ~contains(filt_type,'gait')
            sgtitle('Stim and EKG Filtered LFP');
        else
            sgtitle('EKG Filtered LFP')
        end
        
        figure;
        for i = 1:length(recording_channels)
            subplot(2,2,i);
            pwelch(filt_data.(recording_channels{i}),sampling_rate,round(sampling_rate*0.9),[1:200],sampling_rate);
            ylabel('dB/Hz');
            xlabel('Frequency')
            title(recording_channels{i});
        end
        if contains(filt_type,'stim') && contains(filt_type,'gait')
            sgtitle('Gait, Stim, and EKG Filtered LFP');
        elseif contains(filt_type,'stim') && ~contains(filt_type,'gait')
            sgtitle('Stim and EKG Filtered LFP');
        else
            sgtitle('EKG Filtered LFP')
        end
    end
end

%% Save data if the file path was originally passed in
if exist('file_path','var') && save_data
    file_path2 = strrep(file_path,'TD.mat','TD_original.mat');
    copyfile(file_path,file_path2);
    
    % Copy filt data as timeDomainDataTable to replicate original file
    timeDomainDataTable = filt_data;
    [A,B,C] = fileparts(file_path);
    save_name = fullfile(A,[B,'_filtered',C]);
    save(save_name,'timeDomainDataTable');
end
end