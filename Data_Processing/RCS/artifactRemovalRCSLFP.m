function filt_data = artifactRemovalRCSLFP(input,filt_type,varargin)
% stim_freq,sampling_rate,ekg_template_indexes,save_data
% lfp_data can either be a file path to the raw time domain LFP data or a
% table object
%
% filt_type can be "stimTherapy", "stimSync", "ekg", "gait" and any
% combinations
%
% stim_freq is a float value
%
% sampling_rate is a float value
%
% ekg_template is a 2x1 cell array that contains a 4x2 matrix, where each 
% row is the template EKG to remove for that channel
%
% Example run:
%   Aligned Data with left and right side
%   file_path = '/Volumes/dwang3_shared/Patient Data/RC+S Data/gait_RCS_02/v8_2023-02-28_2nd-side_DbsOpt/Data/Aligned Data/gait_RCS_02_OG_ON_Stim_ON_Med_Trial_1_w_Gait_Events_Julia.mat';
%   filt_type = {'stim_gait_ekg','stim_gait_ekg','stim_gait','stim_gait'} % Left side implant, each column is a key
%   filt_type(2,:) = {'stim_gait','stim_gait_ekg','stim_gait','stim_gait'} % Right side implant, each colum is a key
%   ekg_template_indexes{1,1} = [38866,39177;34147,34461]
%   ekg_template_indexes{2,1} = [nan,nan;42484,42586]
%   artifactRemovalRCSLFP(file_path,filt_type,'sampling_rate',500,'ekg_template_indexes',ekg_template_indexes,'stim_freq',150.6,'save_data',false,'show_plots',true)

%% Parse inputs
for i = 1:2:nargin-2
    switch varargin{i}
        case 'stim_freq'
            stim_freq = varargin{i+1};      % Stimulation frequency
        case 'sampling_rate'
            sampling_rate = varargin{i+1};  % Sampling rate of the LFP data
        case 'wide_band_stim_freq_keys'
            wide_band_stim_freq_keys = varargin{i+1};   % Which recording key contains the stimulation contact. Input is a value between 0-3
        case 'ekg_template_indexes'
            ekg_template_indexes = varargin{i+1};   % 4x2 matrix
        case 'remove_stim_harmonic'
            remove_stim_harm = varargin{i+1};  % Boolean, true or false
        case 'show_plots'
            show_plots = varargin{i+1}; % Boolean, true or false
        case 'save_data'
            save_data = varargin{i+1};  % Boolean, true or false
    end
end

%% Set default values
if isempty(filt_type)
    filt_type = {'none'};   % cell array for the filt type for each channel
end

if ischar(filt_type)
    filt_type = {filt_type};
end

if ~exist('stim_freq','var') || isempty(stim_freq)
    stim_freq = 130.2;
end

if ~exist('sampling_rate','var') || isempty(sampling_rate)
    sampling_rate = 500;
end

if ~exist("wide_band_stim_freq_keys",'var')
    wide_band_stim_freq_keys = [];
end

if ~exist('ekg_template_indexes','var') || isempty(ekg_template_indexes)
    ekg_template_indexes = {nan(4,2)};
end

if ~iscell(ekg_template_indexes)
    ekg_template_indexes = {ekg_template_indexes};
end

if ~exist('remove_stim_harm','var') || isempty(remove_stim_harm)
    remove_stim_harm = false;
end

if ~exist('show_plots','var') || isempty(show_plots)
    show_plots = false;
end

if ~exist('save_data','var') || isempty(save_data)
    save_data = false;
end

%% check to see if the lfp data passed in was a file path or the table
if ~istable(input)
    file_path = input;
    data = load(file_path);

    if isfield(data,'aligned_data')
        % aligned data was passed in
        if isfield(data.aligned_data,'left_LFP_table') || isfield(data.aligned_data,'right_LFP_table')
            % aligned data was passed in
            count = 1;
            if isfield(data.aligned_data,'left_LFP_table')
                lfp_data{count,1} = data.aligned_data.left_LFP_table;
                count = count + 1;
            end

            if isfield(data.aligned_data,'right_LFP_table')
                lfp_data{count,1} = data.aligned_data.right_LFP_table;
            end
        end
    else
        lfp_data{1} = data.timeDomainDataTable;    % The json file path was passed in.
    end
else
    lfp_data{1} = input;
end

%% Check to see that filt_type and ekg_template_index are the correct size
if size(filt_type,1) ~= length(lfp_data)
    filt_type = repmat(filt_type,length(lfp_data),1);
end

if size(filt_type,2) ~= 4
    n_col = size(filt_type,2);
    for i = 1:size(filt_type,1)
        for j = n_col+1:4
            filt_type{i,j} = 'none';
        end
    end
end

if size(ekg_template_indexes,1) ~= length(lfp_data)
    for i = size(ekg_template_indexes,1)+1:length(lfp_data)
        ekg_template_indexes{i,1} = nan(4,2);
    end
end

%% Filter data
for j = 1:length(lfp_data)
    % Filter flags
    stim_sync_filtered = false;
    gait_hf_filtered = false;
    stim_filtered = false;
    ekg_filtered = false;

    % Grab the recording channels indicated by the keyword "key"
    table_variables = lfp_data{j}.Properties.VariableNames;
    recording_channels = table_variables(contains(table_variables,'key'));

    % Setup output variable
    filt_data{j} = lfp_data{j};

    %% Show unfiltered data
    if show_plots
        figure;
        for i = 1:length(recording_channels)
            subplot(2,2,i);
            plot(lfp_data{j}.(recording_channels{i}));
            ylabel('mV');
            xlabel('Data Sample')
            title(recording_channels{i});
        end
        sgtitle('Unfiltered LFP');

        figure;
        for i = 1:length(recording_channels)
            subplot(2,2,i);
            pwelch(lfp_data{j}.(recording_channels{i}),sampling_rate,round(sampling_rate*0.9),[1:200],sampling_rate);
            ylabel('dB/Hz');
            xlabel('Frequency')
            title(recording_channels{i});
        end
        sgtitle('Unfiltered LFP PSD');
    end

    %% Remove walking artifacts (~1-2Hz) and any high frequency noise (+180Hz)
    [b,a] = butter(4,180/(sampling_rate/2),'low');
    [d,c] = butter(6,2/(sampling_rate/2),'high');

    for i = 1:length(recording_channels)
        if contains(filt_type{j,i},'gait')
            filt_data{j}.(recording_channels{i}) = filtfilt(b,a,filt_data{j}.(recording_channels{i}));
            filt_data{j}.(recording_channels{i}) = filtfilt(d,c,filt_data{j}.(recording_channels{i}));
            gait_hf_filtered = true;
        end
    end

    % Show the filtered LFP data and power spectral density
    if show_plots && gait_hf_filtered
        figure;
        for i = 1:length(recording_channels)
            subplot(2,2,i);
            plot(filt_data{j}.(recording_channels{i}));
            ylabel('mV');
            xlabel('Data Sample')
            title(recording_channels{i});
        end
        sgtitle('Gait Filtered LFP');

        figure;
        for i = 1:length(recording_channels)
            subplot(2,2,i);
            pwelch(filt_data{j}.(recording_channels{i}),sampling_rate,round(sampling_rate*0.9),[1:200],sampling_rate);
            ylabel('dB/Hz');
            xlabel('Frequency')
            title(recording_channels{i});
        end
        sgtitle('Gait Filtered LFP PSD');
    end

    %% Filter stim sync pulses (5Hz)
    [b,a] = butter(4,5/(sampling_rate/2),'high');

    for i = 1:length(recording_channels)
        if contains(filt_type{j,i},'stimSync')
            filt_data{j}.(recording_channels{i}) = filtfilt(b,a,filt_data{j}.(recording_channels{i}));
            stim_sync_filtered = true;
        end
    end

    % Show the filtered LFP data and power spectral density
    if show_plots && stim_sync_filtered
        figure;
        for i = 1:length(recording_channels)
            subplot(2,2,i);
            plot(filt_data{j}.(recording_channels{i}));
            ylabel('mV');
            xlabel('Data Sample')
            title(recording_channels{i});
        end
        if gait_hf_filtered
            sgtitle('Gait and Stim Sync Filtered LFP');
        else
            sgtitle('Stim Sync Filtered LFP');
        end
        

        figure;
        for i = 1:length(recording_channels)
            subplot(2,2,i);
            pwelch(filt_data{j}.(recording_channels{i}),sampling_rate,round(sampling_rate*0.9),[1:200],sampling_rate);
            ylabel('dB/Hz');
            xlabel('Frequency')
            title(recording_channels{i});
        end
        if gait_hf_filtered
            sgtitle('Gait and Stim Sync Filtered LFP');
        else
            sgtitle('Stim Sync Filtered LFP');
        end
    end

    %% Filter stimulation artifacts
    stim_filt_done = false;
    count = 0;

    while ~stim_filt_done
        [b,a] = butter(6,[(stim_freq/(2^count))-7.5,(stim_freq/(2^count))+7.5]/(sampling_rate/2),'stop');
        [b2,a2] = butter(6,[(stim_freq/(2^count))-40,(stim_freq/(2^count))+40]/(sampling_rate/2),'stop');

        for i = 1:length(recording_channels)
            if contains(filt_type{j,i},'stimTherapy')
                if ~isempty(wide_band_stim_freq_keys) && sum(wide_band_stim_freq_keys+1 == i) == 1
                    filt_data{j}.(recording_channels{i}) = filtfilt(b2,a2,filt_data{j}.(recording_channels{i}));
                else
                    filt_data{j}.(recording_channels{i}) = filtfilt(b,a,filt_data{j}.(recording_channels{i}));
                end
                stim_filtered = true;
            end
        end

        % Show the raw LFP data and power spectral density
        if show_plots && stim_filtered
            % Show the filtered LFP data and power spectral density
            figure;
            for i = 1:length(recording_channels)
                subplot(2,2,i);
                plot(filt_data{j}.(recording_channels{i}));
                ylabel('mV');
                xlabel('Data Sample')
                title(recording_channels{i});
            end
            if gait_hf_filtered && stim_sync_filtered
                sgtitle('Gait, Stim Sync, and Therapeutic Stim Filtered LFP');
            elseif ~gait_hf_filtered && stim_sync_filtered
                sgtitle('Stim Sync and Therapeutic Stim Filtered LFP');
            elseif gait_hf_filtered && ~stim_sync_filtered
                sgtitle('Gait and Therapeutic Stim Filtered LFP');
            else
                sgtitle('Therapeutic Stim Filtered LFP');
            end

            figure;
            for i = 1:length(recording_channels)
                subplot(2,2,i);
                pwelch(filt_data{j}.(recording_channels{i}),sampling_rate,round(sampling_rate*0.9),[1:200],sampling_rate);
                ylabel('dB/Hz');
                xlabel('Frequency')
                title(recording_channels{i});
            end
            if gait_hf_filtered && stim_sync_filtered
                sgtitle('Gait, Stim Sync, and Therapeutic Stim Filtered LFP');
            elseif ~gait_hf_filtered && stim_sync_filtered
                sgtitle('Stim Sync and Therapeutic Stim Filtered LFP');
            elseif gait_hf_filtered && ~stim_sync_filtered
                sgtitle('Gait and Therapeutic Stim Filtered LFP');
            else
                sgtitle('Therapeutic Stim Filtered LFP');
            end
        end

        if remove_stim_harm && count == 0
            count = count + 1;
        else
            stim_filt_done = true;
        end
    end

    %% Remove EKG artifacts
    for i = 1:size(ekg_template_indexes{j},1)
        if contains(filt_type{j,i},'ekg')
            if ~isnan(ekg_template_indexes{j}(i,1)) && ~isnan(ekg_template_indexes{j}(i,2))
                ekg_template = filt_data{j}.(recording_channels{i})(ekg_template_indexes{j}(i,1):ekg_template_indexes{j}(i,2));
                [sigClean, artRem, artLog, template] = RemoveNoiseTempMatch(filt_data{j}.(recording_channels{i}), sampling_rate, ekg_template, [1 5], [0,0], [], 1, 1, 1);
                filt_data{j}.(recording_channels{i}) = sigClean;
                ekg_filtered = true;
            end
        end
    end

    % Show the EKG and stim filtered LFP data and power spectral density
    if show_plots && ekg_filtered
        figure;
        for i = 1:length(recording_channels)
            subplot(2,2,i);
            plot(filt_data{j}.(recording_channels{i}));
            ylabel('mV');
            xlabel('Data Sample')
            title(recording_channels{i});
        end
        if gait_hf_filtered && stim_sync_filtered && stim_filtered
            sgtitle('Gait, Stim Sync, Therapeutic Stim, and EKG Filtered LFP');
        elseif ~gait_hf_filtered && stim_sync_filtered && stim_filtered
            sgtitle('Stim Sync, Therapeutic Stim, and EKG Filtered LFP');
        elseif gait_hf_filtered && ~stim_sync_filtered && stim_filtered
            sgtitle('Gait, Therapeutic Stim, and EKG Filtered LFP');
        elseif gait_hf_filtered && stim_sync_filtered && ~stim_filtered
            sgtitle('Gait, Stim Sync, and EKG Filtered LFP');
        elseif ~gait_hf_filtered && ~stim_sync_filtered && stim_filtered
            sgtitle('Therapeutic Stim and EKG Filtered LFP');
        elseif ~gait_hf_filtered && stim_sync_filtered && ~stim_filtered
            sgtitle('Stim Sync and EKG Filtered LFP');
        elseif gait_hf_filtered && ~stim_sync_filtered && ~stim_filtered
            sgtitle('Gait and EKG Filtered LFP');
        else
            sgtitle('EKG Filtered LFP');
        end
  

        figure;
        for i = 1:length(recording_channels)
            subplot(2,2,i);
            pwelch(filt_data{j}.(recording_channels{i}),sampling_rate,round(sampling_rate*0.9),[1:200],sampling_rate);
            ylabel('dB/Hz');
            xlabel('Frequency')
            title(recording_channels{i});
        end
        if gait_hf_filtered && stim_sync_filtered && stim_filtered
            sgtitle('Gait, Stim Sync, Therapeutic Stim, and EKG Filtered LFP');
        elseif ~gait_hf_filtered && stim_sync_filtered && stim_filtered
            sgtitle('Stim Sync, Therapeutic Stim, and EKG Filtered LFP');
        elseif gait_hf_filtered && ~stim_sync_filtered && stim_filtered
            sgtitle('Gait, Therapeutic Stim, and EKG Filtered LFP');
        elseif gait_hf_filtered && stim_sync_filtered && ~stim_filtered
            sgtitle('Gait, Stim Sync, and EKG Filtered LFP');
        elseif ~gait_hf_filtered && ~stim_sync_filtered && stim_filtered
            sgtitle('Therapeutic Stim and EKG Filtered LFP');
        elseif ~gait_hf_filtered && stim_sync_filtered && ~stim_filtered
            sgtitle('Stim Sync and EKG Filtered LFP');
        elseif gait_hf_filtered && ~stim_sync_filtered && ~stim_filtered
            sgtitle('Gait and EKG Filtered LFP');
        else
            sgtitle('EKG Filtered LFP');
        end
    end
end

%% Save data if the file path was originally passed in
if save_data
    if istable(input)
        file_path2 = [file_path(1:end-4),'_original.mat'];
        copyfile(file_path,file_path2);

        % Copy filt data as aligned_data to replace original file
        aligned_data = data.aligned_data;
        aligned_data.left_LFP_table = filt_data{1};
        aligned_data.right_LFP_table = filt_data{2};
        save(file_path,"aligned_data");
    else
        file_path2 = strrep(file_path,'TD.mat','TD_original.mat');
        copyfile(file_path,file_path2);

        % Copy filt data as timeDomainDataTable to replicate original file
        timeDomainDataTable = filt_data{1};
        [A,B,C] = fileparts(file_path);
        save_name = fullfile(A,[B,'_filtered',C]);
        save(save_name,'timeDomainDataTable');
    end
end
end