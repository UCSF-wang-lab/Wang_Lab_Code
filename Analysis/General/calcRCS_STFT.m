function calcRCS_STFT(filename,save_path,gapFillType)
% Load data
if ~exist('filename','var')
    [filename,path] = uigetfile('*w_Gait_Events*.mat');
    filename = fullfile(path,filename);
end
load(filename);

% Save location
if ~exist('savepath','var')
%     savepath = uigetdir();
end

% Drop packet fill technique (blank or interp)
if ~exist('gapFillType','var')
    gapFillType = 'blank';
end

A = strfind(filename,'RCS');
subjectID = filename(A(1):A(1)+4);

if contains(filename,"off_stim",'IgnoreCase',true)
    stim_cond = 'OFF Stim';
else
    stim_cond = 'ON Stim';
end

%% Left RCS
if isfield(aligned_data,'left_LFP_table')
    % Spectrogram hyperparameters
    left_sr = aligned_data.DeviceSettings.Left.timeDomainSettings.samplingRate(end);
    WINDOW = left_sr;
    NOVERLAP = WINDOW-1;
    NFFT = 2^nextpow2(WINDOW);
    
    chan_tag_inds = cellfun(@(x) contains(x,'chan'),aligned_data.DeviceSettings.Left.timeDomainSettings.Properties.VariableNames);
    chan_col_names = aligned_data.DeviceSettings.Left.timeDomainSettings.Properties.VariableNames(chan_tag_inds);
    left_chan_names = cellfun(@(x) aligned_data.DeviceSettings.Left.timeDomainSettings.(x){end},chan_col_names,'UniformOutput',false);
    
    left_spect = {};
    left_spect_freq = {};
    left_spect_time = {};
    left_PSD = {};
    remove_ind = [];
    for i = 1:length(left_chan_names)
        same_chan = cellfun(@(x) strcmp(left_chan_names{i},x),left_chan_names(1:i-1));
        if sum(same_chan) == 0  % Not a duplicate channel recording
            [data,time] = addEmptyData(aligned_data.left_taxis,aligned_data.left_LFP_table.(['key',num2str(i-1)]),left_sr,gapFillType);
            [left_spect{end+1},left_spect_freq{end+1},left_spect_time{end+1},left_PSD{end+1}]=spectrogram(data,WINDOW,NOVERLAP,NFFT,left_sr);
        else
            remove_ind = [remove_ind,i];
        end
    end
    left_chan_names(remove_ind) = [];
end

%% Right RCS
if isfield(aligned_data,'right_LFP_table')
    % Spectrogram hyperparameters
    right_sr = aligned_data.DeviceSettings.Right.timeDomainSettings.samplingRate(end);
    WINDOW = right_sr;
    NOVERLAP = WINDOW-1;
    NFFT = 2^nextpow2(WINDOW);
    
    chan_tag_inds = cellfun(@(x) contains(x,'chan'),aligned_data.DeviceSettings.Right.timeDomainSettings.Properties.VariableNames);
    chan_col_names = aligned_data.DeviceSettings.Right.timeDomainSettings.Properties.VariableNames(chan_tag_inds);
    right_chan_names = cellfun(@(x) aligned_data.DeviceSettings.Right.timeDomainSettings.(x){end},chan_col_names,'UniformOutput',false);
    
    right_spect = {};
    right_spect_freq = {};
    right_spect_time = {};
    right_PSD = {};
    remove_ind = [];
    for i = 1:length(right_chan_names)
        same_chan = cellfun(@(x) strcmp(right_chan_names{i},x),right_chan_names(1:i-1));
        if sum(same_chan) == 0  % Not a duplicate channel recording
            [data,time] = addEmptyData(aligned_data.right_taxis,aligned_data.right_LFP_table.(['key',num2str(i-1)]),right_sr,gapFillType);
            [right_spect{end+1},right_spect_freq{end+1},right_spect_time{end+1},right_PSD{end+1}]=spectrogram(data,WINDOW,NOVERLAP,NFFT,right_sr);
        else
            remove_ind = [remove_ind,i];
        end
    end
    right_chan_names(remove_ind) = [];
end
end