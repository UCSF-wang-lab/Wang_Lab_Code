function genSpectrogramGaitEvents(filename)
% Load data
if ~exist('filename','var')
    [filename,path] = uigetfile('*w_Gait_Events.mat');
    filename = fullfile(path,filename);
end
load(filename);

%% Left RCS
if isfield(aligned_data,'left_LFP_table')
    % Spectrogram hyperparameters
    left_sr = aligned_data.DeviceSettings.Left.timeDomainSettings.samplingRate(end);
    WINDOW = left_sr;
    NOVERLAP = WINDOW-1;
    NFFT = 2^nextpow2(WINDOW);
    
    chan_tag_inds = cellfun(@(x) contains(x,'chan'),aligned_data.DeviceSettings.Left.timeDomainSettings.Properties.VariableNames);
    chan_col_names = aligned_data.DeviceSettings.Left.timeDomainSettings.Properties.VariableNames(chan_tag_inds);
    chan_names = cellfun(@(x) aligned_data.DeviceSettings.Left.timeDomainSettings.(x){end},chan_col_names,'UniformOutput',false);
    
    left_spect = {};
    for i = 1:length(chan_names)
        same_chan = cellfun(@(x) strcmp(chan_names{i},x),chan_names(1:i-1));
        if sum(same_chan) == 0  % Not a duplicate channel recording
            [data,time] = addEmptyData(aligned_data.left_taxis,aligned_data.left_LFP_table.(['key',num2str(i-1)]),left_sr);
            [left_spect{end+1},left_spect_freq,left_spect_time]=spectrogram(data,WINDOW,NOVERLAP,NFFT,left_sr);
        end
    end
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
    chan_names = cellfun(@(x) aligned_data.DeviceSettings.Right.timeDomainSettings.(x){end},chan_col_names,'UniformOutput',false);
    
    right_spect = {};
    for i = 1:length(chan_names)
        same_chan = cellfun(@(x) strcmp(chan_names{i},x),chan_names(1:i-1));
        if sum(same_chan) == 0  % Not a duplicate channel recording
            [data,time] = addEmptyData(aligned_data.left_LFP_table.(['key',num2str(i-1)]),right_sr);
            [right_spect{end+1},right_spect_freq,right_spect_time]=spectrogram(data,WINDOW,NOVERLAP,NFFT,right_sr);
        end
    end
end

end

