function [STFT,baseline] = calcRCS_STFT_Baseline(aligned_data,gapFillType,percentOverlap)
if ~exist('gapFillType','var') || isempty(gapFillType)
    gapFillType = 'blank';
end

if ~exist('percentOverlap','var') || isempty(percentOverlap)
    percentOverlap = 0.9;
end

%% Left RCS
if isfield(aligned_data,'left_LFP_table')
    % Spectrogram hyperparameters
    left_sr = aligned_data.DeviceSettings.Left.timeDomainSettings.samplingRate(end);
    WINDOW = left_sr;
    NOVERLAP = round(WINDOW*percentOverlap);
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
    
    figure('Name','Left RCS');
    ax(1) = subplot(3,1,1);
    plot(aligned_data.left_accel_taxis,aligned_data.left_Accel_table.XSamples);
    title({'Select start (first click) and end (second click) baseline time';'X Acce'});
    ax(2) = subplot(3,1,2);
    plot(aligned_data.left_accel_taxis,aligned_data.left_Accel_table.YSamples);
    title('Y Accel');
    ax(3) = subplot(3,1,3);
    plot(aligned_data.left_accel_taxis,aligned_data.left_Accel_table.ZSamples);
    title('Z Accel');
    linkaxes(ax,'x');
    xlim([0,10]);
    [baseline_times,~] = ginput(2);
    
    [~,baseline_start_ind] = cellfun(@(x)min(abs(x-baseline_times(1))),left_spect_time);
    [~,baseline_end_ind] = cellfun(@(x)min(abs(x-baseline_times(2))),left_spect_time);
    
    assert(all(baseline_start_ind==baseline_start_ind(1)),'Sampling rates may different between contacts, which is not possible.');
    assert(all(baseline_end_ind==baseline_end_ind(1)),'Sampling rates may different between contacts, which is not possible.');
    baseline_start_ind = baseline_start_ind(1);
    baseline_end_ind = baseline_end_ind(1);
    
    baseline.Left.Values = cellfun(@(x)x(:,baseline_start_ind:baseline_end_ind),left_spect,'UniformOutput',false);
    baseline.Left.Time = cellfun(@(x)x(:,baseline_start_ind:baseline_end_ind),left_spect_time,'UniformOutput',false);
    baseline.Left.Freq_Values = left_spect_freq;
    baseline.Left.PSD = cellfun(@(x)x(:,baseline_start_ind:baseline_end_ind),left_PSD,'UniformOutput',false);
    baseline.Left.Chan_Names = left_chan_names;
    
    STFT.Left.Values = left_spect;
    STFT.Left.Time = left_spect_time;
    STFT.Left.Freq_Values = left_spect_freq;
    STFT.Left.PSD = left_PSD;
    STFT.Left.Chan_Names = left_chan_names;
end

%% Right RCS
if isfield(aligned_data,'right_LFP_table')
    % Spectrogram hyperparameters
    right_sr = aligned_data.DeviceSettings.Right.timeDomainSettings.samplingRate(end);
    WINDOW = right_sr;
    NOVERLAP = round(WINDOW*percentOverlap);
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
    
    figure('Name','Right IPG');
    ax(1) = subplot(3,1,1);
    plot(aligned_data.right_accel_taxis,aligned_data.right_Accel_table.XSamples);
    title({'Select start (first click) and end (second click) baseline time';'X Acce'});
    ax(2) = subplot(3,1,2);
    plot(aligned_data.right_accel_taxis,aligned_data.right_Accel_table.YSamples);
    title('Y Accel');
    ax(3) = subplot(3,1,3);
    plot(aligned_data.right_accel_taxis,aligned_data.right_Accel_table.ZSamples);
    title('Z Accel');
    linkaxes(ax,'x');
    xlim([0,10]);
    [baseline_times,~] = ginput(2);
    
    [~,baseline_start_ind] = cellfun(@(x)min(abs(x-baseline_times(1))),right_spect_time);
    [~,baseline_end_ind] = cellfun(@(x)min(abs(x-baseline_times(2))),right_spect_time);
    
    assert(all(baseline_start_ind==baseline_start_ind(1)),'Sampling rates may different between contacts, which is not possible.');
    assert(all(baseline_end_ind==baseline_end_ind(1)),'Sampling rates may different between contacts, which is not possible.');
    baseline_start_ind = baseline_start_ind(1);
    baseline_end_ind = baseline_end_ind(1);
    
    baseline.Right.Values = cellfun(@(x)x(:,baseline_start_ind:baseline_end_ind),right_spect,'UniformOutput',false);
    baseline.Right.Time = cellfun(@(x)x(:,baseline_start_ind:baseline_end_ind),right_spect_time,'UniformOutput',false);
    baseline.Right.Freq_Values = left_spect_freq;
    baseline.Right.PSD = cellfun(@(x)x(:,baseline_start_ind:baseline_end_ind),right_PSD,'UniformOutput',false);
    baseline.Right.Chan_Names = right_chan_names;
    
    STFT.Right.Values = right_spect;
    STFT.Right.Time = right_spect_time;
    STFT.Right.Freq_Values = right_spect_freq;
    STFT.Right.PSD = right_PSD;
    STFT.Right.Chan_Names = right_chan_names;
end
end
