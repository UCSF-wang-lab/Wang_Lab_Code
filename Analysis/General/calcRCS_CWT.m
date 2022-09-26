function CWT = calcRCS_CWT(aligned_data,gapFillType,varargin)
if ~exist('gapFillType','var')
    gapFillType = 'blank';
end

for i = 1:2:nargin-3
    switch varargin{i}
        case 'stim'
            if isstring(varargin{i+1})
                if strcmpi(varargin{i+1},'on')
                    stimOn = 1;
                else
                    stimOn = 0;
                end
            elseif isnumeric(varargin{i+1}) || islogical(varargin{i+1})
                stimOn = varargin{i+1};
            else
                stimOn = [];
            end
        case 'stimFreq'
            stimFreq = varargin{i+1};
    end
end

if ~exist('stimOn','var') || isempty(stimOn)
    stimOn = 0;
end

if ~exist('stimFreq','var')
    stimFreq = 130.2;
end

%% Basic filtering
if stimOn
    if isfield(aligned_data,'left_lFP_table')
        left_sr = aligned_data.DeviceSettings.Left.timeDomainSettings.samplingRate(end);
        [b,a] = butter(4,120/(left_sr/2));
        [d,c] = butter(2,stimFreq/(left_sr/2));
    
        chan_names_ind = cellfun(@(x)contains(x,'key'),fields(aligned_data.left_LFP_table));
        chan_names = aligned_data.left_LFP_table.Properties.VariableNames(chan_names_ind);
        
        for i = 1:length(chan_names)
            filt_stage1 = filtfilt(b,a,aligned_data.left_LFP_table.(chan_names{i}));
            filt_stage2 = filtfilt(d,c,filt_stage2);
            aligned_data.left_LFP_table.(chan_names{i}) = filt_stage2;
        end
    end
    
    if isfield(aligned_data,'right_LFP_table')
        right_sr = aligned_data.DeviceSettings.Right.timeDomainSettings.samplingRate(end);
        [b,a] = butter(4,120/(right_sr/2));
        [d,c] = butter(2,stimFreq/(right_sr/2));
    
        chan_names_ind = cellfun(@(x)contains(x,'key'),fields(aligned_data.right_LFP_table));
        chan_names = aligned_data.right_LFP_table.Properties.VariableNames(chan_names_ind);
        
        for i = 1:length(chan_names)
            filt_stage1 = filtfilt(b,a,aligned_data.right_LFP_table.(chan_names{i}));
            filt_stage2 = filtfilt(d,c,filt_stage2);
            aligned_data.right_LFP_table.(chan_names{i}) = filt_stage2;
        end
    end
end

%% Left RCS
if isfield(aligned_data,'left_LFP_table')
    % Spectrogram hyperparameters
    left_sr = aligned_data.DeviceSettings.Left.timeDomainSettings.samplingRate(end);
    if ~isnumeric(left_sr)
        left_sr = 500;
    end
    
    chan_tag_inds = cellfun(@(x) contains(x,'chan'),aligned_data.DeviceSettings.Left.timeDomainSettings.Properties.VariableNames);
    chan_col_names = aligned_data.DeviceSettings.Left.timeDomainSettings.Properties.VariableNames(chan_tag_inds);
    left_chan_names = cellfun(@(x) aligned_data.DeviceSettings.Left.timeDomainSettings.(x){end},chan_col_names,'UniformOutput',false);
    
    left_cwt = {};
    left_cwt_freq = {};
    left_cwt_time = {};
    remove_ind = [];
    for i = 1:length(left_chan_names)
        same_chan = cellfun(@(x) strcmp(left_chan_names{i},x),left_chan_names(1:i-1));
        if sum(same_chan) == 0  % Not a duplicate channel recording
            [data,left_cwt_time{end+1}] = addEmptyData(aligned_data.left_taxis,aligned_data.left_LFP_table.(['key',num2str(i-1)]),left_sr,gapFillType);

            if sum(isnan(data)) == 0
                [left_cwt{end+1},left_cwt_freq{end+1}]=cwt(data,left_sr);
            else
                left_cwt{end+1} = nan(length(left_cwt_freq{end}),length(data));
                left_cwt_freq{end+1} = nan(length(left_cwt_freq{end}),1);
            end
        else
            remove_ind = [remove_ind,i];
        end
    end
    left_chan_names(remove_ind) = [];
    
    CWT.Left.Values = left_cwt;
    CWT.Left.Time = left_cwt_time;
    CWT.Left.Freq_Values = left_cwt_freq;
    CWT.Left.Chan_Names = left_chan_names;
end

%% Right RCS
if isfield(aligned_data,'right_LFP_table')
    % Spectrogram hyperparameters
    right_sr = aligned_data.DeviceSettings.Right.timeDomainSettings.samplingRate(end);
    if ~isnumeric(right_sr)
        right_sr = 500;
    end
    
    chan_tag_inds = cellfun(@(x) contains(x,'chan'),aligned_data.DeviceSettings.Right.timeDomainSettings.Properties.VariableNames);
    chan_col_names = aligned_data.DeviceSettings.Right.timeDomainSettings.Properties.VariableNames(chan_tag_inds);
    right_chan_names = cellfun(@(x) aligned_data.DeviceSettings.Right.timeDomainSettings.(x){end},chan_col_names,'UniformOutput',false);
    
    right_cwt = {};
    right_cwt_freq = {};
    right_cwt_time = {};
    remove_ind = [];
    for i = 1:length(right_chan_names)
        same_chan = cellfun(@(x) strcmp(right_chan_names{i},x),right_chan_names(1:i-1));
        if sum(same_chan) == 0  % Not a duplicate channel recording
            [data,right_cwt_time{end+1}] = addEmptyData(aligned_data.right_taxis,aligned_data.right_LFP_table.(['key',num2str(i-1)]),right_sr,gapFillType);

            if sum(isnan(data)) == 0
                [right_cwt{end+1},right_cwt_freq{end+1}]=cwt(data,right_sr);
            else
                right_cwt{end+1} = nan(length(right_cwt_freq{end}),length(data));
                right_cwt_freq{end+1} = nan(length(right_cwt_freq{end}),1);
            end
        else
            remove_ind = [remove_ind,i];
        end
    end
    right_chan_names(remove_ind) = [];
    
    CWT.Right.Values = right_cwt;
    CWT.Right.Time = right_cwt_time;
    CWT.Right.Freq_Values = right_cwt_freq;
    CWT.Right.Chan_Names = right_chan_names;
end
end