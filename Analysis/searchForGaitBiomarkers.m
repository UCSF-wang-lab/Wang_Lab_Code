function varargout = searchForGaitBiomarkers(aligned_data,signal_analysis_data,freq_lim,event_compare,subjectID,save_flag)
if ~exist('freq_lim','var') || isempty(freq_lim)
    freq_lim = [0 50];
end

if ~exist('event_compare','var') || isempty(event_compare)
    event_compare{1} = {'LTO','RTO'};
    event_compare{2} = {'LHS','RHS'};
end

if ~exist('subjectID','var') || isempty(subjectID)
    subjectID = 'RCSXX';
end

if ~exist('save_flag','var') || isempty(save_flag)
    save_flag = 0;
end

%% Extract data
if isfield(signal_analysis_data,'Left')
    left_sr = uniquetol(aligned_data.DeviceSettings.Left.timeDomainSettings.samplingRate,1);
    time_res_left = uniquetol(diff(signal_analysis_data.Left.Time{1}),1);
    if isfield(signal_analysis_data.Left,'PSD')
        analysis_type = 'FT';
    else
        analysis_type = 'CWT';
    end
    
    event_pairs = nchoosek(aligned_data.gait_events.Properties.VariableNames,2);
    p_val_matrix.Left = cell(1,length(signal_analysis_data.Left.Chan_Names));
    freq_bin_inds = genFreqBinPairs(signal_analysis_data.Left.Freq_Values{1},freq_lim);
    for i = 1:length(signal_analysis_data.Left.Chan_Names)
        temp_mat = nan(sum(diff(freq_bin_inds,1,2)==0),sum(diff(freq_bin_inds,1,2)==0),size(event_pairs,1));
        vals_power = [];
        vals_phase = [];
        for j = 1:length(aligned_data.gait_events.Properties.VariableNames)
            if ~isfield(vals_power,aligned_data.gait_events.Properties.VariableNames{j})
                vals_power.(aligned_data.gait_events.Properties.VariableNames{j}) = [];
            end
            count = 1;
            for k = 1:height(aligned_data.gait_events)
                event_time = aligned_data.gait_events.(aligned_data.gait_events.Properties.VariableNames{j})(k);
                if ~isnan(event_time)
                    [~,event_ind] = min(abs(signal_analysis_data.Left.Time{i}-event_time));
                    if isfield(signal_analysis_data.Left,'PSD')
                        temp = 20*log10(abs(signal_analysis_data.Left.Values{i}(:,event_ind)));
                        temp1 = angle(signal_analysis_data.Left.Values{i}(:,event_ind));
                    else
                        temp = abs(signal_analysis_data.Left.Values{i}(:,event_ind));
                        temp1 = angle(signal_analysis_data.Left.Values{i}(:,event_ind));
                    end
                    
                    if sum(isinf(temp),'all') == 0
                        for m = 1:size(freq_bin_inds,1)
                            vals_power.(aligned_data.gait_events.Properties.VariableNames{j})(count,m) = mean(temp(freq_bin_inds(m,1):freq_bin_inds(m,2)));
                            vals_phase.(aligned_data.gait_events.Properties.VariableNames{j})(count,m) = mean(temp1(freq_bin_inds(m,1):freq_bin_inds(m,2)));
                        end
                        count = count + 1;
                    end
                end
            end
        end
        
        for n = 1:size(event_pairs,1)
            for o = 1:size(freq_bin_inds,1)
                temp_mat(freq_bin_inds(o,1),freq_bin_inds(o,2),n) = vals_power.(event_pairs{n,1})(:,o);
%                 A = vals_power.LHS(:,o);
%                 B = vals_power.RTO(:,o);
%                 C = vals_power.RHS(:,o);
%                 D = vals_power.LTO(:,o);
%                 names = [repelem({'LHS'},length(A),1);repelem({'RTO'},length(B),1);repelem({'RHS'},length(C),1);repelem({'LTO'},length(D),1)];
%                 [~,~,stats] = anova1([A;B;C;D],names);
%                 [c,~,~,gnames] = multcompare(stats);
%                 [gnames(c(:,1)), gnames(c(:,2)), num2cell(c(:,3:6))]
            end
        end
        
    end
end

if isfield(signal_analysis_data,'Right')
    right_sr = uniquetol(aligned_data.DeviceSettings.Right.timeDomainSettings.samplingRate,1);
    time_res_right = uniquetol(diff(signal_analysis_data.Right.Time{1}),1);
    if isfield(signal_analysis_data.Right,'PSD')
        analysis_type = 'FT';
    else
        analysis_type = 'CWT';
    end
    power_at_event.Right = {};
    average_power.Right = {};
    average_power_time.Right = {};
    phase_at_event.Right = {};
    average_phase.Right = {};
    average_phase_time.Right = {};
    range_power.Right = {};
    range_power_time.Right = {};
    range_phase.Right = {};
    range_phase_time.Right = {};
    [~,band_names] = getFreqBandInd(signal_analysis_data.Right.Freq_Values{1});
    for i = 1:length(signal_analysis_data.Right.Chan_Names)
        band_inds = getFreqBandInd(signal_analysis_data.Right.Freq_Values{i});
        if band_inds(1,2) < band_inds(1,1)
            band_inds = fliplr(band_inds);
        end
        for j = 1:length(aligned_data.gait_events.Properties.VariableNames)
            if ~isfield(average_power.Right,aligned_data.gait_events.Properties.VariableNames{j})
                power_at_event.Right.(aligned_data.gait_events.Properties.VariableNames{j}) = cell(1,length(signal_analysis_data.Right.Chan_Names));
                average_power.Right.(aligned_data.gait_events.Properties.VariableNames{j}) = cell(1,length(signal_analysis_data.Right.Chan_Names));
                average_power_time.Right.(aligned_data.gait_events.Properties.VariableNames{j}) = cell(1,length(signal_analysis_data.Right.Chan_Names));
                phase_at_event.Right.(aligned_data.gait_events.Properties.VariableNames{j}) = cell(1,length(signal_analysis_data.Right.Chan_Names));
                average_phase.Right.(aligned_data.gait_events.Properties.VariableNames{j}) = cell(1,length(signal_analysis_data.Right.Chan_Names));
                average_phase_time.Right.(aligned_data.gait_events.Properties.VariableNames{j}) = cell(1,length(signal_analysis_data.Right.Chan_Names));
                range_power.Right.(aligned_data.gait_events.Properties.VariableNames{j}) = cell(1,length(signal_analysis_data.Right.Chan_Names));
                range_power_time.Right.(aligned_data.gait_events.Properties.VariableNames{j}) = cell(1,length(signal_analysis_data.Right.Chan_Names));
                range_phase.Right.(aligned_data.gait_events.Properties.VariableNames{j}) = cell(1,length(signal_analysis_data.Right.Chan_Names));
                range_phase_time.Right.(aligned_data.gait_events.Properties.VariableNames{j}) = cell(1,length(signal_analysis_data.Right.Chan_Names));
            end
            vals_power = [];
            vals_power_time = [];
            vals_phase = [];
            vals_phase_time = [];
            count = 1;
            for k = 1:height(aligned_data.gait_events)
                event_time = aligned_data.gait_events.(aligned_data.gait_events.Properties.VariableNames{j})(k);
                if ~isnan(event_time)
                    [~,min_ind_pre] = min(abs(signal_analysis_data.Right.Time{i}-(event_time-pre_post_time)));
                    [~,min_ind_post] = min(abs(signal_analysis_data.Right.Time{i}-(event_time+pre_post_time)));
                    [~,event_ind] = min(abs(signal_analysis_data.Right.Time{i}-event_time));
                    
                    if isfield(signal_analysis_data.Right,'PSD')
                        temp = 20*log10(abs(signal_analysis_data.Right.Values{i}(:,event_ind)));
                        temp1 = 20*log10(abs(signal_analysis_data.Right.Values{i}(:,min_ind_pre:min_ind_post)));
                        temp2 = angle(signal_analysis_data.Right.Values{i}(:,event_ind));
                        temp3 = angle(signal_analysis_data.Right.Values{i}(:,min_ind_pre:min_ind_post));
                    else
                        temp = 20*log10(abs(signal_analysis_data.Right.Values{i}(:,event_ind)));
                        temp1 = abs(signal_analysis_data.Right.Values{i}(:,min_ind_pre:min_ind_post));
                        temp2 = angle(signal_analysis_data.Right.Values{i}(:,event_ind));
                        temp3 = angle(signal_analysis_data.Right.Values{i}(:,min_ind_pre:min_ind_post));
                    end
                    
                    if size(temp1,2) == round((pre_post_time*2)/time_res_right)+1
                        if sum(isinf(temp),'all') == 0
                            for m = 1:length(band_names)
                                vals_power(count,:,m) = mean(temp(band_inds(m,1):band_inds(m,2),:));
                                vals_power_time(count,:,m) = mean(temp1(band_inds(m,1):band_inds(m,2),:));
                                vals_phase(count,:,m) = mean(temp2(band_inds(m,1):band_inds(m,2),:));
                                vals_phase_time(count,:,m) = mean(temp3(band_inds(m,1):band_inds(m,2),:));
                            end
                            count = count + 1;
                        end
                    end
                end
            end
            power_at_event.Right.(aligned_data.gait_events.Properties.VariableNames{j}){i} = squeeze(vals_power);
            average_power.Right.(aligned_data.gait_events.Properties.VariableNames{j}){i} = squeeze(mean(vals_power,1));
            average_power_time.Right.(aligned_data.gait_events.Properties.VariableNames{j}){i} = squeeze(mean(vals_power_time,1));
            phase_at_event.Right.(aligned_data.gait_events.Properties.VariableNames{j}){i} = squeeze(vals_phase);
            average_phase.Right.(aligned_data.gait_events.Properties.VariableNames{j}){i} = squeeze(mean(vals_phase,1));
            average_phase_time.Right.(aligned_data.gait_events.Properties.VariableNames{j}){i} = squeeze(mean(vals_phase_time,1));
            
            if strcmp(spread_type,'SE')
                range_power.Right.(aligned_data.gait_events.Properties.VariableNames{j}){i} = squeeze(std(vals_power,0,1))./size(vals_power,1);
                range_power_time.Right.(aligned_data.gait_events.Properties.VariableNames{j}){i} = squeeze(std(vals_power_time,0,1))./size(vals_power_time,1);
                range_phase.Right.(aligned_data.gait_events.Properties.VariableNames{j}){i} = squeeze(std(vals_phase,0,1))./size(vals_phase,1);
                range_phase_time.Right.(aligned_data.gait_events.Properties.VariableNames{j}){i} = squeeze(std(vals_phase_time,0,1))./size(vals_phase_time,1);
            else
                range_power.Right.(aligned_data.gait_events.Properties.VariableNames{j}){i} = squeeze(std(vals_power,0,1));
                range_power_time.Right.(aligned_data.gait_events.Properties.VariableNames{j}){i} = squeeze(std(vals_power_time,0,1));
                range_phase.Right.(aligned_data.gait_events.Properties.VariableNames{j}){i} = squeeze(std(vals_phase,0,1));
                range_phase_time.Right.(aligned_data.gait_events.Properties.VariableNames{j}){i} = squeeze(std(vals_phase_time,0,1));
            end
        end
    end
end
end

function freq_bin_inds = genFreqBinPairs(freq_vals,freq_lim)
ind_start = find(freq_vals >= freq_lim(1),1,'first');
ind_end = find(freq_vals <= freq_lim(2),1,'last');
n = ind_end - ind_start + 1;
skips = 1:n-1;
tot_pairs = (n^2-n)/2;
freq_bin_inds = nan(tot_pairs,2);

% Crazy indexing. Look up triangle numbers to recalculate these equations
f = @(x) n*(x+1)-(x^2+x)/2;
freq_bin_inds(1:n,:) = repmat([1:52]',1,2);
for i = 1:length(skips)
    freq_bin_inds(f(i-1)+1:f(i),1) = transpose(1:n-skips(i));
    freq_bin_inds(f(i-1)+1:f(i),2) = freq_bin_inds(f(i-1)+1:f(i),1) + skips(i);
    %         searchForGaitBiomarkers(aligned_data,A,[0,50],{{'LHS','RTO','RHS','LTO'}},'RCS03',0)
end


end