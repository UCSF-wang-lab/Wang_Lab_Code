function varargout = genAvgBandPowerPhaseGaitEvents(aligned_data,signal_analysis_data,pre_post_time,event_compare,spread_type,subjectID,plot_type,save_flag)

if ~exist('pre_post_time','var') || isempty(pre_post_time)
    pre_post_time = 1;
end

if ~exist('event_compare','var') || isempty(event_compare)
    event_compare{1} = {'LTO','RTO'};
    event_compare{2} = {'LHS','RHS'};
end

if ~exist('spread_type','var') || isempty(spread_type)
    spread_type = 'SD';
end

if ~exist('subjectID','var') || isempty(subjectID)
    subjectID = 'RCSXX';
end

if ~exist('plot_type','var') || isempty(plot_type)
    plot_type = 'all';  % all, power, power_boxplot, power_boxplot_scatter, phase, phase_boxplot, phase_boxplot_scatter, phase_polar
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
    power_at_event.Left = {};
    average_power.Left = {};
    average_power_time.Left = {};
    phase_at_event.Left = {};
    average_phase.Left = {};
    average_phase_time.Left = {};
    range_power.Left = {};
    range_power_time.Left = {};
    range_phase.Left = {};
    range_phase_time.Left = {};
    [~,band_names] = getFreqBandInd(signal_analysis_data.Left.Freq_Values{1});
    for i = 1:length(signal_analysis_data.Left.Chan_Names)
        band_inds = getFreqBandInd(signal_analysis_data.Left.Freq_Values{i});
        if band_inds(1,2) < band_inds(1,1)
            band_inds = fliplr(band_inds);
        end
        for j = 1:length(aligned_data.gait_events.Properties.VariableNames)
            if ~isfield(average_power.Left,aligned_data.gait_events.Properties.VariableNames{j})
                power_at_event.Left.(aligned_data.gait_events.Properties.VariableNames{j}) = cell(1,length(signal_analysis_data.Left.Chan_Names));
                average_power.Left.(aligned_data.gait_events.Properties.VariableNames{j}) = cell(1,length(signal_analysis_data.Left.Chan_Names));
                average_power_time.Left.(aligned_data.gait_events.Properties.VariableNames{j}) = cell(1,length(signal_analysis_data.Left.Chan_Names));
                phase_at_event.Left.(aligned_data.gait_events.Properties.VariableNames{j}) = cell(1,length(signal_analysis_data.Left.Chan_Names));
                average_phase.Left.(aligned_data.gait_events.Properties.VariableNames{j}) = cell(1,length(signal_analysis_data.Left.Chan_Names));
                average_phase_time.Left.(aligned_data.gait_events.Properties.VariableNames{j}) = cell(1,length(signal_analysis_data.Left.Chan_Names));
                range_power.Left.(aligned_data.gait_events.Properties.VariableNames{j}) = cell(1,length(signal_analysis_data.Left.Chan_Names));
                range_power_time.Left.(aligned_data.gait_events.Properties.VariableNames{j}) = cell(1,length(signal_analysis_data.Left.Chan_Names));
                range_phase.Left.(aligned_data.gait_events.Properties.VariableNames{j}) = cell(1,length(signal_analysis_data.Left.Chan_Names));
                range_phase_time.Left.(aligned_data.gait_events.Properties.VariableNames{j}) = cell(1,length(signal_analysis_data.Left.Chan_Names));
            end
            vals_power = [];
            vals_power_time = [];
            vals_phase = [];
            vals_phase_time = [];
            count = 1;
            count_time = 1;
            for k = 1:height(aligned_data.gait_events)
                event_time = aligned_data.gait_events.(aligned_data.gait_events.Properties.VariableNames{j})(k);
                if ~isnan(event_time)
                    [~,min_ind_pre] = min(abs(signal_analysis_data.Left.Time{i}-(event_time-pre_post_time)));
                    [~,min_ind_post] = min(abs(signal_analysis_data.Left.Time{i}-(event_time+pre_post_time)));
                    [~,event_ind] = min(abs(signal_analysis_data.Left.Time{i}-event_time));
                    if isfield(signal_analysis_data.Left,'PSD')
                        temp = 20*log10(abs(signal_analysis_data.Left.Values{i}(:,event_ind)));
                        temp1 = 20*log10(abs(signal_analysis_data.Left.Values{i}(:,min_ind_pre:min_ind_post)));
                        temp2 = angle(signal_analysis_data.Left.Values{i}(:,event_ind));
                        temp3 = angle(signal_analysis_data.Left.Values{i}(:,min_ind_pre:min_ind_post));
                    else
                        temp = abs(signal_analysis_data.Left.Values{i}(:,event_ind));
                        temp1 = abs(signal_analysis_data.Left.Values{i}(:,min_ind_pre:min_ind_post));
                        temp2 = angle(signal_analysis_data.Left.Values{i}(:,event_ind));
                        temp3 = angle(signal_analysis_data.Left.Values{i}(:,min_ind_pre:min_ind_post));
                    end
                    
                    if size(temp1,2) == round((pre_post_time*2)/time_res_left)+1
                        if sum(isinf(temp1),'all') == 0
                            for m = 1:length(band_names)
                                vals_power_time(count_time,:,m) = mean(temp1(band_inds(m,1):band_inds(m,2),:));
                                vals_phase_time(count_time,:,m) = mean(temp3(band_inds(m,1):band_inds(m,2),:));
                                count_time = count_time + 1;
                            end
                            
                        end
                        if sum(isinf(temp),'all') == 0
                            for m = 1:length(band_names)
                                vals_power(count,:,m) = mean(temp(band_inds(m,1):band_inds(m,2),:));
                                vals_phase(count,:,m) = mean(temp2(band_inds(m,1):band_inds(m,2),:));
                            end
                            count = count + 1;
                        end
                    end
                end
            end
            power_at_event.Left.(aligned_data.gait_events.Properties.VariableNames{j}){i} = squeeze(vals_power);
            average_power.Left.(aligned_data.gait_events.Properties.VariableNames{j}){i} = squeeze(mean(vals_power,1));
            average_power_time.Left.(aligned_data.gait_events.Properties.VariableNames{j}){i} = squeeze(mean(vals_power_time,1));
            phase_at_event.Left.(aligned_data.gait_events.Properties.VariableNames{j}){i} = squeeze(vals_phase);
            average_phase.Left.(aligned_data.gait_events.Properties.VariableNames{j}){i} = squeeze(mean(vals_phase,1));
            average_phase_time.Left.(aligned_data.gait_events.Properties.VariableNames{j}){i} = squeeze(mean(vals_phase_time,1));
            
            if strcmp(spread_type,'SE')
                range_power.Left.(aligned_data.gait_events.Properties.VariableNames{j}){i} = squeeze(std(vals_power,0,1))./size(vals_power,1);
                range_power_time.Left.(aligned_data.gait_events.Properties.VariableNames{j}){i} = squeeze(std(vals_power_time,0,1))./size(vals_power_time,1);
                range_phase.Left.(aligned_data.gait_events.Properties.VariableNames{j}){i} = squeeze(std(vals_phase,0,1))./size(vals_phase,1);
                range_phase_time.Left.(aligned_data.gait_events.Properties.VariableNames{j}){i} = squeeze(std(vals_phase_time,0,1))./size(vals_phase_time,1);
            else
                range_power.Left.(aligned_data.gait_events.Properties.VariableNames{j}){i} = squeeze(std(vals_power,0,1));
                range_power_time.Left.(aligned_data.gait_events.Properties.VariableNames{j}){i} = squeeze(std(vals_power_time,0,1));
                range_phase.Left.(aligned_data.gait_events.Properties.VariableNames{j}){i} = squeeze(std(vals_phase,0,1));
                range_phase_time.Left.(aligned_data.gait_events.Properties.VariableNames{j}){i} = squeeze(std(vals_phase_time,0,1));
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
            count_time = 1;
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
                        if sum(isinf(temp1),'all') == 0
                            for m = 1:length(band_names)
                                vals_power_time(count_time,:,m) = mean(temp1(band_inds(m,1):band_inds(m,2),:));
                                vals_phase_time(count_time,:,m) = mean(temp3(band_inds(m,1):band_inds(m,2),:));
                            end
                            count_time = count_time + 1;
                        end
                        
                        if sum(isinf(temp),'all') == 0
                            for m = 1:length(band_names)
                                vals_power(count,:,m) = mean(temp(band_inds(m,1):band_inds(m,2),:));
                                vals_phase(count,:,m) = mean(temp2(band_inds(m,1):band_inds(m,2),:));
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

%% Plot
fig_vec = [];
colors = CBMap('GaitEvents',4);

if strcmp(plot_type,'all') || strcmp(plot_type,'power')
    if isfield(average_power_time,'Left')
        time_vec = -1:time_res_left:1;
        for i = 1:length(signal_analysis_data.Left.Chan_Names)
            for j = 1:length(event_compare)
                fig_vec(end+1) = figure;
                for k = 1:length(band_names)
                    ax_hand = subplot(2,4,k);
                    for m = 1:length(event_compare{j})
                        mean_std_plot(time_vec,average_power_time.Left.(event_compare{j}{m}){i}(:,k),range_power_time.Left.(event_compare{j}{m}){i}(:,k),ax_hand,colors.(event_compare{j}{m}),[]);
                        %                     mean_std_plot(time_vec,average_phase.Left.(event_compare{j}{m}){i}(:,k),std_phase.Left.(event_compare{j}{m}){i}(:,k),ax_hand,colors.(event_compare{j}{m}),[]);
                    end
                    hold(ax_hand,'on');
                    xline(0,'--k');
                    hold(ax_hand,'off');
                    xlabel('Time (s)');
                    
                    if isfield(signal_analysis_data.Left,'PSD')
                        ylabel('\muV^2'); %  uV^2
                    else
                        ylabel('Magnitude');
                    end
                    
                    title(band_names{k});
                end
                sgtitle({[subjectID,' Left'];signal_analysis_data.Left.Chan_Names{i};createPlotTitle(event_compare{j})});
            end
        end
    end
    
    if isfield(average_power_time,'Right')
        time_vec = -1:time_res_right:1;
        for i = 1:length(signal_analysis_data.Right.Chan_Names)
            for j = 1:length(event_compare)
                fig_vec(end+1) = figure;
                for k = 1:length(band_names)
                    ax_hand = subplot(2,4,k);
                    for m = 1:length(event_compare{j})
                        mean_std_plot(time_vec,average_power_time.Right.(event_compare{j}{m}){i}(:,k),range_power_time.Right.(event_compare{j}{m}){i}(:,k),ax_hand,colors.(event_compare{j}{m}),[]);
                    end
                    hold(ax_hand,'on');
                    xline(0,'--k');
                    hold(ax_hand,'off');
                    xlabel('Time (s)');
                    
                    if isfield(signal_analysis_data.Right,'PSD')
                        ylabel('\muV^2'); % uV^2
                    else
                        ylabel('Magnitude');
                    end
                    title(band_names{k});
                end
                sgtitle({[subjectID,' Right'];signal_analysis_data.Right.Chan_Names{i};createPlotTitle(event_compare{j})});
            end
        end
    end
end

if strcmp(plot_type,'all') || strcmp(plot_type,'phase')
    if isfield(average_phase_time,'Left')
        time_vec = -1:time_res_left:1;
        for i = 1:length(signal_analysis_data.Left.Chan_Names)
            for j = 1:length(event_compare)
                fig_vec(end+1) = figure;
                for k = 1:length(band_names)
                    ax_hand = subplot(2,4,k);
                    for m = 1:length(event_compare{j})
                        mean_std_plot(time_vec,average_phase_time.Left.(event_compare{j}{m}){i}(:,k),range_phase_time.Left.(event_compare{j}{m}){i}(:,k),ax_hand,colors.(event_compare{j}{m}),[]);
                    end
                    hold(ax_hand,'on');
                    xline(0,'--k');
                    hold(ax_hand,'off');
                    xlabel('Time (s)');
                    
                    if isfield(signal_analysis_data.Left,'PSD')
                        ylabel('Phase');
                    else
                        ylabel('Phase');
                    end
                    
                    title(band_names{k});
                end
                sgtitle({[subjectID,' Left'];signal_analysis_data.Left.Chan_Names{i};createPlotTitle(event_compare{j})});
            end
        end
    end
    
    if isfield(average_phase_time,'Right')
        time_vec = -1:time_res_right:1;
        for i = 1:length(signal_analysis_data.Right.Chan_Names)
            for j = 1:length(event_compare)
                fig_vec(end+1) = figure;
                for k = 1:length(band_names)
                    ax_hand = subplot(2,4,k);
                    for m = 1:length(event_compare{j})
                        mean_std_plot(time_vec,average_phase_time.Right.(event_compare{j}{m}){i}(:,k),range_phase_time.Right.(event_compare{j}{m}){i}(:,k),ax_hand,colors.(event_compare{j}{m}),[]);
                    end
                    hold(ax_hand,'on');
                    xline(0,'--k');
                    hold(ax_hand,'off');
                    xlabel('Time (s)');
                    
                    if isfield(signal_analysis_data.Right,'PSD')
                        ylabel('Phase');
                    else
                        ylabel('Phase');
                    end
                    
                    title(band_names{k});
                end
                sgtitle({[subjectID,' Right'];signal_analysis_data.Right.Chan_Names{i};createPlotTitle(event_compare{j})});
            end
        end
    end
end

if strcmp(plot_type,'all') || contains(plot_type,'power_boxplot')
    if isfield(power_at_event,'Left')
        for i = 1:length(signal_analysis_data.Left.Chan_Names)
            for j = 1:length(event_compare)
                fig_vec(end+1) = figure;
                for k = 1:length(band_names)
                    ax_hand = subplot(2,4,k);
                    box_vals = [];
                    box_group = {};
                    for m = 1:length(event_compare{j})
                        box_vals = [box_vals;power_at_event.Left.(event_compare{j}{m}){i}(:,k)];
                        box_group = [box_group;repelem({event_compare{j}{m}},size(power_at_event.Left.(event_compare{j}{m}){i}(:,k),1),1)];
                    end
                    boxplot(box_vals,box_group,'Colors',[colors.LHS;colors.RTO;colors.RHS;colors.LTO],'Symbol','');
                    hold on;
                    if contains(plot_type,'scatter')
                        for n = 1:length(event_compare{j})
                            inds = strcmp(box_group,event_compare{j}{n});
                            scatter_vals = box_vals(inds);
                            scatter_pos = scatterPointJitter(length(scatter_vals),n-0.3,n+0.3);
                            scatter(scatter_pos,scatter_vals,'o','MarkerFacecolor',colors.(event_compare{j}{n}),'MarkerFaceAlpha',0.5,'MarkerEdgeColor','none');
                        end
                    end
                    
                    if isfield(signal_analysis_data.Left,'PSD')
                        ylabel('\muV^2');
                    else
                        ylabel('Magnitude');
                    end
                    title(band_names{k});
                end
                sgtitle({[subjectID,' Left'];signal_analysis_data.Left.Chan_Names{i};createPlotTitle(event_compare{j})});
            end
        end
    end
    
    if isfield(power_at_event,'Right')
        for i = 1:length(signal_analysis_data.Right.Chan_Names)
            for j = 1:length(event_compare)
                fig_vec(end+1) = figure;
                for k = 1:length(band_names)
                    ax_hand = subplot(2,4,k);
                    box_vals = [];
                    box_group = {};
                    for m = 1:length(event_compare{j})
                        box_vals = [box_vals;power_at_event.Right.(event_compare{j}{m}){i}(:,k)];
                        box_group = [box_group;repelem({event_compare{j}{m}},size(power_at_event.Right.(event_compare{j}{m}){i}(:,k),1),1)];
                    end
                    boxplot(box_vals,box_group,'Colors',[colors.LHS;colors.RTO;colors.RHS;colors.LTO],'Symbol','');
                    hold on;
                    if contains(plot_type,'scatter')
                        for n = 1:length(event_compare{j})
                            inds = strcmp(box_group,event_compare{j}{n});
                            scatter_vals = box_vals(inds);
                            scatter_pos = scatterPointJitter(length(scatter_vals),n-0.3,n+0.3);
                            scatter(scatter_pos,scatter_vals,'o','MarkerFacecolor',colors.(event_compare{j}{n}),'MarkerFaceAlpha',0.5,'MarkerEdgeColor','none');
                        end
                    end
                    if isfield(signal_analysis_data.Right,'PSD')
                        ylabel('\mV^2');
                    else
                        ylabel('Magnitude');
                    end
                    title(band_names{k});
                end
                sgtitle({[subjectID,' Right'];signal_analysis_data.Right.Chan_Names{i};createPlotTitle(event_compare{j})});
            end
        end
    end
end

if strcmp(plot_type,'all') || contains(plot_type,'phase_boxplot')
    if isfield(phase_at_event,'Left')
        for i = 1:length(signal_analysis_data.Left.Chan_Names)
            for j = 1:length(event_compare)
                fig_vec(end+1) = figure;
                for k = 1:length(band_names)
                    ax_hand = subplot(2,4,k);
                    box_vals = [];
                    box_group = {};
                    for m = 1:length(event_compare{j})
                        box_vals = [box_vals;phase_at_event.Left.(event_compare{j}{m}){i}(:,k)];
                        box_group = [box_group;repelem({event_compare{j}{m}},size(phase_at_event.Left.(event_compare{j}{m}){i}(:,k),1),1)];
                    end
                    boxplot(box_vals,box_group,'Colors',[colors.LHS;colors.RTO;colors.RHS;colors.LTO],'Symbol','');
                    hold on;
                    if contains(plot_type,'scatter')
                        for n = 1:length(event_compare{j})
                            inds = strcmp(box_group,event_compare{j}{n});
                            scatter_vals = box_vals(inds);
                            scatter_pos = scatterPointJitter(length(scatter_vals),n-0.3,n+0.3);
                            scatter(scatter_pos,scatter_vals,'o','MarkerFacecolor',colors.(event_compare{j}{n}),'MarkerFaceAlpha',0.5,'MarkerEdgeColor','none');
                        end
                    end
                    ylabel('Phase');
                    ylim([-pi,pi]);
                    yticks([-pi:pi/4:pi]);
                    yticklabels({'-\pi','-3\pi/4','-\pi/2','-\pi/4','0','\pi/4','\pi/2','3\pi/4','\pi'})
                    set(gca,'TickLabelInterpreter','tex')
                    title(band_names{k});
                end
                sgtitle({[subjectID,' Left'];signal_analysis_data.Left.Chan_Names{i};createPlotTitle(event_compare{j})});
            end
        end
    end
    
    if isfield(phase_at_event,'Right')
        for i = 1:length(signal_analysis_data.Right.Chan_Names)
            for j = 1:length(event_compare)
                fig_vec(end+1) = figure;
                for k = 1:length(band_names)
                    ax_hand = subplot(2,4,k);
                    box_vals = [];
                    box_group = {};
                    for m = 1:length(event_compare{j})
                        box_vals = [box_vals;phase_at_event.Right.(event_compare{j}{m}){i}(:,k)];
                        box_group = [box_group;repelem({event_compare{j}{m}},size(phase_at_event.Right.(event_compare{j}{m}){i}(:,k),1),1)];
                    end
                    boxplot(box_vals,box_group,'Colors',[colors.LHS;colors.RTO;colors.RHS;colors.LTO],'Symbol','');
                    hold on;
                    if contains(plot_type,'scatter')
                        for n = 1:length(event_compare{j})
                            inds = strcmp(box_group,event_compare{j}{n});
                            scatter_vals = box_vals(inds);
                            scatter_pos = scatterPointJitter(length(scatter_vals),n-0.3,n+0.3);
                            scatter(scatter_pos,scatter_vals,'o','MarkerFacecolor',colors.(event_compare{j}{n}),'MarkerFaceAlpha',0.5,'MarkerEdgeColor','none');
                        end
                    end
                    ylabel('Phase');
                    ylim([-pi,pi]);
                    yticks([-pi:pi/4:pi]);
                    yticklabels({'-\pi','-3\pi/4','-\pi/2','-\pi/4','0','\pi/4','\pi/2','3\pi/4','\pi'})
                    set(gca,'TickLabelInterpreter','tex')
                    title(band_names{k});
                end
                sgtitle({[subjectID,' Right'];signal_analysis_data.Right.Chan_Names{i};createPlotTitle(event_compare{j})});
            end
        end
    end
end

if strcmp(plot_type,'all') || strcmp(plot_type,'phase_polar')
    if isfield(phase_at_event,'Left')
        for i = 1:length(signal_analysis_data.Left.Chan_Names)
            for j = 1:length(event_compare)
                fig_vec(end+1) = figure;
                for k = 1:length(band_names)
                    ax_hand = subplot(2,4,k);
                    polar_handle = polaraxes('Units',ax_hand.Units,'Position',ax_hand.Position);
                    delete(ax_hand);
                    hold(polar_handle,'on');
                    for m = 1:length(event_compare{j})
                        polarhistogram(polar_handle,wrapTo2Pi(phase_at_event.Left.(event_compare{j}{m}){i}(:,k)),10,'FaceColor',colors.(event_compare{j}{m}),'FaceAlpha',0.3);
                    end
                    hold(polar_handle,'off');
                    title(band_names{k});
                end
                sgtitle({[subjectID,' Left'];signal_analysis_data.Left.Chan_Names{i};createPlotTitle(event_compare{j})});
            end
        end
    end
    
    if isfield(phase_at_event,'Right')
        for i = 1:length(signal_analysis_data.Right.Chan_Names)
            for j = 1:length(event_compare)
                fig_vec(end+1) = figure;
                for k = 1:length(band_names)
                    ax_hand = subplot(2,4,k);
                    polar_handle = polaraxes('Units',ax_hand.Units,'Position',ax_hand.Position);
                    delete(ax_hand);
                    hold(polar_handle,'on');
                    for m = 1:length(event_compare{j})
                        polarhistogram(polar_handle,phase_at_event.Right.(event_compare{j}{m}){i}(:,k),10,'FaceColor',colors.(event_compare{j}{m}),'FaceAlpha',0.3);
                    end
                    hold(polar_handle,'off');
                    title(band_names{k});
                end
                sgtitle({[subjectID,' Right'];signal_analysis_data.Right.Chan_Names{i};createPlotTitle(event_compare{j})});
            end
        end
    end
end

%% Save plots (NEED TO CHANGE; commented out for now)
if save_flag
    save_dir = uigetdir();
    
    figure_format(12,8,12);
    
    % check if saving folders exist
    if ~isfolder(fullfile(save_dir,'AvgEventPower'))
        mkdir(fullfile(save_dir,'AvgEventPower'));
    end
    
    if ~isfolder(fullfile(save_dir,'AvgEventPower',aligned_data.stim_condition))
        mkdir(fullfile(save_dir,'AvgEventPower',aligned_data.stim_condition))
    end
    
    if strcmp(analysis_type,'FT')
        if ~isfolder(fullfile(save_dir,'AvgEventPower',aligned_data.stim_condition,'FT'))
            mkdir(fullfile(save_dir,'AvgEventPower',aligned_data.stim_condition,'FT'))
        end
    elseif strcmp(analysis_type,'CWT')
        if ~isfolder(fullfile(save_dir,'AvgEventPower',aligned_data.stim_condition,'CWT'))
            mkdir(fullfile(save_dir,'AvgEventPower',aligned_data.stim_condition,'CWT'))
        end
    end
    
    folders_to_check = {'FIG_files','PDF_files','TIFF_files'};
    extension = {'.fig','.pdf','.tiff'};
    for n = 1:length(folders_to_check)
        if strcmp(analysis_type,'FT')
            if ~isfolder(fullfile(save_dir,'AvgEventPower',aligned_data.stim_condition,'FT',folders_to_check{n}))
                mkdir(fullfile(save_dir,'AvgEventPower',aligned_data.stim_condition,'FT',folders_to_check{n}));
            end
        elseif strcmp(analysis_type,'CWT')
            if ~isfolder(fullfile(save_dir,'AvgEventPower',aligned_data.stim_condition,'CWT',folders_to_check{n}))
                mkdir(fullfile(save_dir,'AvgEventPower',aligned_data.stim_condition,'CWT',folders_to_check{n}));
            end
        end
    end
    
    for i = 1:length(fig_vec)
        curr_axes = gca(fig_vec(i));
        save_name = [];
        for j = 1:length(curr_axes.Parent.Children(1).String)
            if isempty(save_name)
                save_name = curr_axes.Parent.Children(1).String{j};
            else
                save_name = [save_name,' ', curr_axes.Parent.Children(1).String{j}];
            end
        end
        
        if contains(plot_type,'boxplot') || contains(plot_type,'phase_polar')
            save_name = [save_name,' ', plot_type];
        else
            save_name = [save_name,' ', plot_type, ' ', spread_type];
        end
        
        if strcmp(analysis_type,'FT')
            savefig(fig_vec(i),fullfile(save_dir,'AvgEventPower',aligned_data.stim_condition,'FT',folders_to_check{1},strrep(strrep(save_name,' ','_'),'.','')));
        elseif strcmp(analysis_type,'CWT')
            savefig(fig_vec(i),fullfile(save_dir,'AvgEventPower',aligned_data.stim_condition,'CWT',folders_to_check{1},strrep(strrep(save_name,' ','_'),'.','')));
        end
        
        for k = 2:length(folders_to_check)
            if strcmp(analysis_type,'FT')
                print(fig_vec(i),[fullfile(save_dir,'AvgEventPower',aligned_data.stim_condition,'FT',folders_to_check{k},strrep(strrep(save_name,' ','_'),'.','')),extension{k}],'-r300',['-d',extension{k}(2:end)]);
            elseif strcmp(analysis_type,'CWT')
                print(fig_vec(i),[fullfile(save_dir,'AvgEventPower',aligned_data.stim_condition,'CWT',folders_to_check{k},strrep(strrep(save_name,' ','_'),'.','')),extension{k}],'-r300',['-d',extension{k}(2:end)]);
            end
        end
    end
end

varargout = {power_at_event, average_power, average_power_time,...
    phase_at_event, average_phase, average_phase_time,...
    range_power, range_power_time,...
    range_phase, range_phase_time};
end

function plot_title = createPlotTitle(s)
plot_title = s{1};

for i = 2:length(s)
    plot_title = plot_title + " vs. " + s{i};
end
end