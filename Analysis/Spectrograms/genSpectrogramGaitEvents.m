function genSpectrogramGaitEvents(filename,savepath,gapFillType)
% Load data
if ~exist('filename','var')
    [filename,path] = uigetfile('*w_Gait_Events.mat');
    filename = fullfile(path,filename);
end
load(filename);

% Save location
if ~exist('savepath','var')
    savepath = uigetdir();
end

% Drop packet fill technique (blank or interp)
if ~exist('gapFillType')
    gapFillType = 'blank';
end

A = strfind(filename,'RCS');
subjectID = filename(A(1):A(1)+4);

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
    for i = 1:length(left_chan_names)
        same_chan = cellfun(@(x) strcmp(left_chan_names{i},x),left_chan_names(1:i-1));
        if sum(same_chan) == 0  % Not a duplicate channel recording
            [data,time] = addEmptyData(aligned_data.left_taxis,aligned_data.left_LFP_table.(['key',num2str(i-1)]),left_sr,gapFillType);
            [left_spect{end+1},left_spect_freq{end+1},left_spect_time{end+1},left_PSD{end+1}]=spectrogram(data,WINDOW,NOVERLAP,NFFT,left_sr);
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
    right_chan_names = cellfun(@(x) aligned_data.DeviceSettings.Right.timeDomainSettings.(x){end},chan_col_names,'UniformOutput',false);
    
    right_spect = {};
    right_spect_freq = {};
    right_spect_time = {};
    right_PSD = {};
    for i = 1:length(right_chan_names)
        same_chan = cellfun(@(x) strcmp(right_chan_names{i},x),right_chan_names(1:i-1));
        if sum(same_chan) == 0  % Not a duplicate channel recording
            [data,time] = addEmptyData(aligned_data.right_taxis,aligned_data.right_LFP_table.(['key',num2str(i-1)]),right_sr,gapFillType);
            [right_spect{end+1},right_spect_freq{end+1},right_spect_time{end+1},right_PSD{end+1}]=spectrogram(data,WINDOW,NOVERLAP,NFFT,right_sr);
        end
    end
end

% Colors for plotting; by rows: 1 = black (LHS), 2 = Maroon (LTO), 3 = Red
% (RHS), 4 = Orange (RTO)
colors = CBMap('GaitEvents',4);

%% Full trial plotting with all gait events
% Left
for j = 1:length(left_chan_names)
    figure;
    hold on;
    pcolor(left_spect_time{j},left_spect_freq{j},10*log10(abs(left_PSD{j})))
    shading interp
    colorbar
    caxis([-110,-30])
    title({'Left';left_chan_names{j};'All gait events'});
    
    a1 = [];
    a2 = [];
    a3 = [];
    a4 = [];
    for k = 1:height(aligned_data.gait_events)
        if ~isnan(aligned_data.gait_events.LHS(k))
            if isempty(a1)
                a1 = xline(aligned_data.gait_events.LHS(k),'linewidth',1.5,'Color',colors(1,:),'DisplayName','LHS');
            else
                xline(aligned_data.gait_events.LHS(k),'linewidth',1.5,'Color',colors(1,:));
            end
        end
        
        if ~isnan(aligned_data.gait_events.LTO(k))
            if isempty(a2)
                a2 = xline(aligned_data.gait_events.LTO(k),'linewidth',1.5,'Color',colors(2,:),'DisplayName','LTO');
            else
                xline(aligned_data.gait_events.LTO(k),'linewidth',1.5,'Color',colors(2,:));
            end
        end
        
        if ~isnan(aligned_data.gait_events.RHS(k))
            if isempty(a3)
                a3 = xline(aligned_data.gait_events.RHS(k),'linewidth',1.5,'Color',colors(3,:),'DisplayName','RHS');
            else
                xline(aligned_data.gait_events.RHS(k),'linewidth',1.5,'Color',colors(3,:));
            end
        end
        
        if ~isnan(aligned_data.gait_events.RTO(k))
            if isempty(a4)
                a4 = xline(aligned_data.gait_events.RTO(k),'linewidth',1.5,'Color',colors(4,:),'DisplayName','RTO');
            else
                xline(aligned_data.gait_events.RTO(k),'linewidth',1.5,'Color',colors(4,:));
            end
        end
    end
    legend([a1,a2,a3,a4]);
    ylim([2.5,50]);
end

% Right
for j = 1:length(right_chan_names)
    figure;
    hold on;
    pcolor(right_spect_time{j},right_spect_freq{j},10*log10(abs(right_PSD{j})))
    shading interp
    colorbar
    caxis([-110,-30])
    title({'Right';right_chan_names{j};"All gait events"});
    
    a1 = [];
    a2 = [];
    a3 = [];
    a4 = [];
    for k = 1:height(aligned_data.gait_events)
        if ~isnan(aligned_data.gait_events.LHS(k))
            if isempty(a1)
                a1 = xline(aligned_data.gait_events.LHS(k),'linewidth',1.5,'Color',colors(1,:),'DisplayName','LHS');
            else
                xline(aligned_data.gait_events.LHS(k),'linewidth',1.5,'Color',colors(1,:));
            end
        end
        
        if ~isnan(aligned_data.gait_events.LTO(k))
            if isempty(a2)
                a2 = xline(aligned_data.gait_events.LTO(k),'linewidth',1.5,'Color',colors(2,:),'DisplayName','LTO');
            else
                xline(aligned_data.gait_events.LTO(k),'linewidth',1.5,'Color',colors(2,:));
            end
        end
        
        if ~isnan(aligned_data.gait_events.RHS(k))
            if isempty(a3)
                a3 = xline(aligned_data.gait_events.RHS(k),'linewidth',1.5,'Color',colors(3,:),'DisplayName','RHS');
            else
                xline(aligned_data.gait_events.RHS(k),'linewidth',1.5,'Color',colors(3,:));
            end
        end
        
        if ~isnan(aligned_data.gait_events.RTO(k))
            if isempty(a4)
                a4 = xline(aligned_data.gait_events.RTO(k),'linewidth',1.5,'Color',colors(4,:),'DisplayName','RTO');
            else
                xline(aligned_data.gait_events.RTO(k),'linewidth',1.5,'Color',colors(4,:));
            end
        end
    end
    legend([a1,a2,a3,a4]);
    ylim([2.5,50]);
end

%% Full trial contralateral events
% Left
for j = 1:length(left_chan_names)
    figure;
    hold on;
    pcolor(left_spect_time{j},left_spect_freq{j},10*log10(abs(left_PSD{j})))
    shading interp
    colorbar
    caxis([-110,-30])
    title({'Left';left_chan_names{j};'Contralateral gait events'});
    
    a3 = [];
    a4 = [];
    for k = 1:height(aligned_data.gait_events)
        if ~isnan(aligned_data.gait_events.RHS(k))
            if isempty(a3)
                a3 = xline(aligned_data.gait_events.RHS(k),'linewidth',1.5,'Color',colors(3,:),'DisplayName','RHS');
            else
                xline(aligned_data.gait_events.RHS(k),'linewidth',1.5,'Color',colors(3,:));
            end
        end
        
        if ~isnan(aligned_data.gait_events.RTO(k))
            if isempty(a4)
                a4 = xline(aligned_data.gait_events.RTO(k),'linewidth',1.5,'Color',colors(4,:),'DisplayName','RTO');
            else
                xline(aligned_data.gait_events.RTO(k),'linewidth',1.5,'Color',colors(4,:));
            end
        end
    end
    legend([a3,a4]);
    ylim([2.5,50]);
end

% Right
for j = 1:length(right_chan_names)
    figure;
    hold on;
    pcolor(right_spect_time{j},right_spect_freq{j},10*log10(abs(right_PSD{j})))
    shading interp
    colorbar
    caxis([-110,-30])
    title({'Right';right_chan_names{j};"Contralateral gait events"});
    
    a1 = [];
    a2 = [];
    for k = 1:height(aligned_data.gait_events)
        if ~isnan(aligned_data.gait_events.LHS(k))
            if isempty(a1)
                a1 = xline(aligned_data.gait_events.LHS(k),'linewidth',1.5,'Color',colors(1,:),'DisplayName','LHS');
            else
                xline(aligned_data.gait_events.LHS(k),'linewidth',1.5,'Color',colors(1,:));
            end
        end
        
        if ~isnan(aligned_data.gait_events.LTO(k))
            if isempty(a2)
                a2 = xline(aligned_data.gait_events.LTO(k),'linewidth',1.5,'Color',colors(2,:),'DisplayName','LTO');
            else
                xline(aligned_data.gait_events.LTO(k),'linewidth',1.5,'Color',colors(2,:));
            end
        end
    end
    legend([a1,a2]);
    ylim([2.5,50]);
end

%% Full trial ipsilateral events
% Left
for j = 1:length(left_chan_names)
    figure;
    hold on;
    pcolor(left_spect_time{j},left_spect_freq{j},10*log10(abs(left_PSD{j})))
    shading interp
    colorbar
    caxis([-110,-30])
    title({'Left';left_chan_names{j};"Ipsilateral gait events"});
    
    a1 = [];
    a2 = [];
    for k = 1:height(aligned_data.gait_events)
        if ~isnan(aligned_data.gait_events.LHS(k))
            if isempty(a1)
                a1 = xline(aligned_data.gait_events.LHS(k),'linewidth',1.5,'Color',colors(1,:),'DisplayName','LHS');
            else
                xline(aligned_data.gait_events.LHS(k),'linewidth',1.5,'Color',colors(1,:));
            end
        end
        
        if ~isnan(aligned_data.gait_events.LTO(k))
            if isempty(a2)
                a2 = xline(aligned_data.gait_events.LTO(k),'linewidth',1.5,'Color',colors(2,:),'DisplayName','LTO');
            else
                xline(aligned_data.gait_events.LTO(k),'linewidth',1.5,'Color',colors(2,:));
            end
        end
    end
    legend([a1,a2]);
    ylim([2.5,50]);
end

% Right
for j = 1:length(right_chan_names)
    figure;
    hold on;
    pcolor(right_spect_time{j},right_spect_freq{j},10*log10(abs(right_PSD{j})))
    shading interp
    colorbar
    caxis([-110,-30])
    title({'Right';right_chan_names{j};"Ipsilateral gait events"});
    
    a3 = [];
    a4 = [];
    for k = 1:height(aligned_data.gait_events)
        if ~isnan(aligned_data.gait_events.RHS(k))
            if isempty(a3)
                a3 = xline(aligned_data.gait_events.RHS(k),'linewidth',1.5,'Color',colors(3,:),'DisplayName','RHS');
            else
                xline(aligned_data.gait_events.RHS(k),'linewidth',1.5,'Color',colors(3,:));
            end
        end
        
        if ~isnan(aligned_data.gait_events.RTO(k))
            if isempty(a4)
                a4 = xline(aligned_data.gait_events.RTO(k),'linewidth',1.5,'Color',colors(4,:),'DisplayName','RTO');
            else
                xline(aligned_data.gait_events.RTO(k),'linewidth',1.5,'Color',colors(4,:));
            end
        end
    end
    legend([a3,a4]);
    ylim([2.5,50]);
end

%% PSD at each gait event
PSD_gait_events.Left = {};
PSD_gait_events.Right = {};

% Left
for j = 1:length(left_chan_names)
    figure;
    for k = 1:length(aligned_data.gait_events.Properties.VariableNames)
        PSD_gait_events.Left{j,k} = zeros(length(left_spect_freq{j}),sum(~isnan(aligned_data.gait_events.(aligned_data.gait_events.Properties.VariableNames{k}))));
        count = 1;
        for m = 1:height(aligned_data.gait_events)
            event_time = aligned_data.gait_events.(aligned_data.gait_events.Properties.VariableNames{k})(m);
            if ~isnan(event_time)
                [~,min_ind] = min(abs(left_spect_time{j}-event_time));
                
                PSD_gait_events.Left{j,k}(:,count) = left_PSD{j}(:,min_ind);
                
                count = count + 1;
            end
        end
        
        subplot(2,2,k);
        plot(repmat(left_spect_freq{j},1,size(PSD_gait_events.Left{j,k},2)),10*log10(abs(PSD_gait_events.Left{j,k})),'Color',[0,0,0,0.2]);
        hold on;
        plot(left_spect_freq{j},mean(10*log10(abs(PSD_gait_events.Left{j,k})),2),'-k','linewidth',1.5);
        xlim([2.5,40])
        xlabel('Frequency (Hz)');
        ylabel('db/Hz');
        title(aligned_data.gait_events.Properties.VariableNames{k});
    end
    sgtitle({subjectID;'Left';left_chan_names{j}});
end

% Right
for j = 1:length(right_chan_names)
    figure;
    for k = 1:length(aligned_data.gait_events.Properties.VariableNames)
        PSD_gait_events.Right{j,k} = zeros(length(right_spect_freq{j}),sum(~isnan(aligned_data.gait_events.(aligned_data.gait_events.Properties.VariableNames{k}))));
        count = 1;
        for m = 1:height(aligned_data.gait_events)
            event_time = aligned_data.gait_events.(aligned_data.gait_events.Properties.VariableNames{k})(m);
            if ~isnan(event_time)
                [~,min_ind] = min(abs(right_spect_time{j}-event_time));
                
                PSD_gait_events.Right{j,k}(:,count) = right_PSD{j}(:,min_ind);
                
                count = count + 1;
            end
        end
        
        subplot(2,2,k);
        plot(repmat(right_spect_freq{j},1,size(PSD_gait_events.Right{j,k},2)),10*log10(abs(PSD_gait_events.Right{j,k})),'Color',[0,0,0,0.2]);
        hold on;
        plot(right_spect_freq{j},mean(10*log10(abs(PSD_gait_events.Right{j,k})),2),'-k','linewidth',1.5);
        xlim([2.5,40])
        xlabel('Frequency (Hz)');
        ylabel('db/Hz');
        title(aligned_data.gait_events.Properties.VariableNames{k});
    end
    sgtitle({subjectID;'Right';right_chan_names{j}});
end

%% band power of +-1 second of event averaged with standard deviation
prePostTime = 1;
[~,band_names] = getBandInd(left_spect_freq{1});

% Left
left_average_power_matrix = cell(length(left_chan_names),length(aligned_data.gait_events.Properties.VariableNames));
left_std_power_matrix = cell(length(left_chan_names),length(aligned_data.gait_events.Properties.VariableNames));
for j = 1:length(left_chan_names)
    band_inds = getBandInd(left_spect_freq{j});
    for k = 1:length(aligned_data.gait_events.Properties.VariableNames)
        vals = [];
        count = 1;
        for m = 1:height(aligned_data.gait_events)
            event_time = aligned_data.gait_events.(aligned_data.gait_events.Properties.VariableNames{k})(m);
            if ~isnan(event_time)
                [~,min_ind_pre] = min(abs(left_spect_time{j}-(event_time-prePostTime)));
                [~,min_ind_post] = min(abs(left_spect_time{j}-(event_time+prePostTime)));
                
                for n = 1:length(band_names)
                    vals(count,:,n) = mean(10*log10(abs(left_PSD{j}(band_inds(n,1):band_inds(n,2),min_ind_pre:min_ind_post))));
                end
                count = count + 1;
            end
        end
        left_average_power_matrix{j,k} = mean(vals,1);
        left_std_power_matrix{j,k} = std(vals,0,1);
    end
end

% Right
right_average_power_matrix = cell(length(right_chan_names),length(aligned_data.gait_events.Properties.VariableNames));
right_std_power_matrix = cell(length(right_chan_names),length(aligned_data.gait_events.Properties.VariableNames));
for j = 1:length(right_chan_names)
    band_inds = getBandInd(right_spect_freq{j});
    for k = 1:length(aligned_data.gait_events.Properties.VariableNames)
        vals = [];
        count = 1;
        for m = 1:height(aligned_data.gait_events)
            event_time = aligned_data.gait_events.(aligned_data.gait_events.Properties.VariableNames{k})(m);
            if ~isnan(event_time)
                [~,min_ind_pre] = min(abs(right_spect_time{j}-(event_time-prePostTime)));
                [~,min_ind_post] = min(abs(right_spect_time{j}-(event_time+prePostTime)));
                
                for n = 1:length(band_names)
                    vals(count,:,n) = mean(10*log10(abs(right_PSD{j}(band_inds(n,1):band_inds(n,2),min_ind_pre:min_ind_post))));
                end
                count = count + 1;
            end
        end
        right_average_power_matrix{j,k} = mean(vals,1);
        right_std_power_matrix{j,k} = std(vals,0,1);
    end
end

end

function [band_indicies,band_names] = getBandInd(freq_vec)
band_indicies = zeros(8,2);
band_names = {'Delta','Theta','Alpha','Beta','Low Beta','High Beta','Gamma','High Gamma'};

% delta
[~,band_indicies(1,1)] = min(abs(freq_vec-0.5));
[~,band_indicies(1,2)] = min(abs(freq_vec-3));

% theta
[~,band_indicies(2,1)] = min(abs(freq_vec-4));
[~,band_indicies(2,2)] = min(abs(freq_vec-7));

% alpha
[~,band_indicies(3,1)] = min(abs(freq_vec-8));
[~,band_indicies(3,2)] = min(abs(freq_vec-12));

% beta
[~,band_indicies(4,1)] = min(abs(freq_vec-13));
[~,band_indicies(4,2)] = min(abs(freq_vec-30));

% low beta
[~,band_indicies(5,1)] = min(abs(freq_vec-13));
[~,band_indicies(5,2)] = min(abs(freq_vec-20));

% high beta
[~,band_indicies(6,1)] = min(abs(freq_vec-21));
[~,band_indicies(6,2)] = min(abs(freq_vec-30));

% gamma
[~,band_indicies(7,1)] = min(abs(freq_vec-31));
[~,band_indicies(7,2)] = min(abs(freq_vec-60));

% high gamma
[~,band_indicies(8,1)] = min(abs(freq_vec-61));
[~,band_indicies(8,2)] = min(abs(freq_vec-100));

end