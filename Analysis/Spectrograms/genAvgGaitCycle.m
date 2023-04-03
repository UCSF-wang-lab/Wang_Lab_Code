function genAvgGaitCycle(aligned_data,signal_analysis_data,varargin)
% cycle_start_event,n_percent_bins,baseline_data,subjectID,save_flag)

for i = 1:2:nargin-2
    switch varargin{i}
        case 'cycle_start_event'
            cycle_start_event = varargin{i+1};
        case 'n_percent_bins'
            n_percent_bins = varargin{i+1};
        case 'baseline_data'
            baseline_data = varargin{i+1};
        case 'baseline_normalization'
            baseline_normalization = varargin{i+1};
        case 'subjectID'
            subjectID = varargin{i+1};
        case 'save_flag'
            save_flag = varargin{i+1};
    end
end

if ~exist('cycle_start_event','var')
    cycle_start_event = 'LHS';
end

if ~exist('n_percent_bins','var')
    n_percent_bins = 100;
end

if ~exist('baseline_data','var')
    baseline_data.Left = [];
    baseline_data.Right = [];
end

if ~exist('baseline_normalization','var')
    baseline_normalization = 'none'; % none|percent_change|subtraction
end

if ~exist('subjectID','var')
    subjectID = 'RCSXX';
end

if ~exist('save_flag','var')
    save_flag = 0;
end

%% Average gait cycle value
% re-sort gait events with the starting event
gait_events_sorted = sortGaitEvents(aligned_data.gait_events,cycle_start_event);

% Go through all valid gait cycles
% Left
if isfield(signal_analysis_data,'Left')
    if isfield(signal_analysis_data.Left,'PSD')
        analysis_type = 'FT';
    else
        analysis_type = 'CWT';
    end
    gait_cycle_avg.Left = cell(1,length(signal_analysis_data.Left.Chan_Names));
    for i = 1:length(signal_analysis_data.Left.Chan_Names)
        gait_cycle_mat_left = zeros(length(signal_analysis_data.Left.Freq_Values{i}),n_percent_bins,1);
        count = 1;
        for j = 1:height(gait_events_sorted)-1
            if ~isnan(gait_events_sorted.(cycle_start_event)(j)) && ~isnan(gait_events_sorted.(cycle_start_event)(j+1)) && (diff(gait_events_sorted.(cycle_start_event)(j:j+1)) < 2)
                [~,start_ind] = min(abs(signal_analysis_data.Left.Time{i}-gait_events_sorted.(cycle_start_event)(j)));
                [~,end_ind] = min(abs(signal_analysis_data.Left.Time{i}-gait_events_sorted.(cycle_start_event)(j+1)));
                
                if isfield(signal_analysis_data.Left,'PSD')
                    data_snip = 10*log10(abs(signal_analysis_data.Left.PSD{i}(:,start_ind:end_ind)));
                else
                    data_snip = abs(signal_analysis_data.Left.Values{i}(:,start_ind:end_ind));
                end
                
                if sum(isinf(data_snip),'all') == 0
                    percent_inds = round(linspace(1,size(data_snip,2),n_percent_bins+1));
                    for k = 1:length(percent_inds)-1
                        if k == 1
                            gait_cycle_mat_left(:,k,count) = mean(data_snip(:,percent_inds(k):percent_inds(k+1)),2);
                        else
                            gait_cycle_mat_left(:,k,count) = mean(data_snip(:,percent_inds(k)+1:percent_inds(k+1)),2);
                        end
                    end
                    count = count + 1;
                end
            end
        end
        gait_cycle_avg.Left{i} = mean(gait_cycle_mat_left,3);
    end
end

% Right
if isfield(signal_analysis_data,'Right')
    if isfield(signal_analysis_data.Right,'PSD')
        analysis_type = 'FT';
    else
        analysis_type = 'CWT';
    end
    gait_cycle_avg.Right = cell(1,length(signal_analysis_data.Right.Chan_Names));
    for i = 1:length(signal_analysis_data.Right.Chan_Names)
        gait_cycle_mat_right = zeros(length(signal_analysis_data.Right.Freq_Values{i}),n_percent_bins,1);
        count = 1;
        for j = 1:height(gait_events_sorted)-1
            if ~isnan(gait_events_sorted.(cycle_start_event)(j)) && ~isnan(gait_events_sorted.(cycle_start_event)(j+1)) && (diff(gait_events_sorted.(cycle_start_event)(j:j+1)) < 2)
                [~,start_ind] = min(abs(signal_analysis_data.Right.Time{i}-gait_events_sorted.(cycle_start_event)(j)));
                [~,end_ind] = min(abs(signal_analysis_data.Right.Time{i}-gait_events_sorted.(cycle_start_event)(j+1)));
                if isfield(signal_analysis_data.Right,'PSD')
                    data_snip = 10*log10(abs(signal_analysis_data.Right.PSD{i}(:,start_ind:end_ind)));
                else
                    data_snip = abs(signal_analysis_data.Right.Values{i}(:,start_ind:end_ind));
                end
                if sum(isinf(data_snip),'all') == 0
                    percent_inds = round(linspace(1,size(data_snip,2),n_percent_bins+1));
                    for k = 1:length(percent_inds)-1
                        if k == 1
                            gait_cycle_mat_right(:,k,count) = mean(data_snip(:,percent_inds(k):percent_inds(k+1)),2);
                        else
                            gait_cycle_mat_right(:,k,count) = mean(data_snip(:,percent_inds(k)+1:percent_inds(k+1)),2);
                        end
                    end
                    count = count + 1;
                end
            end
        end
        gait_cycle_avg.Right{i} = mean(gait_cycle_mat_right,3);
    end
end

%% Plot
fig_vec = [];
if isfield(gait_cycle_avg,'Left')
    gait_cycle_avg.Left = normalizeWithBaseline(gait_cycle_avg.Left,baseline_data.Left,baseline_normalization);
    for i = 1:length(signal_analysis_data.Left.Chan_Names)
        fig_vec(end+1) = figure;
        if isfield(signal_analysis_data.Left,'PSD')
            pcolor(1:100,signal_analysis_data.Left.Freq_Values{i},gait_cycle_avg.Left{i});
            ylim([2.5,50]);
        else
            ax = pcolor(1:100,log2(signal_analysis_data.Left.Freq_Values{i}),gait_cycle_avg.Left{i});
            ticks = logspace(log10(2.5),log10(50),10);
            ax.Parent.YTick = log2(ticks);
            ax.Parent.YTickLabel = ticks;
            ylim([log2(2.5),log2(50)]);
        end
        shading interp;
        title({[subjectID,' Left'];signal_analysis_data.Left.Chan_Names{i}});
        xlabel('Gait Cycle %');
        ylabel('Frequency (Hz)');
    end
end

if isfield(gait_cycle_avg,'Right')
    gait_cycle_avg.Right = normalizeWithBaseline(gait_cycle_avg.Right,baseline_data.Right,baseline_normalization);
    for i = 1:length(signal_analysis_data.Right.Chan_Names)
        fig_vec(end+1) = figure;
        if isfield(signal_analysis_data.Right,'PSD')
            pcolor(1:100,signal_analysis_data.Right.Freq_Values{i},gait_cycle_avg.Right{i});
            ylim([2.5,50]);
        else
            ax = pcolor(1:100,log2(signal_analysis_data.Right.Freq_Values{i}),gait_cycle_avg.Right{i});
            ticks = logspace(log10(2.5),log10(50),10);
            ax.Parent.YTick = log2(ticks);
            ax.Parent.YTickLabel = ticks;
            ylim([log2(2.5),log2(50)]);
        end
        shading interp;
        title({[subjectID,' Right'];signal_analysis_data.Right.Chan_Names{i}});
        xlabel('Gait Cycle %');
        ylabel('Frequency (Hz)');
    end
end

%% Save plots
if save_flag
    save_dir = uigetdir();
    
    figure_format(6,6,10);
    
    % check if saving folders exist
    if ~isfolder(fullfile(save_dir,'AvgGaitCycleSpec'))
        mkdir(fullfile(save_dir,'AvgGaitCycleSpec'));
    end
    
    if ~isfolder(fullfile(save_dir,'AvgGaitCycleSpec',[aligned_data.stim_condition,'_STIM']))
        mkdir(fullfile(save_dir,'AvgGaitCycleSpec',[aligned_data.stim_condition,'_STIM']))
    end
    
    if strcmp(analysis_type,'FT')
        if ~isfolder(fullfile(save_dir,'AvgGaitCycleSpec',[aligned_data.stim_condition,'_STIM'],'FT'))
            mkdir(fullfile(save_dir,'AvgGaitCycleSpec',[aligned_data.stim_condition,'_STIM'],'FT'))
        end
    elseif strcmp(analysis_type,'CWT')
        if ~isfolder(fullfile(save_dir,'AvgGaitCycleSpec',[aligned_data.stim_condition,'_STIM'],'CWT'))
            mkdir(fullfile(save_dir,'AvgGaitCycleSpec',[aligned_data.stim_condition,'_STIM'],'CWT'))
        end
    end
    
    folders_to_check = {'FIG_files','PDF_files','TIFF_files'};
    extension = {'.fig','.pdf','.tiff'};
    for n = 1:length(folders_to_check)
        if strcmp(analysis_type,'FT')
            if ~isfolder(fullfile(save_dir,'AvgGaitCycleSpec',[aligned_data.stim_condition,'_STIM'],'FT',folders_to_check{n}))
                mkdir(fullfile(save_dir,'AvgGaitCycleSpec',[aligned_data.stim_condition,'_STIM'],'FT',folders_to_check{n}));
            end
        elseif strcmp(analysis_type,'CWT')
            if ~isfolder(fullfile(save_dir,'AvgGaitCycleSpec',[aligned_data.stim_condition,'_STIM'],'CWT',folders_to_check{n}))
                mkdir(fullfile(save_dir,'AvgGaitCycleSpec',[aligned_data.stim_condition,'_STIM'],'CWT',folders_to_check{n}));
            end
        end
    end
    
    for i = 1:length(fig_vec)
        curr_axes = gca(fig_vec(i));
        save_name = [];
        for j = 1:length(curr_axes.Title.String)
            if isempty(save_name)
                save_name = curr_axes.Title.String{j};
            else
                save_name = [save_name,' ', curr_axes.Title.String{j}];
            end
        end
        
        if strcmp(analysis_type,'FT')
            savefig(fig_vec(i),fullfile(save_dir,'AvgGaitCycleSpec',[aligned_data.stim_condition,'_STIM'],'FT',folders_to_check{1},strrep(save_name,' ','_')));
        elseif strcmp(analysis_type,'CWT')
            savefig(fig_vec(i),fullfile(save_dir,'AvgGaitCycleSpec',[aligned_data.stim_condition,'_STIM'],'CWT',folders_to_check{1},strrep(save_name,' ','_')));
        end
        
        for k = 2:length(folders_to_check)
            if strcmp(analysis_type,'FT')
                print(fig_vec(i),[fullfile(save_dir,'AvgGaitCycleSpec',[aligned_data.stim_condition,'_STIM'],'FT',folders_to_check{k},strrep(save_name,' ','_')),extension{k}],'-r300',['-d',extension{k}(2:end)]);
            elseif strcmp(analysis_type,'CWT')
                print(fig_vec(i),[fullfile(save_dir,'AvgGaitCycleSpec',[aligned_data.stim_condition,'_STIM'],'CWT',folders_to_check{k},strrep(save_name,' ','_')),extension{k}],'-r300',['-d',extension{k}(2:end)]);
            end
        end
    end
end
end

function normalized_data = normalizeWithBaseline(gait_cycle_avg,baseline_data,normalize_type)
normalized_data = gait_cycle_avg;
if ~strcmp(normalize_type,'none')
    baseline_vals = cellfun(@(x) mean(10*log10(abs(x)),2),baseline_data.Values,'UniformOutput',false);
    assert(length(baseline_data.Values)==length(gait_cycle_avg));
    
    switch normalize_type
        case 'percent_change'
            normalized_data = cellfun(@(x,y)(x-y)./abs(y),gait_cycle_avg,baseline_vals,'UniformOutput',false);
        case 'subtraction'
            normalized_data = cellfun(@(x,y)x-y,gait_cycle_avg,baseline_vals,'UniformOutput',false);
    end
end
end