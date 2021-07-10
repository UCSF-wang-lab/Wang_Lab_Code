function genPSDGaitEvents(aligned_data,signal_analysis_data,subjectID,save_flag)
if ~exist('save_flag','var')
    save_flag = 0;
end

if ~exist('subjectID','var')
    subjectID = 'RCSXX';
end

%% PSD at each gait event
% Left
if isfield(signal_analysis_data,'Left')
    PSD_gait_events.Left = {};
    for i = 1:length(signal_analysis_data.Left.Chan_Names)
        for j = 1:length(aligned_data.gait_events.Properties.VariableNames)
            if ~isfield(PSD_gait_events.Left,aligned_data.gait_events.Properties.VariableNames{j})
                PSD_gait_events.Left.(aligned_data.gait_events.Properties.VariableNames{j}) = cell(1,length(signal_analysis_data.Left.Chan_Names));
            end
            for k = 1:height(aligned_data.gait_events)
                event_time = aligned_data.gait_events.(aligned_data.gait_events.Properties.VariableNames{j})(k);
                if ~isnan(event_time)
                    [~,min_ind] = min(abs(signal_analysis_data.Left.Time{i}-event_time));
                    power_values = signal_analysis_data.Left.PSD{i}(:,min_ind);
                    if sum(power_values == 0) == 0
                        if isempty(PSD_gait_events.Left.(aligned_data.gait_events.Properties.VariableNames{j}))
                            PSD_gait_events.Left.(aligned_data.gait_events.Properties.VariableNames{j}){i} = power_values;
                        else
                            PSD_gait_events.Left.(aligned_data.gait_events.Properties.VariableNames{j}){i}(:,end+1) = power_values;
                        end
                    end
                end
            end
        end
    end
end

% Right
if isfield(signal_analysis_data,'Right')
    PSD_gait_events.Right = {};
    for i = 1:length(signal_analysis_data.Right.Chan_Names)
        for j = 1:length(aligned_data.gait_events.Properties.VariableNames)
            if ~isfield(PSD_gait_events.Right,aligned_data.gait_events.Properties.VariableNames{j})
                PSD_gait_events.Right.(aligned_data.gait_events.Properties.VariableNames{j}) = cell(1,length(signal_analysis_data.Right.Chan_Names));
            end
            for k = 1:height(aligned_data.gait_events)
                event_time = aligned_data.gait_events.(aligned_data.gait_events.Properties.VariableNames{j})(k);
                if ~isnan(event_time)
                    [~,min_ind] = min(abs(signal_analysis_data.Right.Time{i}-event_time));
                    power_values = signal_analysis_data.Right.PSD{i}(:,min_ind);
                    if sum(power_values == 0) == 0
                        if isempty(PSD_gait_events.Right.(aligned_data.gait_events.Properties.VariableNames{j}))
                            PSD_gait_events.Right.(aligned_data.gait_events.Properties.VariableNames{j}){i} = power_values;
                        else
                            PSD_gait_events.Right.(aligned_data.gait_events.Properties.VariableNames{j}){i}(:,end+1) = power_values;
                        end
                    end
                end
            end
        end
    end
end

%% Plot
fig_vec = [];
if isfield(PSD_gait_events,'Left')
    for i = 1:length(signal_analysis_data.Left.Chan_Names)
        fig_vec(end+1) = figure;
        for j = 1:length(aligned_data.gait_events.Properties.VariableNames)
            curr_event = aligned_data.gait_events.Properties.VariableNames{j};
            freq_mat = repmat(signal_analysis_data.Left.Freq_Values{i},1,size(PSD_gait_events.Left.(curr_event){i},2));
            avg_power = mean(10*log10(abs(PSD_gait_events.Left.(curr_event){i})),2);
            
            switch curr_event
                case 'LTO'
                    subplot(2,2,1);
                case 'LHS'
                    subplot(2,2,2);
                case 'RTO'
                    subplot(2,2,3);
                case 'RHS'
                    subplot(2,2,4);
            end
            
            plot(freq_mat,10*log10(abs(PSD_gait_events.Left.(curr_event){i})),'Color',[0,0,0,0.2]);
            hold on;
            plot(freq_mat(:,1),avg_power,'-k','linewidth',1.5);
            hold off;
            title(aligned_data.gait_events.Properties.VariableNames{j});
            xlim([2.5,40])
            xlabel('Frequency (Hz)');
            ylabel('db/Hz');
        end
        sgtitle({[subjectID,' Left'];signal_analysis_data.Left.Chan_Names{i}});
    end
end

% Right
if isfield(PSD_gait_events,'Right')
    for i = 1:length(signal_analysis_data.Right.Chan_Names)
        fig_vec(end+1) = figure;
        for j = 1:length(aligned_data.gait_events.Properties.VariableNames)
            curr_event = aligned_data.gait_events.Properties.VariableNames{j};
            freq_mat = repmat(signal_analysis_data.Right.Freq_Values{i},1,size(PSD_gait_events.Right.(curr_event){i},2));
            avg_power = mean(10*log10(abs(PSD_gait_events.Right.(curr_event){i})),2);
            
            switch curr_event
                case 'LTO'
                    subplot(2,2,1);
                case 'LHS'
                    subplot(2,2,2);
                case 'RTO'
                    subplot(2,2,3);
                case 'RHS'
                    subplot(2,2,4);
            end
            
            plot(freq_mat,10*log10(abs(PSD_gait_events.Right.(curr_event){i})),'Color',[0,0,0,0.2]);
            hold on;
            plot(freq_mat(:,1),avg_power,'-k','linewidth',1.5);
            hold off;
            title(aligned_data.gait_events.Properties.VariableNames{j});
            xlim([2.5,40])
            xlabel('Frequency (Hz)');
            ylabel('db/Hz');
        end
        sgtitle({[subjectID,' Right'];signal_analysis_data.Right.Chan_Names{i}});
    end
end

if save_flag
    save_dir = uigetdir();
    
    figure_format(6,6,10);
    
    % check if saving folders exist
    if ~isfolder(fullfile(save_dir,'PSD'))
        mkdir(fullfile(save_dir,'PSD'));
    end
    
    if ~isfolder(fullfile(save_dir,'PSD',[aligned_data.stim_condition,'_STIM']))
        mkdir(fullfile(save_dir,'PSD',[aligned_data.stim_condition,'_STIM']))
    end
    
    if ~isfolder(fullfile(save_dir,'PSD',[aligned_data.stim_condition,'_STIM'],'FT'))
        mkdir(fullfile(save_dir,'PSD',[aligned_data.stim_condition,'_STIM'],'FT'))
    end
    
    folders_to_check = {'FIG_files','PDF_files','TIFF_files'};
    extension = {'.fig','.pdf','.tiff'};
    for n = 1:length(folders_to_check)
        if ~isfolder(fullfile(save_dir,'PSD',[aligned_data.stim_condition,'_STIM'],'FT',folders_to_check{n}))
            mkdir(fullfile(save_dir,'PSD',[aligned_data.stim_condition,'_STIM'],'FT',folders_to_check{n}));
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
        
        savefig(fig_vec(i),fullfile(save_dir,'PSD',[aligned_data.stim_condition,'_STIM'],'FT',folders_to_check{1},strrep(save_name,' ','_')));
        for k = 2:length(folders_to_check)
            print(fig_vec(i),[fullfile(save_dir,'PSD',[aligned_data.stim_condition,'_STIM'],'FT',folders_to_check{k},strrep(save_name,' ','_')),extension{k}],'-r300',['-d',extension{k}(2:end)]);
        end
        
    end
end
end