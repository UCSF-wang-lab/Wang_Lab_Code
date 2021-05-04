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

% searchForGaitBiomarkers(aligned_data,A,[0,50],{{'LHS','RTO','RHS','LTO'}},'RCS03',0)

%% Extract data
tic
if isfield(signal_analysis_data,'Left')
    left_sr = uniquetol(aligned_data.DeviceSettings.Left.timeDomainSettings.samplingRate,1);
    time_res_left = uniquetol(diff(signal_analysis_data.Left.Time{1}),1);
    if isfield(signal_analysis_data.Left,'PSD')
        analysis_type = 'FT';
    else
        analysis_type = 'CWT';
    end
    
    event_pairs = nchoosek(aligned_data.gait_events.Properties.VariableNames,2);
    [freq_bin_inds.Left,freq_vals.Left] = genFreqBinPairs(signal_analysis_data.Left.Freq_Values{1},freq_lim);
    
    channel_power.Left = cell(1,length(signal_analysis_data.Left.Chan_Names));
    channel_anova_matrix.Left = cell(1,length(signal_analysis_data.Left.Chan_Names));
    channel_mult_compare_matrix.Left = cell(1,size(event_pairs,1));
    for i = 1:length(signal_analysis_data.Left.Chan_Names)
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
                        for m = 1:size(freq_bin_inds.Left,1)
                            vals_power.(aligned_data.gait_events.Properties.VariableNames{j})(count,m) = mean(temp(freq_bin_inds.Left(m,1):freq_bin_inds.Left(m,2)));
                            vals_phase.(aligned_data.gait_events.Properties.VariableNames{j})(count,m) = mean(temp1(freq_bin_inds.Left(m,1):freq_bin_inds.Left(m,2)));
                        end
                        count = count + 1;
                    end
                end
            end
        end
        channel_power.Left{i} = vals_power;
        
        channel_anova_matrix.Left{i} = nan(sum(diff(freq_bin_inds.Left,1,2)==0),sum(diff(freq_bin_inds.Left,1,2)==0));
        for n = 1:size(freq_bin_inds.Left,1)
            X = [];
            groups = {};
            for o = 1:length(fields(vals_power))
                vals = vals_power.(aligned_data.gait_events.Properties.VariableNames{o})(:,n);
                X = [X;vals];
                groups = [groups;repelem(aligned_data.gait_events.Properties.VariableNames(o),length(vals),1)];
            end
            [p,~,stats] = anova1(X,groups,'off');
            [c,~,~,gnames] = multcompare(stats,'Display','off');
            channel_anova_matrix.Left{i}(freq_bin_inds.Left(n,1),freq_bin_inds.Left(n,2)) = p;
            for r = 1:size(c,1)
                mat_ind = find(and(strcmp(gnames(c(r,1)),event_pairs(:,1)),strcmp(gnames(c(r,2)),event_pairs(:,2))));
                if isempty(channel_mult_compare_matrix.Left{mat_ind})
                    channel_mult_compare_matrix.Left{mat_ind} = cellmat(6,1,sum(diff(freq_bin_inds.Left,1,2)==0),sum(diff(freq_bin_inds.Left,1,2)==0),nan);
                end
                channel_mult_compare_matrix.Left{i}{mat_ind}(freq_bin_inds.Left(n,1),freq_bin_inds.Left(n,2)) = c(r,6);
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

    event_pairs = nchoosek(aligned_data.gait_events.Properties.VariableNames,2);
    [freq_bin_inds.Right,freq_vals.Right] = genFreqBinPairs(signal_analysis_data.Right.Freq_Values{1},freq_lim);

    channel_power.Right = cell(1,length(signal_analysis_data.Right.Chan_Names));
    channel_anova_matrix.Right = cell(1,length(signal_analysis_data.Right.Chan_Names));
    channel_mult_compare_matrix.Right = cell(1,size(event_pairs,1));
    for i = 1:length(signal_analysis_data.Right.Chan_Names)
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
                    [~,event_ind] = min(abs(signal_analysis_data.Right.Time{i}-event_time));
                    if isfield(signal_analysis_data.Right,'PSD')
                        temp = 20*log10(abs(signal_analysis_data.Right.Values{i}(:,event_ind)));
                        temp1 = angle(signal_analysis_data.Right.Values{i}(:,event_ind));
                    else
                        temp = abs(signal_analysis_data.Right.Values{i}(:,event_ind));
                        temp1 = angle(signal_analysis_data.Right.Values{i}(:,event_ind));
                    end

                    if sum(isinf(temp),'all') == 0
                        for m = 1:size(freq_bin_inds.Right,1)
                            vals_power.(aligned_data.gait_events.Properties.VariableNames{j})(count,m) = mean(temp(freq_bin_inds.Right(m,1):freq_bin_inds.Right(m,2)));
                            vals_phase.(aligned_data.gait_events.Properties.VariableNames{j})(count,m) = mean(temp1(freq_bin_inds.Right(m,1):freq_bin_inds.Right(m,2)));
                        end
                        count = count + 1;
                    end
                end
            end
        end
        channel_power.Right{i} = vals_power;

        channel_anova_matrix.Right{i} = nan(sum(diff(freq_bin_inds.Right,1,2)==0),sum(diff(freq_bin_inds.Right,1,2)==0));
        for n = 1:size(freq_bin_inds.Right,1)
            X = [];
            groups = {};
            for o = 1:length(fields(vals_power))
                vals = vals_power.(aligned_data.gait_events.Properties.VariableNames{o})(:,n);
                X = [X;vals];
                groups = [groups;repelem(aligned_data.gait_events.Properties.VariableNames(o),length(vals),1)];
            end
            [p,~,stats] = anova1(X,groups,'off');
            [c,~,~,gnames] = multcompare(stats,'Display','off');
            channel_anova_matrix.Right{i}(freq_bin_inds.Right(n,1),freq_bin_inds.Right(n,2)) = p;
            for r = 1:size(c,1)
                mat_ind = find(and(strcmp(gnames(c(r,1)),event_pairs(:,1)),strcmp(gnames(c(r,2)),event_pairs(:,2))));
                if isempty(channel_mult_compare_matrix.Right{mat_ind})
                    channel_mult_compare_matrix.Right{mat_ind} = cellmat(6,1,sum(diff(freq_bin_inds.Left,1,2)==0),sum(diff(freq_bin_inds.Left,1,2)==0),nan);
                end
                channel_mult_compare_matrix.Right{i}{mat_ind}(freq_bin_inds.Left(n,1),freq_bin_inds.Left(n,2)) = c(r,6);
            end
        end
    end
end

%% Plot to visualize p-vals
fig_vec = [];
if isfield(channel_anova_matrix,'Left')
    for i = 1:length(channel_anova_matrix.Left)
        fig_vec(end+1) = figure;
        ax_hand = pcolor(1:size(channel_anova_matrix.Left{i},1),1:size(channel_anova_matrix.Left{i},2),channel_anova_matrix.Left{i});
        colormap(flipud(jet));
        caxis([0,1]);
        set(gca,'YDir',"reverse");
        title({subjectID;'Left';signal_analysis_data.Left.Chan_Names{i};'ANOVA'});
        xlabel('End bin');
        ylabel('Start bin');
        c_hand = colorbar;
        c_hand.Title.String = 'p-val';
    end
end

if isfield(channel_anova_matrix,'Right')
    for i = 1:length(channel_anova_matrix.Right)
        fig_vec(end+1) = figure;
        pcolor(1:size(channel_anova_matrix.Right{i},1),1:size(channel_anova_matrix.Right{i},2),channel_anova_matrix.Right{i});
        colormap(flipud(jet));
        caxis([0,1]);
        set(gca,'YDir',"reverse");
        title({subjectID;'Right';signal_analysis_data.Right.Chan_Names{i};'ANOVA'});
        xlabel('End bin');
        ylabel('Start bin');
        c_hand = colorbar;
        c_hand.Title.String = 'p-val';
    end
end

if isfield(channel_mult_compare_matrix,'Left')
    for i = 1:length(signal_analysis_data.Left.Chan_Names)
        for j = 1:size(event_pairs,1)
            fig_vec(end+1) = figure;
            pcolor(1:size(channel_mult_compare_matrix.Left{i}{j},1),1:size(channel_mult_compare_matrix.Left{i}{j},2),channel_mult_compare_matrix.Left{i}{j});
            colormap(flipud(jet));
            caxis([0,1]);
            set(gca,'YDir',"reverse");
            title({subjectID;'Left';signal_analysis_data.Left.Chan_Names{i};[event_pairs{j,1},' vs. ', event_pairs{j,2},' Multiple Comparison']});
            xlabel('End bin');
            ylabel('Start bin');
            c_hand = colorbar;
            c_hand.Title.String = 'p-val';
        end
    end
end

if isfield(channel_mult_compare_matrix,'Right')
    for i = 1:length(signal_analysis_data.Right.Chan_Names)
        for j = 1:size(event_pairs,1)
            fig_vec(end+1) = figure;
            pcolor(1:size(channel_mult_compare_matrix.Right{i}{j},1),1:size(channel_mult_compare_matrix.Right{i}{j},2),channel_mult_compare_matrix.Right{i}{j});
            colormap(flipud(jet));
            caxis([0,1]);
            set(gca,'YDir',"reverse");
            title({subjectID;'Right';signal_analysis_data.Right.Chan_Names{i};[event_pairs{j,1},' vs. ', event_pairs{j,2},' Multiple Comparison']});
            xlabel('End bin');
            ylabel('Start bin');
            c_hand = colorbar;
            c_hand.Title.String = 'p-val';
        end
    end
end

%% Save plots
if save_flag
    save_dir = uigetdir();
    
%     figure_format(12,8,12);
    
    % check if saving folders exist
    if ~isfolder(fullfile(save_dir,'GaitBiomarkerSearch'))
        mkdir(fullfile(save_dir,'GaitBiomarkerSearch'));
    end
    
    if ~isfolder(fullfile(save_dir,'GaitBiomarkerSearch',aligned_data.stim_condition))
        mkdir(fullfile(save_dir,'GaitBiomarkerSearch',aligned_data.stim_condition))
    end
    
    if strcmp(analysis_type,'FT')
        if ~isfolder(fullfile(save_dir,'GaitBiomarkerSearch',aligned_data.stim_condition,'FT'))
            mkdir(fullfile(save_dir,'GaitBiomarkerSearch',aligned_data.stim_condition,'FT'))
        end
    elseif strcmp(analysis_type,'CWT')
        if ~isfolder(fullfile(save_dir,'GaitBiomarkerSearch',aligned_data.stim_condition,'CWT'))
            mkdir(fullfile(save_dir,'GaitBiomarkerSearch',aligned_data.stim_condition,'CWT'))
        end
    end
    
    folders_to_check = {'FIG_files','PDF_files','TIFF_files'};
    extension = {'.fig','.pdf','.tiff'};
    for n = 1:length(folders_to_check)
        if strcmp(analysis_type,'FT')
            if ~isfolder(fullfile(save_dir,'GaitBiomarkerSearch',aligned_data.stim_condition,'FT',folders_to_check{n}))
                mkdir(fullfile(save_dir,'GaitBiomarkerSearch',aligned_data.stim_condition,'FT',folders_to_check{n}));
            end
        elseif strcmp(analysis_type,'CWT')
            if ~isfolder(fullfile(save_dir,'GaitBiomarkerSearch',aligned_data.stim_condition,'CWT',folders_to_check{n}))
                mkdir(fullfile(save_dir,'GaitBiomarkerSearch',aligned_data.stim_condition,'CWT',folders_to_check{n}));
            end
        end
    end
    
    for i = 1:length(fig_vec)
        curr_axes = gca(fig_vec(i));
        save_name = [];
        for j = 1:length(curr_axes.Parent.Children(2).Title.String)
            if isempty(save_name)
                save_name = curr_axes.Parent.Children(2).Title.String{j};
            else
                save_name = [save_name,' ', curr_axes.Parent.Children(2).Title.String{j}];
            end
        end
        
        if strcmp(analysis_type,'FT')
            savefig(fig_vec(i),fullfile(save_dir,'GaitBiomarkerSearch',aligned_data.stim_condition,'FT',folders_to_check{1},strrep(strrep(save_name,' ','_'),'.','')));
        elseif strcmp(analysis_type,'CWT')
            savefig(fig_vec(i),fullfile(save_dir,'GaitBiomarkerSearch',aligned_data.stim_condition,'CWT',folders_to_check{1},strrep(strrep(save_name,' ','_'),'.','')));
        end
        
        for k = 2:length(folders_to_check)
            if strcmp(analysis_type,'FT')
                print(fig_vec(i),[fullfile(save_dir,'GaitBiomarkerSearch',aligned_data.stim_condition,'FT',folders_to_check{k},strrep(strrep(save_name,' ','_'),'.','')),extension{k}],'-r300',['-d',extension{k}(2:end)]);
            elseif strcmp(analysis_type,'CWT')
                print(fig_vec(i),[fullfile(save_dir,'GaitBiomarkerSearch',aligned_data.stim_condition,'CWT',folders_to_check{k},strrep(strrep(save_name,' ','_'),'.','')),extension{k}],'-r300',['-d',extension{k}(2:end)]);
            end
        end
    end
end

toc
varargout = {channel_anova_matrix,channel_mult_compare_matrix,event_pairs,freq_vals};
end

function [freq_bin_inds,freqs] = genFreqBinPairs(freq_vals,freq_lim)
ind_start = find(freq_vals >= freq_lim(1),1,'first');
ind_end = find(freq_vals <= freq_lim(2),1,'last');
n = ind_end - ind_start + 1;
skips = 1:n-1;
tot_pairs = (n^2-n)/2;
freqs = freq_vals(ind_start:ind_end);
freq_bin_inds = nan(tot_pairs,2);

% Crazy indexing. Look up triangle numbers to recalculate these equations
f = @(x) n*(x+1)-(x^2+x)/2;
freq_bin_inds(1:n,:) = repmat([1:n]',1,2);
for i = 1:length(skips)
    freq_bin_inds(f(i-1)+1:f(i),1) = transpose(1:n-skips(i));
    freq_bin_inds(f(i-1)+1:f(i),2) = freq_bin_inds(f(i-1)+1:f(i),1) + skips(i);
end
end