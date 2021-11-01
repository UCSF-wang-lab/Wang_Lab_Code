function varargout = searchForGaitBiomarkers(aligned_data,signal_analysis_data,freq_lim,event_compare,subjectID,save_dir,save_type,save_flag)
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

if ~exist('save_dir','var')
    save_dir = [];
end

if ~exist('save_type','var') || isempty(save_type)
    save_type = 'all'; % all|anova|mult_compare|summed|summed_all|summed_exclude|stat_table
end

if ~exist('save_flag','var') || isempty(save_flag)
    save_flag = 0;
end

% searchForGaitBiomarkers(aligned_data,A,[0,50],{{'LHS','RTO','RHS','LTO'}},'RCS07',[],'summed',0);

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
    channel_f_matrix.Left = cell(1,length(signal_analysis_data.Left.Chan_Names));
    channel_df_between_matrix.Left = nan(1,length(signal_analysis_data.Right.Chan_Names));
    channel_df_within_matrix.Left = nan(1,length(signal_analysis_data.Right.Chan_Names));
%     channel_mult_compare_matrix.Left = cell(1,size(event_pairs,1));
    channel_mult_compare_matrix.Left = cell(1,length(signal_analysis_data.Left.Chan_Names));
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
            [p,anova_table,stats] = anova1(X,groups,'off');
            [c,~,~,gnames] = multcompare(stats,'Display','off');
            channel_anova_matrix.Left{i}(freq_bin_inds.Left(n,1),freq_bin_inds.Left(n,2)) = p;
            channel_f_matrix.Left{i}(freq_bin_inds.Left(n,1),freq_bin_inds.Left(n,2)) = anova_table{2,5};
            if isnan(channel_df_between_matrix.Left(i))
                channel_df_between_matrix.Left(i) = anova_table{2,3};
            end
            if isnan(channel_df_within_matrix.Left(i))
                channel_df_within_matrix.Left(i) = anova_table{3,3};
            end
            for r = 1:size(c,1)
                mat_ind = find(and(strcmp(gnames(c(r,1)),event_pairs(:,1)),strcmp(gnames(c(r,2)),event_pairs(:,2))));
                if isempty(channel_mult_compare_matrix.Left{i})
                    channel_mult_compare_matrix.Left{i} = cellmat(size(event_pairs,1),1,sum(diff(freq_bin_inds.Left,1,2)==0),sum(diff(freq_bin_inds.Left,1,2)==0),nan);
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
    channel_f_matrix.Right = cell(1,length(signal_analysis_data.Right.Chan_Names));
    channel_df_between_matrix.Right = nan(1,length(signal_analysis_data.Right.Chan_Names));
    channel_df_within_matrix.Right = nan(1,length(signal_analysis_data.Right.Chan_Names));
%     channel_mult_compare_matrix.Right = cell(1,size(event_pairs,1));
    channel_mult_compare_matrix.Right = cell(1,length(signal_analysis_data.Right.Chan_Names));
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
            [p,anova_table,stats] = anova1(X,groups,'off');
            [c,~,~,gnames] = multcompare(stats,'Display','off');
            channel_anova_matrix.Right{i}(freq_bin_inds.Right(n,1),freq_bin_inds.Right(n,2)) = p;
            channel_f_matrix.Right{i}(freq_bin_inds.Right(n,1),freq_bin_inds.Right(n,2)) = anova_table{2,5};
            if isnan(channel_df_between_matrix.Right(i))
                channel_df_between_matrix.Right(i) = anova_table{2,3};
            end
            if isnan(channel_df_within_matrix.Right(i))
                channel_df_within_matrix.Right(i) = anova_table{3,3};
            end
            for r = 1:size(c,1)
                mat_ind = find(and(strcmp(gnames(c(r,1)),event_pairs(:,1)),strcmp(gnames(c(r,2)),event_pairs(:,2))));
                if isempty(channel_mult_compare_matrix.Right{i})
                    channel_mult_compare_matrix.Right{i} = cellmat(size(event_pairs,1),1,sum(diff(freq_bin_inds.Right,1,2)==0),sum(diff(freq_bin_inds.Right,1,2)==0),nan);
                end
                channel_mult_compare_matrix.Right{i}{mat_ind}(freq_bin_inds.Right(n,1),freq_bin_inds.Right(n,2)) = c(r,6);
            end
        end
    end
end

%% Create table of statistic output
% multiple-comparision table
subject_ID = [];
side =[];
contact = [];
freq1 = [];
freq2 = [];
event1 = []; %repelem({event_pairs{:,1}},n_vals)';
event2 = []; %repelem({event_pairs{:,2}},n_vals)';
pVals = [];
if isfield(signal_analysis_data,'Left')
    n_vals = size(freq_bin_inds.Left,1);
    for x = 1:length(signal_analysis_data.Left.Chan_Names)
        inds = strfind(signal_analysis_data.Left.Chan_Names{x},' ');
        chan_name = {signal_analysis_data.Left.Chan_Names{x}(1:inds(1)-1)};
        for y = 1:size(event_pairs,1)
            temp = nan(n_vals,1);
            for z = 1:n_vals
                temp(z) = channel_mult_compare_matrix.Left{x}{y}(freq_bin_inds.Left(z,1),freq_bin_inds.Left(z,2));
            end
            contact = [contact;repmat(chan_name,n_vals,1)];
            side = [side;repmat('L',n_vals,1)];
            event1 = [event1;repelem({event_pairs{y,1}},n_vals)'];
            event2 = [event2;repelem({event_pairs{y,2}},n_vals)'];
            freq1 = [freq1;signal_analysis_data.Left.Freq_Values{x}(freq_bin_inds.Left(:,1))];
            freq2 = [freq2;signal_analysis_data.Left.Freq_Values{x}(freq_bin_inds.Left(:,2))];
            pVals = [pVals;temp];
        end
    end
    subject_ID = [subject_ID;repmat({subjectID},n_vals*size(event_pairs,1)*length(signal_analysis_data.Left.Chan_Names),1)];
end

if isfield(signal_analysis_data,'Right')
    n_vals = size(freq_bin_inds.Right,1);
    for x = 1:length(signal_analysis_data.Right.Chan_Names)
        inds = strfind(signal_analysis_data.Right.Chan_Names{x},' ');
        chan_name = {signal_analysis_data.Right.Chan_Names{x}(1:inds(1)-1)};
        for y = 1:size(event_pairs,1)
            temp = nan(n_vals,1);
            for z = 1:n_vals
                temp(z) = channel_mult_compare_matrix.Right{x}{y}(freq_bin_inds.Right(z,1),freq_bin_inds.Right(z,2));
            end
            contact = [contact;repmat(chan_name,n_vals,1)];
            side = [side;repmat('R',n_vals,1)];
            event1 = [event1;repelem({event_pairs{y,1}},n_vals)'];
            event2 = [event2;repelem({event_pairs{y,2}},n_vals)'];
            freq1 = [freq1;signal_analysis_data.Right.Freq_Values{x}(freq_bin_inds.Right(:,1))];
            freq2 = [freq2;signal_analysis_data.Right.Freq_Values{x}(freq_bin_inds.Right(:,2))];
            pVals = [pVals;temp];
        end
    end
    subject_ID = [subject_ID;repmat({subjectID},n_vals*size(event_pairs,1)*length(signal_analysis_data.Right.Chan_Names),1)];
end
multiple_comp_table = table(subject_ID,side,contact,event1,event2,freq1,freq2,pVals,'VariableNames',{'SubjectID','Side','Contact','GaitEvent1','GaitEvent2','Freq1','Freq2','pVal'});

% ANOVA table
subject_ID = [];
side =[];
contact = [];
freq1 = [];
freq2 = [];
df_between = [];
df_within = [];
fStat = [];
pVals = [];

if isfield(signal_analysis_data,'Left')
    n_vals = size(freq_bin_inds.Left,1);
    for x = 1:length(signal_analysis_data.Left.Chan_Names)
        inds = strfind(signal_analysis_data.Left.Chan_Names{x},' ');
        chan_name = {signal_analysis_data.Left.Chan_Names{x}(1:inds(1)-1)};
        temp = nan(n_vals,1);
        temp2 = nan(n_vals,1);
        for z = 1:n_vals
            temp(z) = channel_anova_matrix.Left{x}(freq_bin_inds.Left(z,1),freq_bin_inds.Left(z,2));
            temp2(z) = channel_f_matrix.Left{x}(freq_bin_inds.Left(z,1),freq_bin_inds.Left(z,2));
        end
        contact = [contact;repmat(chan_name,n_vals,1)];
        side = [side;repmat('L',n_vals,1)];
        freq1 = [freq1;signal_analysis_data.Left.Freq_Values{x}(freq_bin_inds.Left(:,1))];
        freq2 = [freq2;signal_analysis_data.Left.Freq_Values{x}(freq_bin_inds.Left(:,2))];
        pVals = [pVals;temp];
        fStat = [fStat;temp2];
        df_between = [df_between;repelem(channel_df_between_matrix(x),n_vals)'];
        df_within = [df_within;repelem(channel_df_within_matrix(x),n_vals)'];
    end
    subject_ID = [subject_ID;repmat({subjectID},n_vals*length(signal_analysis_data.Left.Chan_Names),1)];
end

if isfield(signal_analysis_data,'Right')
    n_vals = size(freq_bin_inds.Right,1);
    for x = 1:length(signal_analysis_data.Right.Chan_Names)
        inds = strfind(signal_analysis_data.Right.Chan_Names{x},' ');
        chan_name = {signal_analysis_data.Right.Chan_Names{x}(1:inds(1)-1)};
        for y = 1:size(event_pairs,1)
            temp = nan(n_vals,1);
            temp2 = nan(n_vals,1);
            for z = 1:n_vals
                temp(z) = channel_anova_matrix.Right{x}(freq_bin_inds.Right(z,1),freq_bin_inds.Right(z,2));
                temp2(z) = channel_f_matrix.Right{x}(freq_bin_inds.Right(z,1),freq_bin_inds.Right(z,2));
            end
            contact = [contact;repmat(chan_name,n_vals,1)];
            side = [side;repmat('R',n_vals,1)];
            freq1 = [freq1;signal_analysis_data.Right.Freq_Values{x}(freq_bin_inds.Right(:,1))];
            freq2 = [freq2;signal_analysis_data.Right.Freq_Values{x}(freq_bin_inds.Right(:,2))];
            pVals = [pVals;temp];
            fStat = [fStat;temp2];
            df_between = [df_between;repelem(channel_df_between_matrix(x),n_vals)'];
            df_within = [df_within;repelem(channel_df_within_matrix(x),n_vals)'];
        end
    end
    subject_ID = [subject_ID;repmat({subjectID},n_vals*length(signal_analysis_data.Right.Chan_Names),1)];
end
anova_table = table(subject_ID,side,contact,freq1,freq2,df_between,df_within,fStat,pVals,'VariableNames',{'SubjectID','Side','Contact','Freq1','Freq2','DF_Between','DF_Within','fStat','pVal'});

%% Plot to visualize p-vals
fig_vec = [];
% ANOVA Left
if strcmp(save_type,'all') || strcmp(save_type,'anova')
    if isfield(channel_anova_matrix,'Left')
        for i = 1:length(channel_anova_matrix.Left)
            fig_vec(end+1) = figure;
            ax_hand = pcolor(1:size(channel_anova_matrix.Left{i},1),1:size(channel_anova_matrix.Left{i},2),channel_anova_matrix.Left{i});
            colormap(flipud(jet));
            caxis([0,1]);
            set(gca,'YDir',"reverse");
            title({subjectID;'Left';signal_analysis_data.Left.Chan_Names{i};'ANOVA'});
%             title({subjectID;['Left Brain: Contact ',signal_analysis_data.Left.Chan_Names{i}(1:4),' ANOVA']}); % For paper figure
            xlabel('End bin');
            ylabel('Start bin');
            c_hand = colorbar;
            c_hand.Title.String = 'p-val';
        end
    end
    
    % ANOVA Right
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
end

% Mult-compare Left
if strcmp(save_type,'all') || strcmp(save_type,'mult_compare')
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
    
    % Mult-Compare Right
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
end

% Summed Mult-Compare Left (excluding gait events close in time)
if strcmp(save_type,'all') || strcmp(save_type,'summed') || strcmp(save_type,'summed_exclude')
    if isfield(channel_mult_compare_matrix,'Left')
        exclude_inds = findExclusionInds(event_pairs);
        for i = 1:length(signal_analysis_data.Left.Chan_Names)
            fig_vec(end+1) = figure;
            ind_array = 1:length(channel_mult_compare_matrix.Left{i});
            ind_array(exclude_inds) = [];
            extracted_matrix = cell2mat(arrayfun(@(x)permute(x{:},[3 1 2]),channel_mult_compare_matrix.Left{i}(ind_array),'UniformOutput',false));
            summed_matrix = squeeze(sum(extracted_matrix,1));
            pcolor(1:size(summed_matrix,1),1:size(summed_matrix,2),summed_matrix);
            colormap(flipud(jet));
            set(gca,'YDir',"reverse");
            xlabel('End bin');
            ylabel('Start bin');
            caxis([0.02,1].*length(ind_array));
            set(gca,'ColorScale','log');
            c_hand = colorbar;
            c_hand.Title.String = '\Sigmap-val';
%             c_hand.Ticks = logspace(log10(0.05*length(ind_array)),log10(1*length(ind_array)),10);
            c_hand.Ticks = [0.1,0.2,0.4,0.8,1.5,2.5,4];
            title({subjectID;'Left';signal_analysis_data.Left.Chan_Names{i};'Summed Multiple Comparison Excluding Close Events'});
        end
    end
    
    % Summed Mult-Compare Right (excluding gait events close in time)
    if isfield(channel_mult_compare_matrix,'Right')
        exclude_inds = findExclusionInds(event_pairs);
        for i = 1:length(signal_analysis_data.Right.Chan_Names)
            fig_vec(end+1) = figure;
            ind_array = 1:length(channel_mult_compare_matrix.Left{i});
            ind_array(exclude_inds) = [];
            extracted_matrix = cell2mat(arrayfun(@(x)permute(x{:},[3 1 2]),channel_mult_compare_matrix.Right{i}(ind_array),'UniformOutput',false));
            summed_matrix = squeeze(sum(extracted_matrix,1));
            pcolor(1:size(summed_matrix,1),1:size(summed_matrix,2),summed_matrix);
            colormap(flipud(jet));
            set(gca,'YDir',"reverse");
            xlabel('End bin');
            ylabel('Start bin');
            caxis([0.02,1].*length(ind_array));
            set(gca,'ColorScale','log');
            c_hand = colorbar;
            c_hand.Title.String = '\Sigmap-val';
%             c_hand.Ticks = logspace(log10(0.05*length(ind_array)),log10(1*length(ind_array)),10);
            c_hand.Ticks = [0.1,0.2,0.4,0.8,1.5,2.5,4];
            title({subjectID;'Right';signal_analysis_data.Right.Chan_Names{i};'Summed Multiple Comparison Excluding Close Events'});
        end
    end
end

% Summed Mult-Compare Left (all events)
if strcmp(save_type,'all') || strcmp(save_type,'summed') || strcmp(save_type,'summed_all')
    if isfield(channel_mult_compare_matrix,'Left')
        for i = 1:length(signal_analysis_data.Left.Chan_Names)
            fig_vec(end+1) = figure;
            ind_array = 1:length(channel_mult_compare_matrix.Left{i});
            extracted_matrix = cell2mat(arrayfun(@(x)permute(x{:},[3 1 2]),channel_mult_compare_matrix.Left{i}(ind_array),'UniformOutput',false));
            summed_matrix = squeeze(sum(extracted_matrix,1));
            pcolor(1:size(summed_matrix,1),1:size(summed_matrix,2),summed_matrix);
            colormap(flipud(jet));
            set(gca,'YDir',"reverse");
            xlabel('End bin');
            ylabel('Start bin');
            caxis([0.022,1].*length(ind_array));
            set(gca,'ColorScale','log');
            c_hand = colorbar;
            c_hand.Title.String = '\Sigmap-val';
%             c_hand.Ticks = logspace(log10(0.05*length(ind_array)),log10(1*length(ind_array)),10);
            c_hand.Ticks = [0.2,0.3,0.5,0.8,1.5,2.5,4,6];
            title({subjectID;'Left';signal_analysis_data.Left.Chan_Names{i};'Summed Multiple Comparison All Events'});
        end
    end
    
    % Summed Mult-Compare Right (all events)
    if isfield(channel_mult_compare_matrix,'Right')
        for i = 1:length(signal_analysis_data.Right.Chan_Names)
            fig_vec(end+1) = figure;
            ind_array = 1:length(channel_mult_compare_matrix.Left{i});
            extracted_matrix = cell2mat(arrayfun(@(x)permute(x{:},[3 1 2]),channel_mult_compare_matrix.Right{i}(ind_array),'UniformOutput',false));
            summed_matrix = squeeze(sum(extracted_matrix,1));
            pcolor(1:size(summed_matrix,1),1:size(summed_matrix,2),summed_matrix);
            colormap(flipud(jet));
            set(gca,'YDir',"reverse");
            xlabel('End bin');
            ylabel('Start bin');
            caxis([0.022,1].*length(ind_array));
            set(gca,'ColorScale','log');
            c_hand = colorbar;
            c_hand.Title.String = '\Sigmap-val';
%             c_hand.Ticks = logspace(log10(0.05*length(ind_array)),log10(1*length(ind_array)),10);
            c_hand.Ticks = [0.2,0.3,0.5,0.8,1.5,2.5,4,6];
            title({subjectID;'Right';signal_analysis_data.Right.Chan_Names{i};'Summed Multiple Comparison All Events'});
        end
    end
end

%% Save stat table
if save_flag && strcmp(save_type,'stat_table')
    if isempty(save_dir)
        save_dir = uigetdir();
    end
    
    % check if saving folders exist
    if ~isfolder(fullfile(save_dir,'GaitBiomarkerSearch'))
        mkdir(fullfile(save_dir,'GaitBiomarkerSearch'));
    end
    
    if ~isfolder(fullfile(save_dir,'GaitBiomarkerSearch',aligned_data.stim_condition))
        mkdir(fullfile(save_dir,'GaitBiomarkerSearch',aligned_data.stim_condition))
    end
    
    if strcmp(analysis_type,'FT')
        if ~isfolder(fullfile(save_dir,'GaitBiomarkerSearch',aligned_data.stim_condition,'FT'))
            mkdir(fullfile(save_dir,'GaitBiomarkerSearch',aligned_data.stim_condition,'FT'));
            writetable(multiple_comp_table,[fullfile(save_dir,'GaitBiomarkerSearch',aligned_data.stim_condition,'FT'),'/',subjectID,'_multi-comp_stat_table.csv']);
            writetable(anova_table,[fullfile(save_dir,'GaitBiomarkerSearch',aligned_data.stim_condition,'FT'),'/',subjectID,'_anova_stat_table.csv']);
        else
            writetable(multiple_comp_table,[fullfile(save_dir,'GaitBiomarkerSearch',aligned_data.stim_condition,'FT'),'/',subjectID,'_multi-comp_stat_table.csv']);
            writetable(anova_table,[fullfile(save_dir,'GaitBiomarkerSearch',aligned_data.stim_condition,'FT'),'/',subjectID,'_anova_stat_table.csv']);
        end
    elseif strcmp(analysis_type,'CWT')
        if ~isfolder(fullfile(save_dir,'GaitBiomarkerSearch',aligned_data.stim_condition,'CWT'))
            mkdir(fullfile(save_dir,'GaitBiomarkerSearch',aligned_data.stim_condition,'CWT'));
            writetable(multiple_comp_table,[fullfile(save_dir,'GaitBiomarkerSearch',aligned_data.stim_condition,'CWT'),'/',subjectID,'_multi-comp_stat_table.csv']);
            writetable(anova_table,[fullfile(save_dir,'GaitBiomarkerSearch',aligned_data.stim_condition,'CWT'),'/',subjectID,'_anova_stat_table.csv']);
        else
            writetable(multiple_comp_table,[fullfile(save_dir,'GaitBiomarkerSearch',aligned_data.stim_condition,'CWT'),'/',subjectID,'_multi-comp_stat_table.csv']);
            writetable(anova_table,[fullfile(save_dir,'GaitBiomarkerSearch',aligned_data.stim_condition,'CWT'),'/',subjectID,'_anova_stat_table.csv']);
        end
    end
end


%% Save plots
if save_flag && ~strcmp(save_type,'stat_table')
    if isempty(save_dir)
        save_dir = uigetdir();
    end
    figure_format(17.4,15.2,10,'centimeters');
%     figure_format(12.9,11.2,10,'centimeters'); % For paper figure
    
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
if freq_vals(2) < freq_vals(1)
    freq_vals = flipud(freq_vals);
end

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

function [exclude_inds] = findExclusionInds(event_pairs)
exclude_inds(1) = find((strcmp(event_pairs(:,1),'LHS') | strcmp(event_pairs(:,1),'RTO')) & (strcmp(event_pairs(:,2),'LHS') | strcmp(event_pairs(:,2),'RTO')));
exclude_inds(2) = find((strcmp(event_pairs(:,1),'RHS') | strcmp(event_pairs(:,1),'LTO')) & (strcmp(event_pairs(:,2),'RHS') | strcmp(event_pairs(:,2),'LTO')));
end