function calcGrandAvgGaitCycle(fileList,varargin)

%% Option variables
for i = 1:2:nargin-2
    switch varargin{i}
        case 'subjectID'
            subjectID = varargin{i+1};
        case 'nPercentBins'
            nPercentBins = varargin{i+1};
        case 'gcStartEvent'
            gcStartEvent = varargin{i+1};
        case 'gaitEventRange'
            gaitEventRange = varargin{i+1};
    end
end

% Set default options if not passed in by user
if ~exist('subjectID','var') || isempty(subjectID)
    subjectID = {};             % TODO, SET UP DEFAULT NAMES (I.E. SUBJECT 1, 2, ...)
end

if ~exist('nPercentBins','var') || isempty(nPercentBins)
    nPercentBins = 100;
end

if ~exist('gcStartEvent','var') || isempty(gcStartEvent)
    gcStartEvent = 'LHS';
end

if ~exist('gaitEventRange','var') || isempty(gaitEventRange)
    gaitEventRange = [];                                        % TODO, SET UP DEFAULT MATRIX
end

%% Set up variables for each recording key

% To hold gait cycle data
gc.Left.key0 = cell(1,length(fileList));  % typically ventral DBS contact
gc.Left.key1 = cell(1,length(fileList));  % typically dorsal DBS contact
gc.Left.key2 = cell(1,length(fileList));  % typically ventral cortical contact
gc.Left.key3 = cell(1,length(fileList));  % typically dorsal cortical contact

gc.Right.key0 = cell(1,length(fileList));  % typically ventral DBS contact
gc.Right.key1 = cell(1,length(fileList));  % typically dorsal DBS contact
gc.Right.key2 = cell(1,length(fileList));  % typically ventral cortical contact
gc.Right.key3 = cell(1,length(fileList));  % typically dorsal cortical contact

% To hold normalization data
normalization.Left.key0 = cell(1,length(fileList));
normalization.Left.key0 = cell(1,length(fileList));
normalization.Left.key0 = cell(1,length(fileList));
normalization.Left.key0 = cell(1,length(fileList));

normalization.Right.key0 = cell(1,length(fileList));
normalization.Right.key1 = cell(1,length(fileList));
normalization.Right.key2 = cell(1,length(fileList));
normalization.Right.key3 = cell(1,length(fileList));

% To hold output frequency of wavelet ouptut (fc stands for frequency
% cell)
fc.Left.key0 = cell(1,length(fileList));
fc.Left.key1 = cell(1,length(fileList));
fc.Left.key2 = cell(1,length(fileList));
fc.Left.key3 = cell(1,length(fileList));

fc.Right.key0 = cell(1,length(fileList));
fc.Right.key1 = cell(1,length(fileList));
fc.Right.key2 = cell(1,length(fileList));
fc.Right.key3 = cell(1,length(fileList));

% Used to index correctly
left_count = 1;
right_count = 1;

% Grab gait cycles
for i = 1:length(fileList)
    load(fileList{i});
    gaitEventsSorted = sortGaitEvents(aligned_data.gait_events,gcStartEvent);
    signalAnalysisData = calcRCS_CWT(aligned_data);
    
    gaitEventRange = [];
    if contains(stn_patients{i},'RCS07')
        gaitEventRange(1) = 1;
        gaitEventRange(2) = height(gaitEventsSorted);
    elseif contains(stn_patients{i},'RCS12')
        gaitEventRange(1) = find(gaitEventsSorted.LHS > 51,1,'first');
        gaitEventRange(2) = find(gaitEventsSorted.LHS < 234,1,'last');
    elseif contains(stn_patients{i},'RCS15')
        gaitEventRange(1) = 1;
        gaitEventRange(2) = height(gaitEventsSorted);
    end
    
    if isfield(signalAnalysisData,'Left')
        walking_start_ind = find(signalAnalysisData.Left.Time{1} >= min(gaitEventsSorted{gaitEventRange(1),:})-1,1,'first');
        walking_end_ind = find(signalAnalysisData.Left.Time{1} <= max(gaitEventsSorted{gaitEventRange(2),:}),1,'last');
        gait_cycle_mat_left = zeros(length(signalAnalysisData.Left.Freq_Values{STN_contact_pairs_ind(i,1)}),nPercentBins,1);
        count = 1;
        for k = gaitEventRange(1):gaitEventRange(2)-1
            if ~isnan(gaitEventsSorted.(gcStartEvent)(k)) && ~isnan(gaitEventsSorted.(gcStartEvent)(k+1)) && (diff(gaitEventsSorted.(gcStartEvent)([k,k+1])) < 2)
                [~,start_ind] = min(abs(signalAnalysisData.Left.Time{STN_contact_pairs_ind(i,1)}-gaitEventsSorted.(gcStartEvent)(k)));
                [~,end_ind] = min(abs(signalAnalysisData.Left.Time{STN_contact_pairs_ind(i,1)}-gaitEventsSorted.(gcStartEvent)(k+1)));
                data_snip = abs(signalAnalysisData.Left.Values{STN_contact_pairs_ind(i,1)}(:,start_ind:end_ind));
                
                if sum(isinf(data_snip),'all') == 0
                    percent_inds = round(linspace(1,size(data_snip,2),nPercentBins+1));
                    for m = 1:length(percent_inds)-1
                        if m == 1
                            gait_cycle_mat_left(:,m,count) = mean(data_snip(:,percent_inds(m):percent_inds(m+1)),2);
                        else
                            gait_cycle_mat_left(:,m,count) = mean(data_snip(:,percent_inds(m)+1:percent_inds(m+1)),2);
                        end
                    end
                    count = count + 1;
                end
            end
        end
        gait_cycles_STN.Left{left_count} = gait_cycle_mat_left;
        normalization_STN.Left{left_count} = [mean(abs(signalAnalysisData.Left.Values{STN_contact_pairs_ind(i,1)}(:,walking_start_ind:walking_end_ind)),2,'omitnan'),...
            std(abs(signalAnalysisData.Left.Values{STN_contact_pairs_ind(i,1)}(:,walking_start_ind:walking_end_ind)),0,2,'omitnan'),...
            median(abs(signalAnalysisData.Left.Values{STN_contact_pairs_ind(i,1)}(:,walking_start_ind:walking_end_ind)),2,'omitnan')];
        fc.Left{left_count} = signalAnalysisData.Left.Freq_Values{STN_contact_pairs_ind(i,1)};
        left_count = left_count + 1;
    end
    
    if isfield(signalAnalysisData,'Right')
        walking_start_ind = find(signalAnalysisData.Right.Time{1} >= min(gaitEventsSorted{gaitEventRange(1),:})-1,1,'first');
        walking_end_ind = find(signalAnalysisData.Right.Time{1} <= max(gaitEventsSorted{gaitEventRange(2),:}),1,'last');
        gait_cycle_mat_right = zeros(length(signalAnalysisData.Right.Freq_Values{STN_contact_pairs_ind(i,2)}),nPercentBins,1);
        count = 1;
        for k = gaitEventRange(1):gaitEventRange(2)-1
            if ~isnan(gaitEventsSorted.(gcStartEvent)(k)) && ~isnan(gaitEventsSorted.(gcStartEvent)(k+1)) && (diff(gaitEventsSorted.(gcStartEvent)([k,k+1])) < 2)
                [~,start_ind] = min(abs(signalAnalysisData.Right.Time{STN_contact_pairs_ind(i,2)}-gaitEventsSorted.(gcStartEvent)(k)));
                [~,end_ind] = min(abs(signalAnalysisData.Right.Time{STN_contact_pairs_ind(i,2)}-gaitEventsSorted.(gcStartEvent)(k+1)));
                data_snip = abs(signalAnalysisData.Right.Values{STN_contact_pairs_ind(i,2)}(:,start_ind:end_ind));
                
                if sum(isinf(data_snip),'all') == 0
                    percent_inds = round(linspace(1,size(data_snip,2),nPercentBins+1));
                    for m = 1:length(percent_inds)-1
                        if m == 1
                            gait_cycle_mat_right(:,m,count) = mean(data_snip(:,percent_inds(m):percent_inds(m+1)),2);
                        else
                            gait_cycle_mat_right(:,m,count) = mean(data_snip(:,percent_inds(m)+1:percent_inds(m+1)),2);
                        end
                    end
                    count = count + 1;
                end
            end
        end
        gait_cycles_STN.Right{right_count} = gait_cycle_mat_right;
        normalization_STN.Right{right_count} = [mean(abs(signalAnalysisData.Right.Values{STN_contact_pairs_ind(i,2)}(:,walking_start_ind:walking_end_ind)),2,'omitnan'),...
            std(abs(signalAnalysisData.Right.Values{STN_contact_pairs_ind(i,2)}(:,walking_start_ind:walking_end_ind)),0,2,'omitnan'),...
            median(abs(signalAnalysisData.Right.Values{STN_contact_pairs_ind(i,2)}(:,walking_start_ind:walking_end_ind)),2,'omitnan')];
        fc.Right{right_count} = signalAnalysisData.Right.Freq_Values{STN_contact_pairs_ind(i,2)};
        right_count = right_count + 1;
    end
end

% Normalize and average all together
normalized_data_STN.Left = [];
normalized_data_STN.Right = [];

for i = 1:length(gait_cycles_STN.Left)
    start_ind = find(fc.Left{i}<=100,1,'first');
    end_ind = find(fc.Left{i}>=0.1,1,'last');
    normalization_mat = nan(size(gait_cycles_STN.Left{i}));
    for j = 1:size(gait_cycles_STN.Left{i},3)
%         normalization_mat(:,:,j) = (gait_cycles_S1.Left{i}(:,:,j)-normalization_S1.Left{i}(:,1))./normalization_S1.Left{i}(:,1);    % deltaF/F (percent change)
%         normalization_mat(:,:,j) = (gait_cycles_S1.Left{i}(:,:,j)-normalization_S1.Left{i}(:,1))./normalization_S1.Left{i}(:,2);    % Z-score
        normalization_mat(:,:,j) = normalizeData(gait_cycles_STN.Left{i}(:,:,j),'zscore',normalization_STN.Left{i});
    end
    normalized_data_STN.Left = cat(3,normalized_data_STN.Left,normalization_mat(start_ind:end_ind,:,:));
end

for i = 1:length(gait_cycles_STN.Right)
    start_ind = find(fc.Right{i}<=100,1,'first');
    end_ind = find(fc.Right{i}>=0.1,1,'last');
    normalization_mat = nan(size(gait_cycles_STN.Right{i}));
    for j = 1:size(gait_cycles_STN.Right{i},3)
%         normalization_mat(:,:,j) = (gait_cycles_S1.Right{i}(:,:,j)-normalization_S1.Right{i}(:,1))./normalization_S1.Right{i}(:,1);    % deltaF/F (percent change)
%         normalization_mat(:,:,j) = (gait_cycles_S1.Right{i}(:,:,j)-normalization_S1.Right{i}(:,1))./normalization_S1.Right{i}(:,2);    % Z-score
        normalization_mat(:,:,j) = normalizeData(gait_cycles_STN.Right{i}(:,:,j),'zscore',normalization_STN.Right{i});
    end
    normalized_data_STN.Right = cat(3,normalized_data_STN.Right,normalization_mat(start_ind:end_ind,:,:));
end

grand_average_STN.Left = mean(normalized_data_STN.Left,3);
grand_average_STN.Right = mean(normalized_data_STN.Right,3);
start_ind = find(fc.Left{1}<=100,1,'first');
end_ind = find(fc.Left{1}>=0.1,1,'last');
freq_vec = fc.Left{1}(start_ind:end_ind);

figure(1); clf;
if(~isempty(grand_average_STN.Left))
    a(1) = subplot(4,2,1);
    ax(1) = pcolor(1:100,log2(freq_vec),grand_average_STN.Left);
    ticks = logspace(log10(2.5),log10(50),5);
    ax(1).Parent.YTick = log2(ticks);
    ax(1).Parent.YTickLabel = sprintf('%1.2f\n',round(ticks,2));
    ylim([log2(2.5),log2(50)]);
%     xticks([1,10,20,30,40,50,60,70,80,90,100]);
%     xticklabels({'0','10','20','30','40','50','60','70','80','90','100'});
    shading interp;
    caxis(a(1),[-0.4,0.4]);
    ylabel({'Left';'Frequency (Hz)'});
    title('Ventral STN');
    cb_obj_vSTN_gc = colorbar;
    cb_obj_vSTN_gc.Label.String = 'Normalized Power';
    cb_obj_vSTN_gc.Ticks = [-0.40:0.10:0.40];
%     cb_obj_gc.Label.String = 'Z-Score';
end
% 
if(~isempty(grand_average_STN.Right))
    a(2) = subplot(4,2,3);
    ax(2) = pcolor(1:100,log2(freq_vec),grand_average_STN.Right);
    ticks = logspace(log10(2.5),log10(50),5);
    ax(2).Parent.YTick = log2(ticks);
    ax(2).Parent.YTickLabel = sprintf('%1.2f\n',round(ticks,2));
    ylim([log2(2.5),log2(50)]);
    xticks([1,10,20,30,40,50,60,70,80,90,100]);
    xticklabels({'0','10','20','30','40','50','60','70','80','90','100'});
    shading interp;
    caxis(a(2),[-0.4,0.4]);
    ylabel({'Right';'Frequency (Hz)'});
    xlabel('% Gait Cycle');
end

%% dorsal STN
STN_contact_pairs_ind = [2,2;2,2;2,2];
gait_cycles_STN.Left = cell(1,length(stn_patients));
gait_cycles_STN.Right = cell(1,length(stn_patients));

normalization_STN.Left = cell(1,length(stn_patients));
normalization_STN.Right = cell(1,length(stn_patients));

fc.Left = cell(1,length(stn_patients));
fc.Right = cell(1,length(stn_patients));

left_count = 1;
right_count = 1;

% Grab gait cycles
for i = 1:length(stn_patients)
    load(stn_patients{i});
    gaitEventsSorted = sortGaitEvents(aligned_data.gait_events,gcStartEvent);
    signalAnalysisData = calcRCS_CWT(aligned_data);
    
    gaitEventRange = [];
    if contains(stn_patients{i},'RCS07')
        gaitEventRange(1) = 1;
        gaitEventRange(2) = height(gaitEventsSorted);
    elseif contains(stn_patients{i},'RCS12')
        gaitEventRange(1) = find(gaitEventsSorted.LHS > 51,1,'first');
        gaitEventRange(2) = find(gaitEventsSorted.LHS < 234,1,'last');
    elseif contains(stn_patients{i},'RCS15')
        gaitEventRange(1) = 1;
        gaitEventRange(2) = height(gaitEventsSorted);
    end
    
    if isfield(signalAnalysisData,'Left')
        walking_start_ind = find(signalAnalysisData.Left.Time{1} >= min(gaitEventsSorted{gaitEventRange(1),:})-1,1,'first');
        walking_end_ind = find(signalAnalysisData.Left.Time{1} <= max(gaitEventsSorted{gaitEventRange(2),:}),1,'last');
        gait_cycle_mat_left = zeros(length(signalAnalysisData.Left.Freq_Values{STN_contact_pairs_ind(i,1)}),nPercentBins,1);
        count = 1;
        for k = gaitEventRange(1):gaitEventRange(2)-1
            if ~isnan(gaitEventsSorted.(gcStartEvent)(k)) && ~isnan(gaitEventsSorted.(gcStartEvent)(k+1)) && (diff(gaitEventsSorted.(gcStartEvent)([k,k+1])) < 2)
                [~,start_ind] = min(abs(signalAnalysisData.Left.Time{STN_contact_pairs_ind(i,1)}-gaitEventsSorted.(gcStartEvent)(k)));
                [~,end_ind] = min(abs(signalAnalysisData.Left.Time{STN_contact_pairs_ind(i,1)}-gaitEventsSorted.(gcStartEvent)(k+1)));
                data_snip = abs(signalAnalysisData.Left.Values{STN_contact_pairs_ind(i,1)}(:,start_ind:end_ind));
                
                if sum(isinf(data_snip),'all') == 0
                    percent_inds = round(linspace(1,size(data_snip,2),nPercentBins+1));
                    for m = 1:length(percent_inds)-1
                        if m == 1
                            gait_cycle_mat_left(:,m,count) = mean(data_snip(:,percent_inds(m):percent_inds(m+1)),2);
                        else
                            gait_cycle_mat_left(:,m,count) = mean(data_snip(:,percent_inds(m)+1:percent_inds(m+1)),2);
                        end
                    end
                    count = count + 1;
                end
            end
        end
        gait_cycles_STN.Left{left_count} = gait_cycle_mat_left;
        normalization_STN.Left{left_count} = [mean(abs(signalAnalysisData.Left.Values{STN_contact_pairs_ind(i,1)}(:,walking_start_ind:walking_end_ind)),2,'omitnan'),...
            std(abs(signalAnalysisData.Left.Values{STN_contact_pairs_ind(i,1)}(:,walking_start_ind:walking_end_ind)),0,2,'omitnan'),...
            median(abs(signalAnalysisData.Left.Values{STN_contact_pairs_ind(i,1)}(:,walking_start_ind:walking_end_ind)),2,'omitnan')];
        fc.Left{left_count} = signalAnalysisData.Left.Freq_Values{STN_contact_pairs_ind(i,1)};
        left_count = left_count + 1;
    end
    
    if isfield(signalAnalysisData,'Right')
        walking_start_ind = find(signalAnalysisData.Right.Time{1} >= min(gaitEventsSorted{gaitEventRange(1),:})-1,1,'first');
        walking_end_ind = find(signalAnalysisData.Right.Time{1} <= max(gaitEventsSorted{gaitEventRange(2),:}),1,'last');
        gait_cycle_mat_right = zeros(length(signalAnalysisData.Right.Freq_Values{STN_contact_pairs_ind(i,2)}),nPercentBins,1);
        count = 1;
        for k = gaitEventRange(1):gaitEventRange(2)-1
            if ~isnan(gaitEventsSorted.(gcStartEvent)(k)) && ~isnan(gaitEventsSorted.(gcStartEvent)(k+1)) && (diff(gaitEventsSorted.(gcStartEvent)([k,k+1])) < 2)
                [~,start_ind] = min(abs(signalAnalysisData.Right.Time{STN_contact_pairs_ind(i,2)}-gaitEventsSorted.(gcStartEvent)(k)));
                [~,end_ind] = min(abs(signalAnalysisData.Right.Time{STN_contact_pairs_ind(i,2)}-gaitEventsSorted.(gcStartEvent)(k+1)));
                data_snip = abs(signalAnalysisData.Right.Values{STN_contact_pairs_ind(i,2)}(:,start_ind:end_ind));
                
                if sum(isinf(data_snip),'all') == 0
                    percent_inds = round(linspace(1,size(data_snip,2),nPercentBins+1));
                    for m = 1:length(percent_inds)-1
                        if m == 1
                            gait_cycle_mat_right(:,m,count) = mean(data_snip(:,percent_inds(m):percent_inds(m+1)),2);
                        else
                            gait_cycle_mat_right(:,m,count) = mean(data_snip(:,percent_inds(m)+1:percent_inds(m+1)),2);
                        end
                    end
                    count = count + 1;
                end
            end
        end
        gait_cycles_STN.Right{right_count} = gait_cycle_mat_right;
        normalization_STN.Right{right_count} = [mean(abs(signalAnalysisData.Right.Values{STN_contact_pairs_ind(i,2)}(:,walking_start_ind:walking_end_ind)),2,'omitnan'),...
            std(abs(signalAnalysisData.Right.Values{STN_contact_pairs_ind(i,2)}(:,walking_start_ind:walking_end_ind)),0,2,'omitnan'),...
            median(abs(signalAnalysisData.Right.Values{STN_contact_pairs_ind(i,2)}(:,walking_start_ind:walking_end_ind)),2,'omitnan')];
        fc.Right{right_count} = signalAnalysisData.Right.Freq_Values{STN_contact_pairs_ind(i,2)};
        right_count = right_count + 1;
    end
end

% Normalize and average all together
normalized_data_STN.Left = [];
normalized_data_STN.Right = [];

for i = 1:length(gait_cycles_STN.Left)
    start_ind = find(fc.Left{i}<=100,1,'first');
    end_ind = find(fc.Left{i}>=0.1,1,'last');
    normalization_mat = nan(size(gait_cycles_STN.Left{i}));
    for j = 1:size(gait_cycles_STN.Left{i},3)
%         normalization_mat(:,:,j) = (gait_cycles_S1.Left{i}(:,:,j)-normalization_S1.Left{i}(:,1))./normalization_S1.Left{i}(:,1);    % deltaF/F (percent change)
%         normalization_mat(:,:,j) = (gait_cycles_S1.Left{i}(:,:,j)-normalization_S1.Left{i}(:,1))./normalization_S1.Left{i}(:,2);    % Z-score
        normalization_mat(:,:,j) = normalizeData(gait_cycles_STN.Left{i}(:,:,j),'zscore',normalization_STN.Left{i});
    end
    normalized_data_STN.Left = cat(3,normalized_data_STN.Left,normalization_mat(start_ind:end_ind,:,:));
end

for i = 1:length(gait_cycles_STN.Right)
    start_ind = find(fc.Right{i}<=100,1,'first');
    end_ind = find(fc.Right{i}>=0.1,1,'last');
    normalization_mat = nan(size(gait_cycles_STN.Right{i}));
    for j = 1:size(gait_cycles_STN.Right{i},3)
%         normalization_mat(:,:,j) = (gait_cycles_S1.Right{i}(:,:,j)-normalization_S1.Right{i}(:,1))./normalization_S1.Right{i}(:,1);    % deltaF/F (percent change)
%         normalization_mat(:,:,j) = (gait_cycles_S1.Right{i}(:,:,j)-normalization_S1.Right{i}(:,1))./normalization_S1.Right{i}(:,2);    % Z-score
        normalization_mat(:,:,j) = normalizeData(gait_cycles_STN.Right{i}(:,:,j),'zscore',normalization_STN.Right{i});
    end
    normalized_data_STN.Right = cat(3,normalized_data_STN.Right,normalization_mat(start_ind:end_ind,:,:));
end

grand_average_STN.Left = mean(normalized_data_STN.Left,3);
grand_average_STN.Right = mean(normalized_data_STN.Right,3);
start_ind = find(fc.Left{1}<=100,1,'first');
end_ind = find(fc.Left{1}>=0.1,1,'last');
freq_vec = fc.Left{1}(start_ind:end_ind);

if(~isempty(grand_average_STN.Left))
    a(3) = subplot(4,2,2);
    ax(3) = pcolor(1:100,log2(freq_vec),grand_average_STN.Left);
    ticks = logspace(log10(2.5),log10(50),5);
    ax(3).Parent.YTick = log2(ticks);
    ax(3).Parent.YTickLabel = sprintf('%1.2f\n',round(ticks,2));
    ylim([log2(2.5),log2(50)]);
%     xticks([1,10,20,30,40,50,60,70,80,90,100]);
%     xticklabels({'0','10','20','30','40','50','60','70','80','90','100'});
    shading interp;
    caxis(a(3),[-0.2,0.2]);
    ylabel({'Left';'Frequency (Hz)'});
    title('Dorsal STN');
    cb_obj_dSTN_gc = colorbar;
    cb_obj_dSTN_gc.Label.String = 'Normalized Power';
    cb_obj_dSTN_gc.Ticks = [-0.25:0.05:0.25];
%     cb_obj_gc.Label.String = 'Z-Score';
end
% 
if(~isempty(grand_average_STN.Right))
    a(4) = subplot(4,2,4);
    ax(4) = pcolor(1:100,log2(freq_vec),grand_average_STN.Right);
    ticks = logspace(log10(2.5),log10(50),5);
    ax(4).Parent.YTick = log2(ticks);
    ax(4).Parent.YTickLabel = sprintf('%1.2f\n',round(ticks,2));
    ylim([log2(2.5),log2(50)]);
    xticks([1,10,20,30,40,50,60,70,80,90,100]);
    xticklabels({'0','10','20','30','40','50','60','70','80','90','100'});
    shading interp;
    caxis(a(4),[-0.2,0.2]);
    ylabel({'Right';'Frequency (Hz)'});
    xlabel('% Gait Cycle');
end


%% S1 Grand Average Gait Cycle
gait_cycles_S1.Left = cell(1,length(stn_patients));
gait_cycles_S1.Right = cell(1,length(stn_patients));

normalization_S1.Left = cell(1,length(stn_patients));
normalization_S1.Right = cell(1,length(stn_patients));

fc.Left = cell(1,length(stn_patients));
fc.Right = cell(1,length(stn_patients));

left_count = 1;
right_count = 1;

movement_start = 10;

% Grab coherence by gait cycle
% Grab gait cycles
for i = 1:length(stn_patients)
    load(stn_patients{i});
    gaitEventsSorted = sortGaitEvents(aligned_data.gait_events,gcStartEvent);
    signalAnalysisData = calcRCS_CWT(aligned_data);
    
    gaitEventRange = [];
    if contains(stn_patients{i},'RCS07')
        gaitEventRange(1) = 1;
        gaitEventRange(2) = height(gaitEventsSorted);
    elseif contains(stn_patients{i},'RCS12')
        gaitEventRange(1) = find(gaitEventsSorted.LHS > 51,1,'first');
        gaitEventRange(2) = find(gaitEventsSorted.LHS < 234,1,'last');
    elseif contains(stn_patients{i},'RCS15')
        gaitEventRange(1) = 1;
        gaitEventRange(2) = height(gaitEventsSorted);
    end
    
    if isfield(signalAnalysisData,'Left')
        walking_start_ind = find(signalAnalysisData.Left.Time{3} >= min(gaitEventsSorted{gaitEventRange(1),:})-1,1,'first');
        walking_end_ind = find(signalAnalysisData.Left.Time{3} <= max(gaitEventsSorted{gaitEventRange(2),:}),1,'last');
        gait_cycle_mat_left = zeros(length(signalAnalysisData.Left.Freq_Values{3}),nPercentBins,1);
        count = 1;
        for k = gaitEventRange(1):gaitEventRange(2)-1
            if ~isnan(gaitEventsSorted.(gcStartEvent)(k)) && ~isnan(gaitEventsSorted.(gcStartEvent)(k+1)) && (diff(gaitEventsSorted.(gcStartEvent)([k,k+1])) < 2)
                [~,start_ind] = min(abs(signalAnalysisData.Left.Time{3}-gaitEventsSorted.(gcStartEvent)(k)));
                [~,end_ind] = min(abs(signalAnalysisData.Left.Time{3}-gaitEventsSorted.(gcStartEvent)(k+1)));
                data_snip = abs(signalAnalysisData.Left.Values{3}(:,start_ind:end_ind));
                
                if sum(isinf(data_snip),'all') == 0
                    percent_inds = round(linspace(1,size(data_snip,2),nPercentBins+1));
                    for m = 1:length(percent_inds)-1
                        if m == 1
                            gait_cycle_mat_left(:,m,count) = mean(data_snip(:,percent_inds(m):percent_inds(m+1)),2);
                        else
                            gait_cycle_mat_left(:,m,count) = mean(data_snip(:,percent_inds(m)+1:percent_inds(m+1)),2);
                        end
                    end
                    count = count + 1;
                end
            end
        end
        gait_cycles_S1.Left{left_count} = gait_cycle_mat_left;
        normalization_S1.Left{left_count} = [mean(abs(signalAnalysisData.Left.Values{3}(:,walking_start_ind:walking_end_ind)),2,'omitnan'),...
            std(abs(signalAnalysisData.Left.Values{3}(:,walking_start_ind:walking_end_ind)),0,2,'omitnan'),...
            median(abs(signalAnalysisData.Left.Values{3}(:,walking_start_ind:walking_end_ind)),2,'omitnan')];
        
        
        fc.Left{left_count} = signalAnalysisData.Left.Freq_Values{3};
        left_count = left_count + 1;
    end
    
    if isfield(signalAnalysisData,'Right')
        walking_start_ind = find(signalAnalysisData.Right.Time{3} >= min(gaitEventsSorted{gaitEventRange(1),:})-1,1,'first');
        walking_end_ind = find(signalAnalysisData.Right.Time{3} <= max(gaitEventsSorted{gaitEventRange(2),:}),1,'last');
        gait_cycle_mat_right = zeros(length(signalAnalysisData.Right.Freq_Values{3}),nPercentBins,1);
        count = 1;
        for k = gaitEventRange(1):gaitEventRange(2)-1
            if ~isnan(gaitEventsSorted.(gcStartEvent)(k)) && ~isnan(gaitEventsSorted.(gcStartEvent)(k+1)) && (diff(gaitEventsSorted.(gcStartEvent)([k,k+1])) < 2)
                [~,start_ind] = min(abs(signalAnalysisData.Right.Time{3}-gaitEventsSorted.(gcStartEvent)(k)));
                [~,end_ind] = min(abs(signalAnalysisData.Right.Time{3}-gaitEventsSorted.(gcStartEvent)(k+1)));
                data_snip = abs(signalAnalysisData.Right.Values{3}(:,start_ind:end_ind));
                
                if sum(isinf(data_snip),'all') == 0
                    percent_inds = round(linspace(1,size(data_snip,2),nPercentBins+1));
                    for m = 1:length(percent_inds)-1
                        if m == 1
                            gait_cycle_mat_right(:,m,count) = mean(data_snip(:,percent_inds(m):percent_inds(m+1)),2);
                        else
                            gait_cycle_mat_right(:,m,count) = mean(data_snip(:,percent_inds(m)+1:percent_inds(m+1)),2);
                        end
                    end
                    count = count + 1;
                end
            end
        end
        gait_cycles_S1.Right{right_count} = gait_cycle_mat_right;
        normalization_S1.Right{right_count} = [mean(abs(signalAnalysisData.Right.Values{3}(:,walking_start_ind:walking_end_ind)),2,'omitnan'),...
            std(abs(signalAnalysisData.Right.Values{3}(:,walking_start_ind:walking_end_ind)),0,2,'omitnan'),...
            median(abs(signalAnalysisData.Right.Values{3}(:,walking_start_ind:walking_end_ind)),2,'omitnan')];
        
        fc.Right{right_count} = signalAnalysisData.Right.Freq_Values{3};
        right_count = right_count + 1;
    end
end

% Normalize and average all together
normalized_data_S1.Left = [];
normalized_data_S1.Right = [];

for i = 1:length(gait_cycles_S1.Left)
    start_ind = find(fc.Left{i}<=100,1,'first');
    end_ind = find(fc.Left{i}>=0.1,1,'last');
    normalization_mat = nan(size(gait_cycles_S1.Left{i}));
    for j = 1:size(gait_cycles_S1.Left{i},3)
        normalization_mat(:,:,j) = normalizeData(gait_cycles_S1.Left{i}(:,:,j),'zscore',normalization_S1.Left{i});
    end
    normalized_data_S1.Left = cat(3,normalized_data_S1.Left,normalization_mat(start_ind:end_ind,:,:));
end

for i = 1:length(gait_cycles_S1.Right)
    start_ind = find(fc.Right{i}<=100,1,'first');
    end_ind = find(fc.Right{i}>=0.1,1,'last');
    normalization_mat = nan(size(gait_cycles_S1.Right{i}));
    for j = 1:size(gait_cycles_S1.Right{i},3)
        normalization_mat(:,:,j) = normalizeData(gait_cycles_S1.Right{i}(:,:,j),'zscore',normalization_S1.Right{i});
    end
    normalized_data_S1.Right = cat(3,normalized_data_S1.Right,normalization_mat(start_ind:end_ind,:,:));
end

grand_average_S1.Left = mean(normalized_data_S1.Left,3);
grand_average_S1.Right = mean(normalized_data_S1.Right,3);
start_ind = find(fc.Left{1}<=100,1,'first');
end_ind = find(fc.Left{1}>=0.1,1,'last');
freq_vec = fc.Left{1}(start_ind:end_ind);

if(~isempty(grand_average_S1.Left))
    a(5) = subplot(4,2,6);
    ax(5) = pcolor(1:100,log2(freq_vec),grand_average_S1.Left);
    ticks = logspace(log10(2.5),log10(50),5);
    ax(5).Parent.YTick = log2(ticks);
    ax(5).Parent.YTickLabel = sprintf('%1.2f\n',round(ticks,2));
    ylim([log2(2.5),log2(50)]);
%     xticks([1,10,20,30,40,50,60,70,80,90,100]);
%     xticklabels({'0','10','20','30','40','50','60','70','80','90','100'});
    shading interp;
    caxis(a(5),[-0.3,0.3]);
    ylabel({'Left';'Frequency (Hz)'});
    title('S1');
    cb_obj_S1_gc = colorbar;
    cb_obj_S1_gc.Label.String = 'Normalized Power';
    cb_obj_S1_gc.Ticks = [-0.3:0.1:0.3];
%     cb_obj_gc.Label.String = 'Z-Score';
end
% 
if(~isempty(grand_average_S1.Right))
    a(6) = subplot(4,2,8);
    ax(6) = pcolor(1:100,log2(freq_vec),grand_average_S1.Right);
    ticks = logspace(log10(2.5),log10(50),5);
    ax(6).Parent.YTick = log2(ticks);
    ax(6).Parent.YTickLabel = sprintf('%1.2f\n',round(ticks,2));
    ylim([log2(2.5),log2(50)]);
    xticks([1,10,20,30,40,50,60,70,80,90,100]);
    xticklabels({'0','10','20','30','40','50','60','70','80','90','100'});
    shading interp;
    caxis(a(6),[-0.3,0.3]);
    ylabel({'Right';'Frequency (Hz)'});
    xlabel('% Gait Cycle');
end

%% M1 Grand Average Gait Cycle
gait_cycles_M1.Left = cell(1,length(stn_patients));
gait_cycles_M1.Right = cell(1,length(stn_patients));

normalization_M1.Left = cell(1,length(stn_patients));
normalization_M1.Right = cell(1,length(stn_patients));

fc.Left = cell(1,length(stn_patients));
fc.Right = cell(1,length(stn_patients));

left_count = 1;
right_count = 1;

movement_start = 10;

% Grab coherence by gait cycle
% Grab gait cycles
for i = 1:length(stn_patients)
    load(stn_patients{i});
    gaitEventsSorted = sortGaitEvents(aligned_data.gait_events,gcStartEvent);
    signalAnalysisData = calcRCS_CWT(aligned_data);
    
    gaitEventRange = [];
    if contains(stn_patients{i},'RCS07')
        gaitEventRange(1) = 1;
        gaitEventRange(2) = height(gaitEventsSorted);
    elseif contains(stn_patients{i},'RCS12')
        gaitEventRange(1) = find(gaitEventsSorted.LHS > 51,1,'first');
        gaitEventRange(2) = find(gaitEventsSorted.LHS < 234,1,'last');
    elseif contains(stn_patients{i},'RCS15')
        gaitEventRange(1) = 1;
        gaitEventRange(2) = height(gaitEventsSorted);
    end
    
    if isfield(signalAnalysisData,'Left')
        walking_start_ind = find(signalAnalysisData.Left.Time{4} >= min(gaitEventsSorted{gaitEventRange(1),:})-1,1,'first');
        walking_end_ind = find(signalAnalysisData.Left.Time{4} <= max(gaitEventsSorted{gaitEventRange(2),:}),1,'last');
        gait_cycle_mat_left = zeros(length(signalAnalysisData.Left.Freq_Values{4}),nPercentBins,1);
        count = 1;
        for k = gaitEventRange(1):gaitEventRange(2)-1
            if ~isnan(gaitEventsSorted.(gcStartEvent)(k)) && ~isnan(gaitEventsSorted.(gcStartEvent)(k+1)) && (diff(gaitEventsSorted.(gcStartEvent)([k,k+1])) < 2)
                [~,start_ind] = min(abs(signalAnalysisData.Left.Time{4}-gaitEventsSorted.(gcStartEvent)(k)));
                [~,end_ind] = min(abs(signalAnalysisData.Left.Time{4}-gaitEventsSorted.(gcStartEvent)(k+1)));
                data_snip = abs(signalAnalysisData.Left.Values{4}(:,start_ind:end_ind));
                
                if sum(isinf(data_snip),'all') == 0
                    percent_inds = round(linspace(1,size(data_snip,2),nPercentBins+1));
                    for m = 1:length(percent_inds)-1
                        if m == 1
                            gait_cycle_mat_left(:,m,count) = mean(data_snip(:,percent_inds(m):percent_inds(m+1)),2);
                        else
                            gait_cycle_mat_left(:,m,count) = mean(data_snip(:,percent_inds(m)+1:percent_inds(m+1)),2);
                        end
                    end
                    count = count + 1;
                end
            end
        end
        gait_cycles_M1.Left{left_count} = gait_cycle_mat_left;
        normalization_M1.Left{left_count} = [mean(abs(signalAnalysisData.Left.Values{4}(:,walking_start_ind:walking_end_ind)),2,'omitnan'),...
            std(abs(signalAnalysisData.Left.Values{4}(:,walking_start_ind:walking_end_ind)),0,2,'omitnan'),...
            median(abs(signalAnalysisData.Left.Values{4}(:,walking_start_ind:walking_end_ind)),2,'omitnan')];
        
        
        fc.Left{left_count} = signalAnalysisData.Left.Freq_Values{4};
        left_count = left_count + 1;
    end
    
    if isfield(signalAnalysisData,'Right')
        walking_start_ind = find(signalAnalysisData.Right.Time{4} >= min(gaitEventsSorted{gaitEventRange(1),:})-1,1,'first');
        walking_end_ind = find(signalAnalysisData.Right.Time{4} <= max(gaitEventsSorted{gaitEventRange(2),:}),1,'last');
        gait_cycle_mat_right = zeros(length(signalAnalysisData.Right.Freq_Values{4}),nPercentBins,1);
        count = 1;
        for k = gaitEventRange(1):gaitEventRange(2)-1
            if ~isnan(gaitEventsSorted.(gcStartEvent)(k)) && ~isnan(gaitEventsSorted.(gcStartEvent)(k+1)) && (diff(gaitEventsSorted.(gcStartEvent)([k,k+1])) < 2)
                [~,start_ind] = min(abs(signalAnalysisData.Right.Time{4}-gaitEventsSorted.(gcStartEvent)(k)));
                [~,end_ind] = min(abs(signalAnalysisData.Right.Time{4}-gaitEventsSorted.(gcStartEvent)(k+1)));
                data_snip = abs(signalAnalysisData.Right.Values{4}(:,start_ind:end_ind));
                
                if sum(isinf(data_snip),'all') == 0
                    percent_inds = round(linspace(1,size(data_snip,2),nPercentBins+1));
                    for m = 1:length(percent_inds)-1
                        if m == 1
                            gait_cycle_mat_right(:,m,count) = mean(data_snip(:,percent_inds(m):percent_inds(m+1)),2);
                        else
                            gait_cycle_mat_right(:,m,count) = mean(data_snip(:,percent_inds(m)+1:percent_inds(m+1)),2);
                        end
                    end
                    count = count + 1;
                end
            end
        end
        gait_cycles_M1.Right{right_count} = gait_cycle_mat_right;
        normalization_M1.Right{right_count} = [mean(abs(signalAnalysisData.Right.Values{4}(:,walking_start_ind:walking_end_ind)),2,'omitnan'),...
            std(abs(signalAnalysisData.Right.Values{4}(:,walking_start_ind:walking_end_ind)),0,2,'omitnan'),...
            median(abs(signalAnalysisData.Right.Values{4}(:,walking_start_ind:walking_end_ind)),2,'omitnan')];
        
        fc.Right{right_count} = signalAnalysisData.Right.Freq_Values{4};
        right_count = right_count + 1;
    end
end

% Normalize and average all together
normalized_data_M1.Left = [];
normalized_data_M1.Right = [];

for i = 1:length(gait_cycles_M1.Left)
    start_ind = find(fc.Left{i}<=100,1,'first');
    end_ind = find(fc.Left{i}>=0.1,1,'last');
    normalization_mat = nan(size(gait_cycles_M1.Left{i}));
    for j = 1:size(gait_cycles_M1.Left{i},3)
%         normalization_mat(:,:,j) = (gait_cycles_M1.Left{i}(:,:,j)-normalization_M1.Left{i}(:,1))./normalization_M1.Left{i}(:,1);    % deltaF/F (percent change)
%         normalization_mat(:,:,j) = (gait_cycles_M1.Left{i}(:,:,j)-normalization_M1.Left{i}(:,1))./normalization_M1.Left{i}(:,2);    % Z-score
        normalization_mat(:,:,j) = normalizeData(gait_cycles_M1.Left{i}(:,:,j),'zscore',normalization_M1.Left{i});
    end
    normalized_data_M1.Left = cat(3,normalized_data_M1.Left,normalization_mat(start_ind:end_ind,:,:));
end

for i = 1:length(gait_cycles_M1.Right)
    start_ind = find(fc.Right{i}<=100,1,'first');
    end_ind = find(fc.Right{i}>=0.1,1,'last');
    normalization_mat = nan(size(gait_cycles_M1.Right{i}));
    for j = 1:size(gait_cycles_M1.Right{i},3)
%         normalization_mat(:,:,j) = (gait_cycles_M1.Right{i}(:,:,j)-normalization_M1.Right{i}(:,1))./normalization_M1.Right{i}(:,1);    % deltaF/F (percent change)
%         normalization_mat(:,:,j) = (gait_cycles_M1.Right{i}(:,:,j)-normalization_M1.Right{i}(:,1))./normalization_M1.Right{i}(:,2);    % Z-score
        normalization_mat(:,:,j) = normalizeData(gait_cycles_M1.Right{i}(:,:,j),'zscore',normalization_M1.Right{i});
    end
    normalized_data_M1.Right = cat(3,normalized_data_M1.Right,normalization_mat(start_ind:end_ind,:,:));
end

grand_average_M1.Left = mean(normalized_data_M1.Left,3);
grand_average_M1.Right = mean(normalized_data_M1.Right,3);
start_ind = find(fc.Left{1}<=100,1,'first');
end_ind = find(fc.Left{1}>=0.1,1,'last');
freq_vec = fc.Left{1}(start_ind:end_ind);

if(~isempty(grand_average_M1.Left))
    a(7) = subplot(4,2,5);
    ax(7) = pcolor(1:100,log2(freq_vec),grand_average_M1.Left);
    ticks = logspace(log10(2.5),log10(50),5);
    ax(7).Parent.YTick = log2(ticks);
    ax(7).Parent.YTickLabel = sprintf('%1.2f\n',round(ticks,2));
    ylim([log2(2.5),log2(50)]);
%     xticks([1,10,20,30,40,50,60,70,80,90,100]);
%     xticklabels({'0','10','20','30','40','50','60','70','80','90','100'});
    shading interp;
    caxis(a(7),[-0.3,0.3]);
    ylabel({'Left';'Frequency (Hz)'});
    title('M1');
    cb_obj_M1_gc = colorbar;
    cb_obj_M1_gc.Label.String = 'Normalized Power';
    cb_obj_M1_gc.Ticks = [-0.3:0.1:0.3];
%     cb_obj_gc.Label.String = 'Z-Score';
end
% 
if(~isempty(grand_average_M1.Right))
    a(8) = subplot(4,2,7);
    ax(8) = pcolor(1:100,log2(freq_vec),grand_average_M1.Right);
    ticks = logspace(log10(2.5),log10(50),5);
    ax(8).Parent.YTick = log2(ticks);
    ax(8).Parent.YTickLabel = sprintf('%1.2f\n',round(ticks,2));
    ylim([log2(2.5),log2(50)]);
    xticks([1,10,20,30,40,50,60,70,80,90,100]);
    xticklabels({'0','10','20','30','40','50','60','70','80','90','100'});
    shading interp;
    caxis(a(8),[-0.3,0.3]);
    ylabel({'Right';'Frequency (Hz)'});
    xlabel('% Gait Cycle');
end

%% Adjust single gait cycle plots
figure_format(18.5,9,4,'centimeters','painters');

% vSTN
curr_axes = subplot(4,2,1);
curr_axes.Position = [0.06,0.78,0.37,0.17];
curr_axes.Title.FontSize = 5;

curr_axes = subplot(4,2,3);
curr_axes.Position = [0.06,0.56,0.37,0.17];

% dSTN
curr_axes = subplot(4,2,2);
curr_axes.Position = [0.56,0.78,0.37,0.17];
curr_axes.Title.FontSize = 5;

curr_axes = subplot(4,2,4);
curr_axes.Position = [0.56,0.56,0.37,0.17];

% M1
curr_axes = subplot(4,2,5);
curr_axes.Position = [0.06,0.29,0.37,0.17];
curr_axes.Title.FontSize = 5;

curr_axes = subplot(4,2,7);
curr_axes.Position = [0.06,0.07,0.37,0.17];

% S1
curr_axes = subplot(4,2,6);
curr_axes.Position = [0.56,0.29,0.37,0.17];
curr_axes.Title.FontSize = 5;

curr_axes = subplot(4,2,8);
curr_axes.Position = [0.56,0.07,0.37,0.17];

% Colorbars
cb_obj_vSTN_gc.Position(2) = 0.56;
cb_obj_vSTN_gc.Position(3) = 0.01;
cb_obj_vSTN_gc.Position(4) = 0.39;
cb_obj_vSTN_gc.FontSize = 4;

cb_obj_dSTN_gc.Position(2) = 0.56;
cb_obj_dSTN_gc.Position(3) = 0.01;
cb_obj_dSTN_gc.Position(4) = 0.39;
cb_obj_dSTN_gc.FontSize = 4;

cb_obj_S1_gc.Position(2) = 0.07;
cb_obj_S1_gc.Position(3) = 0.01;
cb_obj_S1_gc.Position(4) = 0.39;
cb_obj_S1_gc.FontSize = 4;

cb_obj_M1_gc.Position(2) = 0.07;
cb_obj_M1_gc.Position(3) = 0.01;
cb_obj_M1_gc.Position(4) = 0.39;
cb_obj_M1_gc.FontSize = 4;

%% Add A, B annotations
% a_textbox = annotation('textbox',[0.01,0.98,0.1,0.03],'String','A','Fontsize',9,'EdgeColor','none');
% b_textbox = annotation('textbox',[0.01,0.48,0.1,0.03],'String','B','Fontsize',9,'EdgeColor','none');

%% Print
print('fig3_STN_M1_grand_average_gc_v7.pdf','-r300','-dpdf')
end