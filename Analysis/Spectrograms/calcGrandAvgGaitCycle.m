function calcGrandAvgGaitCycle(fileList,varargin)
%% calcGrandAvgGaitCycle
% Calculates the grand average gait cycle average spectrogram from all the
% files passed into this function. Multiple optional inputs allows
% customization of the computation including resolution of the gait cycle,
% and normalization strategy.
%
% INPUTS:  Required
%               fileList    [=] list of strings in a cell array. Each
%                               string is the file path to the .mat file of
%                               the aligned data with gait events.
%          Optional
%               nPercentBins[=] The number of bins to break up the gait
%                               cycle in. Default 100, meaning each bin is
%                               1% of the gait cycle.
%
%               gcStartEvent[=] Which gait event to be the start of the
%                               gait cycle. Can be left heel strike (LHS),
%                               right heel strike (RHS), left toe off
%                               (LTO), right toe off (RTO). Default is LHS.

%               normBy      [=] What to normalize the gait cycle average
%                               too. Can be 3 options: none, baseline,
%                               average_during_walking. Default is none.

%               normType    [=] How to calculate the normalization. Can be
%                               3 options: none, percent_change, or zscore.
%                               Default is none.

%               geRangeTable[=] Limits the number of gait events to
%                               consider for the gait cycle average. This
%                               is because the gait event markings may
%                               include when the patient is walking to a
%                               chair at the end of the trial. Default is
%                               to consider all gait events.
%               
%               keys        [=] Cell array of keys to analyze. Default is
%                               to analyze all keys.
%
%               swapKeys    [=] 2x<number of files> cell array that
%                               indicates if keys need to be swapped.
%                               This is so the code can be backward
%                               compatible with Phil's RCS patient's
%                               recording settings.
%
%               savePlot    [=] Boolean option to save the resulting plot.
%                               Default is false.
%
%   Example call:
%
%
% Date:     05/16/2022
% Author:   Kenneth H. Louie (kenneth.louie@ucsf.edu)
% Project:  MJFF aDBS Gait

%% Option variables
for i = 1:2:nargin-2
    switch varargin{i}
        case 'nPercentBins'
            nPercentBins = varargin{i+1};
        case 'gcStartEvent'
            gcStartEvent = varargin{i+1};
        case 'normBy'
            normBy = varargin{i+1};
        case 'normType'
            normType = varargin{i+1};
        case 'geRangeTable'
            geRangeTable = varargin{i+1};
        case 'keys'
            keys = varargin{i+1};
        case 'swapKeys'
            swapKeys = varargin{i+1};
        case 'savePlot'
            savePlot = varargin{i+1};
    end
end

% Set default options if not passed in by user
if ~exist('nPercentBins','var') || isempty(nPercentBins)
    nPercentBins = 100;
end

if ~exist('gcStartEvent','var') || isempty(gcStartEvent)
    gcStartEvent = 'LHS';
end

if ~exist('normBy','var') || isempty(normBy)
    normBy = 'none';
end

if ~exist('normType','var') || isempty(normType)
    normType = 'none';
end

if ~exist('geRangeTable','var') || isempty(geRangeTable)
    geRangeTable = [];
end

if ~exist('keys','var') || isempty(keys)
    keys = {'key0','key1','key2','key3'};
end

if ~exist('swapKeys','var') || isempty('swapKeys')
    swapKeys = [];
end

if ~exist('savePlot','var') || isempty(savePlot)
    savePlot = 0;   % Does not save plot by default.
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
normalization.Left.key1 = cell(1,length(fileList));
normalization.Left.key2 = cell(1,length(fileList));
normalization.Left.key3 = cell(1,length(fileList));

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

% Hold normalized data
normalized.Left.key0 = [];
normalized.Left.key1 = [];
normalized.Left.key2 = [];
normalized.Left.key3 = [];

normalized.Right.key0 = [];
normalized.Right.key1 = [];
normalized.Right.key2 = [];
normalized.Right.key3 = [];

% Hold grand averages
grandAverage.Left.key0 = [];
grandAverage.Left.key1 = [];
grandAverage.Left.key2 = [];
grandAverage.Left.key3 = [];

grandAverage.Right.key0 = [];
grandAverage.Right.key1 = [];
grandAverage.Right.key2 = [];
grandAverage.Right.key3 = [];

leftImplant = zeros(1,length(fileList));
rightImplant = zeros(1,length(fileList));

%% Go through all files and extract data
for i = 1:length(fileList)
    load(fileList{i});
    gaitEventsSorted = sortGaitEvents(aligned_data.gait_events,gcStartEvent);
    signalAnalysisData = calcRCS_CWT(aligned_data);
    
    % Check if neural data keys need to be swapped
    if ~isempty(swapKeys) && ~isempty(swapKeys{i})
        switch swapKeys{i}{1}
            case 'key0'
                ind1 = 1;
            case 'key1'
                ind1 = 2;
            case 'key2'
                ind1 = 3;
            case 'key3'
                ind1 = 4;
        end

        switch swapKeys{i}{2}
            case 'key0'
                ind2 = 1;
            case 'key1'
                ind2 = 2;
            case 'key2'
                ind2 = 3;
            case 'key3'
                ind2 = 4;
        end

        signalAnalysisData = swapData(signalAnalysisData,ind1,ind2);
    end
    
    if isempty(geRangeTable)
        geRange = [1,height(gaitEventsSorted)];
    elseif isinf(geRangeTable(i,2)) && geRangeTable(i,1) == 1
        geRange = [1,height(gaitEventsSorted)];
    else
        start_ind = find(gaitEventsSorted.(gcStartEvent) > geRangeTable(i,1),1,'first');
        end_ind = find(gaitEventsSorted.(gcStartEvent) < geRangeTable(i,2),1,'last');
        geRange = [start_ind,end_ind];
    end
    
    if isfield(signalAnalysisData,'Left')
        leftImplant(i) = 1;
        nChannels = length(signalAnalysisData.Left.Time);
        for j = 1:nChannels
            if sum(contains(keys,sprintf('key%i',j-1))) ~= 0
                walking_start_ind = find(signalAnalysisData.Left.Time{j} >= min(gaitEventsSorted{geRange(1),:})-1,1,'first');
                walking_end_ind = find(signalAnalysisData.Left.Time{j} <= max(gaitEventsSorted{geRange(2),:}),1,'last');

                gait_cycle_mat_left = zeros(length(signalAnalysisData.Left.Freq_Values{j}),nPercentBins,1);
                count = 1;
                for k = geRange(1):geRange(2)-1
                    if ~isnan(gaitEventsSorted.(gcStartEvent)(k)) && ~isnan(gaitEventsSorted.(gcStartEvent)(k+1)) && (diff(gaitEventsSorted.(gcStartEvent)([k,k+1])) < 2)
                        [~,start_ind] = min(abs(signalAnalysisData.Left.Time{j}-gaitEventsSorted.(gcStartEvent)(k)));
                        [~,end_ind] = min(abs(signalAnalysisData.Left.Time{j}-gaitEventsSorted.(gcStartEvent)(k+1)));
                        data_snip = abs(signalAnalysisData.Left.Values{j}(:,start_ind:end_ind));

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

                gc.Left.(['key',num2str(j-1)]){i} = gait_cycle_mat_left;
                if strcmp(normBy,'average_during_walking')
                    normalization.Left.(['key',num2str(j-1)]){i} = [mean(abs(signalAnalysisData.Left.Values{j}(:,walking_start_ind:walking_end_ind)),2,'omitnan'),...
                        std(abs(signalAnalysisData.Left.Values{j}(:,walking_start_ind:walking_end_ind)),0,2,'omitnan'),...
                        median(abs(signalAnalysisData.Left.Values{j}(:,walking_start_ind:walking_end_ind)),2,'omitnan')];
                elseif strcmp(normBy,'baseline')
                    % TODO
                end
                fc.Left.(['key',num2str(j-1)]){i} = signalAnalysisData.Left.Freq_Values{j};
            end
        end
    end
    
    if isfield(signalAnalysisData,'Right')
        rightImplant(i) = 1;
        nChannels = length(signalAnalysisData.Right.Time);
        for j = 1:nChannels
            if sum(contains(keys,sprintf('key%i',j-1))) ~= 0
                walking_start_ind = find(signalAnalysisData.Right.Time{j} >= min(gaitEventsSorted{geRange(1),:})-1,1,'first');
                walking_end_ind = find(signalAnalysisData.Right.Time{j} <= max(gaitEventsSorted{geRange(2),:}),1,'last');

                gait_cycle_mat_right = zeros(length(signalAnalysisData.Right.Freq_Values{j}),nPercentBins,1);
                count = 1;
                for k = geRange(1):geRange(2)-1
                    if ~isnan(gaitEventsSorted.(gcStartEvent)(k)) && ~isnan(gaitEventsSorted.(gcStartEvent)(k+1)) && (diff(gaitEventsSorted.(gcStartEvent)([k,k+1])) < 2)
                        [~,start_ind] = min(abs(signalAnalysisData.Right.Time{j}-gaitEventsSorted.(gcStartEvent)(k)));
                        [~,end_ind] = min(abs(signalAnalysisData.Right.Time{j}-gaitEventsSorted.(gcStartEvent)(k+1)));
                        data_snip = abs(signalAnalysisData.Right.Values{j}(:,start_ind:end_ind));

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

                gc.Right.(['key',num2str(j-1)]){i} = gait_cycle_mat_right;
                if strcmp(normBy,'average_during_walking')
                    normalization.Right.(['key',num2str(j-1)]){i} = [mean(abs(signalAnalysisData.Right.Values{j}(:,walking_start_ind:walking_end_ind)),2,'omitnan'),...
                        std(abs(signalAnalysisData.Right.Values{j}(:,walking_start_ind:walking_end_ind)),0,2,'omitnan'),...
                        median(abs(signalAnalysisData.Right.Values{j}(:,walking_start_ind:walking_end_ind)),2,'omitnan')];
                elseif strcmp(normBy,'baseline')
                    % TODO
                end
                fc.Right.(['key',num2str(j-1)]){i} = signalAnalysisData.Right.Freq_Values{j};
            end
        end
    end
end

% Remove fields of non analyzed keys
allKeys = {'key0','key1','key2','key3'};
removeInds = ~contains(allKeys,keys);
keys2Remove = allKeys(removeInds);

for i = 1:length(keys2Remove)
    gc.Left = rmfield(gc.Left,keys2Remove{i});
    gc.Right = rmfield(gc.Right,keys2Remove{i});
    normalization.Left = rmfield(normalization.Left,keys2Remove{i});
    normalization.Right = rmfield(normalization.Right,keys2Remove{i});
    fc.Left = rmfield(fc.Left,keys2Remove{i});
    fc.Right = rmfield(fc.Right,keys2Remove{i});
    normalized.Left = rmfield(normalized.Left,keys2Remove{i});
    normalized.Right = rmfield(normalized.Right,keys2Remove{i});
    grandAverage.Left = rmfield(grandAverage.Left,keys2Remove{i});
    grandAverage.Right = rmfield(grandAverage.Right,keys2Remove{i});
end


% Normalize and average all together
for i = 1:length(fields(gc.Left))
    for j = 1:length(gc.Left.(['key',num2str(i-1)]))
        if ~isempty(gc.Left.(['key',num2str(i-1)]){j})
            start_ind = find(fc.Left.(['key',num2str(i-1)]){j}<=100,1,'first');
            end_ind = find(fc.Left.(['key',num2str(i-1)]){j}>=0.1,1,'last');
            normalization_mat = nan(size(gc.Left.(['key',num2str(i-1)]){j}));
            for k = 1:size(gc.Left.(['key',num2str(i-1)]){j},3)
                if strcmp(normType,'zscore')
                    normalization_mat(:,:,k) = normalizeData(gc.Left.(['key',num2str(i-1)]){j}(:,:,k),'zscore',normalization.Left.(['key',num2str(i-1)]){j});
                elseif strcmp(normType,'percent_change')
                    % TODO
                end
            end
            normalized.Left.(['key',num2str(i-1)]) = cat(3,normalized.Left.(['key',num2str(i-1)]),normalization_mat(start_ind:end_ind,:,:));
        end
    end
end

for i = 1:length(fields(gc.Right))
    for j = 1:length(gc.Right.(['key',num2str(i-1)]))
        if ~isempty(gc.Right.(['key',num2str(i-1)]){j})
            start_ind = find(fc.Right.(['key',num2str(i-1)]){j}<=100,1,'first');
            end_ind = find(fc.Right.(['key',num2str(i-1)]){j}>=0.1,1,'last');
            normalization_mat = nan(size(gc.Right.(['key',num2str(i-1)]){j}));
            for k = 1:size(gc.Right.(['key',num2str(i-1)]){j},3)
                normalization_mat(:,:,k) = normalizeData(gc.Right.(['key',num2str(i-1)]){j}(:,:,k),'zscore',normalization.Right.(['key',num2str(i-1)]){j});
            end
            normalized.Right.(['key',num2str(i-1)]) = cat(3,normalized.Right.(['key',num2str(i-1)]),normalization_mat(start_ind:end_ind,:,:));
        end
    end
end

% calculate grand averages
for i = 1:length(fields(normalized.Left))
    grandAverage.Left.(['key',num2str(i-1)]) = mean(normalized.Left.(['key',num2str(i-1)]),3);
end

for i = 1:length(fields(normalized.Right))
    grandAverage.Right.(['key',num2str(i-1)]) = mean(normalized.Right.(['key',num2str(i-1)]),3);
end

%% Plotting
for i = 1:length(fields(grandAverage.Left))
    start_ind = find(fc.Left.(['key',num2str(i-1)]){1}<=100,1,'first');
    end_ind = find(fc.Left.(['key',num2str(i-1)]){1}>=0.1,1,'last');
    freq_vec = fc.Left.(['key',num2str(i-1)]){1}(start_ind:end_ind);

    figure;
    ax = pcolor(1:100,log2(freq_vec),grandAverage.Left.(['key',num2str(i-1)]));
    ticks = logspace(log10(2.5),log10(50),5);
    ax.Parent.YTick = log2(ticks);
    ax.Parent.YTickLabel = sprintf('%1.2f\n',round(ticks,2));
    ylim([log2(2.5),log2(50)]);
    xticks([1,10,20,30,40,50,60,70,80,90,100]);
    xticklabels({'0','10','20','30','40','50','60','70','80','90','100'});
    shading interp;
    ylabel('Frequency (Hz)');
    title({'Left';sprintf('key%i Grand Average',i-1)});
end

for i = 1:length(fields(grandAverage.Right))
    start_ind = find(fc.Right.(['key',num2str(i-1)]){1}<=100,1,'first');
    end_ind = find(fc.Right.(['key',num2str(i-1)]){1}>=0.1,1,'last');
    freq_vec = fc.Right.(['key',num2str(i-1)]){1}(start_ind:end_ind);

    figure;
    ax = pcolor(1:100,log2(freq_vec),grandAverage.Right.(['key',num2str(i-1)]));
    ticks = logspace(log10(2.5),log10(50),5);
    ax.Parent.YTick = log2(ticks);
    ax.Parent.YTickLabel = sprintf('%1.2f\n',round(ticks,2));
    ylim([log2(2.5),log2(50)]);
    xticks([1,10,20,30,40,50,60,70,80,90,100]);
    xticklabels({'0','10','20','30','40','50','60','70','80','90','100'});
    shading interp;
    ylabel('Frequency (Hz)');
    title({'Right';sprintf('key%i Grand Average',i-1)});
end
end

function newData = swapData(oldData,swapInd1,swapInd2)
newData = oldData;
topFields = fields(oldData);

for i = 1:length(topFields)
    subFields = fields(oldData.(topFields{i}));
    
    for j = 1:length(subFields)
        swapData = oldData.(topFields{i}).(subFields{j}){swapInd1};
        newData.(topFields{i}).(subFields{j}){swapInd1} = newData.(topFields{i}).(subFields{j}){swapInd2};
        newData.(topFields{i}).(subFields{j}){swapInd2} = swapData;
    end
end
end