function gaitEventTable = getGaitEvents(data,varargin)

% Parse optional inputs
for i = 1:2:nargin-1
    switch varargin{i}
        case 'threshold'
            detectionThreshold = varargin{i+1};
        case 'markingTimeRange'
            markingTimeRange = varargin{i+1};
        case 'nextGaitCycleThreshold'
            nextGaitCycleThreshold = varargin{i+1};
        case 'firstMeta'
            firstMeta = varargin{i+1};
        case 'fifthMeta'
            fifthMeta = varargin{i+1};
    end
end

% Default values
if ~exist('detectionThreshold','var') || isempty(detectionThreshold)
    detectionThreshold = 5;
end

if ~exist('markingTimeRange','var') || isempty(markingTimeRange)
    markingTimeRange = [0,inf];
end

if ~exist('nextGaitCycleThreshold','var') || isempty(nextGaitCycleThreshold)
    nextGaitCycleThreshold = 1.3;
end

if ~exist('firstMeta','var') || isempty(firstMeta)
    firstMeta = false;
end

if ~exist('fifthMeta','var') || isempty(fifthMeta)
    fifthMeta = false;
end

%% Determine gait events
% Heel strike
if sum(~cellfun(@isempty,regexp(data.Chan_names,'FSR.*Analog'))) > 0
    LHS_chan_name = 'FSR_adapter_Left_Analog_1';
    LHS_chan_name2 = 'FSR_adapter_Left_Analog_3';
    RHS_chan_name = 'FSR_adapter_Right_Analog_1';
    RHS_chan_name2 = 'FSR_adapter_Right_Analog_3';
else
    LHS_chan_name = 'FSR_adapter_15_Left_FSRA_15';
    LHS_chan_name2 = 'FSR_adapter_15_Left_FSRC_15';
    RHS_chan_name = 'FSR_adapter_16_Right_FSRA_16';
    RHS_chan_name2 = 'FSR_adapter_16_Right_FSRC_16';
end

LHS_times = data.Time.(LHS_chan_name)((data.Data.(LHS_chan_name)(1:end-1) < detectionThreshold) & (data.Data.(LHS_chan_name)(2:end) > detectionThreshold))';
LHS_times = LHS_times(LHS_times>=markingTimeRange(1) & LHS_times<=markingTimeRange(2));
LHS_source = zeros([length(LHS_times),1]);
remove_inds = find(diff(LHS_times)<0.75)+1;
LHS_times(remove_inds) = [];
LHS_source(remove_inds) = [];

RHS_times = data.Time.(RHS_chan_name)((data.Data.(RHS_chan_name)(1:end-1) < detectionThreshold) & (data.Data.(RHS_chan_name)(2:end) > detectionThreshold))';
RHS_times = RHS_times(RHS_times>=markingTimeRange(1) & RHS_times<=markingTimeRange(2));
RHS_source = zeros([length(RHS_times),1]);
remove_inds = find(diff(RHS_times)<0.75)+1;
RHS_times(remove_inds) = [];
RHS_source(remove_inds) = [];

if fifthMeta
    LHS_times2 = data.Time.(LHS_chan_name2)((data.Data.(LHS_chan_name2)(1:end-1) < detectionThreshold) & (data.Data.(LHS_chan_name2)(2:end) > detectionThreshold))';
    LHS_times2 = LHS_times2(LHS_times2>=markingTimeRange(1) & LHS_times2<=markingTimeRange(2));
    LHS_source2 = ones([length(LHS_times2),1]);
    remove_inds = find(diff(LHS_times2)<0.75)+1;
    LHS_times2(remove_inds) = [];
    LHS_source2(remove_inds) = [];
    LHS_times = combineAndSort([LHS_times;LHS_times2],[LHS_source;LHS_source2]);

    RHS_times2 = data.Time.(RHS_chan_name2)((data.Data.(RHS_chan_name2)(1:end-1) < detectionThreshold) & (data.Data.(RHS_chan_name2)(2:end) > detectionThreshold))';
    RHS_times2 = RHS_times2(RHS_times2>=markingTimeRange(1) & RHS_times2<=markingTimeRange(2));
    RHS_source2 = ones([length(RHS_times2),1]);
    remove_inds = find(diff(RHS_times2)<0.75)+1;
    RHS_times2(remove_inds) = [];
    RHS_source2(remove_inds) = [];
    RHS_times = combineAndSort([RHS_times;RHS_times2],[RHS_source;RHS_source2]);
end

% Toe off
if sum(~cellfun(@isempty,regexp(data.Chan_names,'FSR.*Analog'))) > 0
    LTO_chan_name = 'FSR_adapter_Left_Analog_4';
    LTO_chan_name2 = 'FSR_adapter_Left_Analog_2';
    RTO_chan_name = 'FSR_adapter_Right_Analog_4';
    RTO_chan_name2 = 'FSR_adapter_Right_Analog_2';
else
    LTO_chan_name = 'FSR_adapter_15_Left_FSRD_15';
    LTO_chan_name2 = 'FSR_adapter_15_Left_FSRB_15';
    RTO_chan_name = 'FSR_adapter_16_Right_FSRD_16';
    RTO_chan_name2 = 'FSR_adapter_16_Right_FSRB_16';
end

LTO_times = data.Time.(LTO_chan_name)((data.Data.(LTO_chan_name)(1:end-1) > detectionThreshold) & (data.Data.(LTO_chan_name)(2:end) < detectionThreshold))';
LTO_times = LTO_times(LTO_times>=markingTimeRange(1) & LTO_times<=markingTimeRange(2));
LTO_source = zeros([length(LTO_times),1]);
remove_inds = find(diff(LTO_times)<0.75)+1;
LTO_times(remove_inds) = [];
LTO_source(remove_inds) = [];

RTO_times = data.Time.(RTO_chan_name)((data.Data.(RTO_chan_name)(1:end-1) > detectionThreshold) & (data.Data.(RTO_chan_name)(2:end) < detectionThreshold))';
RTO_times = RTO_times(RTO_times>=markingTimeRange(1) & RTO_times<=markingTimeRange(2));
RTO_source = zeros([length(RTO_times),1]);
remove_inds = find(diff(RTO_times)<0.75)+1;
RTO_times(remove_inds) = [];
RTO_source(remove_inds) = [];

if firstMeta
    LTO_times2 = data.Time.(LTO_chan_name2)((data.Data.(LTO_chan_name2)(1:end-1) > detectionThreshold) & (data.Data.(LTO_chan_name2)(2:end) < detectionThreshold))';
    LTO_times2 = LTO_times2(LTO_times2>=markingTimeRange(1) & LTO_times2<=markingTimeRange(2));
    LTO_source2 = ones([length(LTO_times2),1]);
    remove_inds = find(diff(LTO_times2)<0.75)+1;
    LTO_times2(remove_inds) = [];
    LTO_source2(remove_inds) = [];
    LTO_times = combineAndSort([LTO_times';LTO_times2'],[LTO_source;LTO_source2]);

    RTO_times2 = data.Time.(RTO_chan_name2)((data.Data.(RTO_chan_name2)(1:end-1) > detectionThreshold) & (data.Data.(RTO_chan_name2)(2:end) < detectionThreshold))';
    RTO_times2 = RTO_times2(RTO_times2>=markingTimeRange(1) & RTO_times2<=markingTimeRange(2));
    RTO_source2 = ones([length(RTO_times2),1]);
    remove_inds = find(diff(RTO_times2)<0.75)+1;
    RTO_times2(remove_inds) = [];
    RTO_source2(remove_inds) = [];
    RTO_times = combineAndSort([RTO_times';RTO_times2'],[RTO_source;RTO_source2]);
end

% % Remove gait events that are before marking start time
% LHS_times = LHS_times(LHS_times>=markingTimeRange(1) & LHS_times<=markingTimeRange(2));
% RTO_times = RTO_times(RTO_times>=markingTimeRange(1) & RTO_times<=markingTimeRange(2));
% RHS_times = RHS_times(RHS_times>=markingTimeRange(1) & RHS_times<=markingTimeRange(2));
% LTO_times = LTO_times(LTO_times>=markingTimeRange(1) & LTO_times<=markingTimeRange(2));
% 
% % Remove gait events that are too close in time to each other
% LHS_times(find(diff(LHS_times)<0.75)+1) = [];
% RTO_times(find(diff(RTO_times)<0.75)+1) = [];
% RHS_times(find(diff(RHS_times)<0.75)+1) = [];
% LTO_times(find(diff(LTO_times)<0.75)+1) = [];

% Setup output variable
gaitEventMat = [];
gaitEventMat_count = 0;

% State variables
LHS_count = 1;
RTO_count = 1;
RHS_count = 1;
LTO_count = 1;

% Determine start event
try
    [~,startEventInd] = min([LHS_times(1),RTO_times(1),RHS_times(1),LTO_times(1)]);
catch
    A = zeros(4,1);
    if isempty(LHS_times)
        A(1) = nan;
    else
        A(1) = LHS_times(1);
    end

    if isempty(RTO_times)
        A(2) = nan;
    else
        A(2) = RTO_times(1);
    end

    if isempty(RHS_times)
        A(3) = nan;
    else
        A(3) = RHS_times(1);
    end

    if isempty(LTO_times)
        A(4) = nan;
    else
        A(4) = LTO_times(1);
    end

    [~,startEventInd] = min(A,[],'omitnan');

end
if startEventInd == 1
    startEvent = 'LHS';
elseif startEventInd == 2
    startEvent = 'RTO';
elseif startEventInd == 3
    startEvent = 'RHS';
else
    startEvent = 'LTO';
end
storingIndexFun = @(x) (x-startEventInd)+1;
referenceIndexFun = @(x) ceil(mod(x,4.1));

% Go through all events detected and order them
while (LHS_count + RTO_count + RHS_count + LTO_count) ~= sum([length(LHS_times),length(RTO_times),length(RHS_times),length(LTO_times)])+4
    gaitEventMat = [gaitEventMat;nan(1,4)];
    gaitEventMat_count = gaitEventMat_count + 1;
    
    if length(LHS_times)>=LHS_count && length(RTO_times) >= RTO_count && length(RHS_times) >= RHS_count && length(LTO_times) >= LTO_count
        currGaitEvents = [LHS_times(LHS_count),RTO_times(RTO_count),RHS_times(RHS_count),LTO_times(LTO_count)];
    else
        currGaitEvents = nan(1,4);
        if length(LHS_times)>=LHS_count
            currGaitEvents(1) = LHS_times(LHS_count);
        end
        if length(RTO_times)>=RTO_count
            currGaitEvents(2) = RTO_times(RTO_count);
        end
        if length(RHS_times)>=RHS_count
            currGaitEvents(3) = RHS_times(RHS_count);
        end
        if length(LTO_times)>=LTO_count
            currGaitEvents(4) = LTO_times(LTO_count);
        end
    end
    [minVal,minInd] = min(currGaitEvents);
    
    for i = minInd:minInd+3
        if sum(isnan(gaitEventMat(gaitEventMat_count,:))) ~= 4
            storingInd = storingIndexFun(i);
            if storingInd > 4
                storingInd = storingInd-4;
            elseif storingInd <= 0
                storingInd = storingInd+4;
            end
            if abs(currGaitEvents(referenceIndexFun(i))-min(gaitEventMat(gaitEventMat_count,:))) < nextGaitCycleThreshold && isnan(gaitEventMat(gaitEventMat_count,storingInd))
                gaitEventMat(gaitEventMat_count,storingInd) = currGaitEvents(referenceIndexFun(i));
                switch referenceIndexFun(i)
                    case 1
                        LHS_count = LHS_count+1;
                    case 2
                        RTO_count = RTO_count+1;
                    case 3
                        RHS_count = RHS_count+1;
                    case 4
                        LTO_count = LTO_count+1;
                end
            end
        else
            storingInd = storingIndexFun(minInd);
            if storingInd > 4
                storingInd = storingInd-4;
            elseif storingInd <= 0
                storingInd = storingInd+4;
            end
            gaitEventMat(gaitEventMat_count,storingInd) = minVal;
            switch referenceIndexFun(minInd)
                case 1
                    LHS_count = LHS_count+1;
                case 2
                    RTO_count = RTO_count+1;
                case 3
                    RHS_count = RHS_count+1;
                case 4
                    LTO_count = LTO_count+1;
            end
        end
    end
    
    if gaitEventMat_count>1 && any(diff(gaitEventMat(gaitEventMat_count,:))<0)
        neg_diff_inds = find(diff(gaitEventMat(gaitEventMat_count,:))<0);
        gaitEventMat = [gaitEventMat;nan(1,4)];
        for i = 1:length(neg_diff_inds)
            gaitEventMat(gaitEventMat_count+1,neg_diff_inds(i)) = gaitEventMat(gaitEventMat_count,neg_diff_inds(i));
            gaitEventMat(gaitEventMat_count,neg_diff_inds(i)) = nan;
        end
    end
end

% Get rid of rows with all nan
removeInd = [];
for i = 1:size(gaitEventMat,1)
    if sum(isnan(gaitEventMat(i,:))) == 4
        removeInd = [removeInd;i];
    end
end
gaitEventMat(removeInd,:) = [];

gaitEventOrder = {startEvent};
for i = 1:3
    gaitEventOrder{i+1} = getNextGaitEventName(gaitEventOrder{i});
end

gaitEventTable = array2table(gaitEventMat,'VariableNames',gaitEventOrder);

end

function nextGaitEventName = getNextGaitEventName(currEvent)

if strcmp(currEvent,'LHS')
    nextGaitEventName = 'RTO';
elseif strcmp(currEvent,'RTO')
    nextGaitEventName = 'RHS';
elseif strcmp(currEvent,'RHS')
    nextGaitEventName = 'LTO';
else
    nextGaitEventName = 'LHS';
end
end

function outGaitEvents = combineAndSort(gait_events,gait_event_source)
    [gait_events_sorted,sorted_inds] = sort(gait_events);
    gait_event_source_sorted = gait_event_source(sorted_inds);

    remove_inds = [];
    for i = 1:length(gait_events_sorted)-1
        if (gait_events_sorted(i+1)-gait_events_sorted(i)) < 0.2
            if gait_event_source_sorted(i+1) == 0 && gait_event_source_sorted(i) == 1
                remove_inds(end+1) = i;
            elseif gait_event_source_sorted(i+1) == 1 && gait_event_source_sorted(i) == 0
                remove_inds(end+1) = i+1;
            end
        end
    end

    gait_events_sorted(remove_inds) = [];
    outGaitEvents = gait_events_sorted;
end