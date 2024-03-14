function gaitEventTable = getGaitEventsGSLT(data,varargin)

% Parse optional inputs
for i = 1:2:nargin-1
    switch varargin{i}
        case 'detectionThreshold'
            detectionThreshold = varargin{i+1};
        case 'markingTimeRange'
            markingTimeRange = varargin{i+1};
        case 'startTrimTime'
            start_trim_time = varargin{i+1};
    end
end

% Load data
if ~exist('data','var') || isempty(data)
    [filename,filepath]=uigetfile('*.mat');
    load(fullfile(filepath,filename));
    data = out_struct;
end

% Default values
if ~exist('detectionThreshold','var') || isempty(detectionThreshold)
    detectionThreshold = 5;
end

if ~exist('markingTimeRange','var') || isempty(markingTimeRange)
    markingTimeRange = [0,inf];
end

if ~exist('start_trim_time','var') || ~isnumeric(start_trim_time)|| isempty(start_trim_time)
    start_trim_time = -10;
end

% Determine gait events
% LHS_times = data.Time.FSR_adapter_15_Left_FSRA_15((data.Data.FSR_adapter_15_Left_FSRA_15(1:end-1) < detectionThreshold) & (data.Data.FSR_adapter_15_Left_FSRA_15(2:end) > detectionThreshold))';
LHS_times = data.Time.FSR_adapter_15_Left_FSRA_15((data.Data.FSR_adapter_15_Left_FSRA_15(1:end-1) < 5) & (data.Data.FSR_adapter_15_Left_FSRA_15(2:end) > 5))';

% RTO_times = data.Time.FSR_adapter_16_Right_FSRD_16((data.Data.FSR_adapter_16_Right_FSRD_16(1:end-1) > detectionThreshold) & (data.Data.FSR_adapter_16_Right_FSRD_16(2:end) < detectionThreshold))';
RTO_times = data.Time.FSR_adapter_16_Right_FSRD_16((data.Data.FSR_adapter_16_Right_FSRD_16(1:end-1) > 1) & (data.Data.FSR_adapter_16_Right_FSRD_16(2:end) < 1))';

% RHS_times = data.Time.FSR_adapter_16_Right_FSRA_16((data.Data.FSR_adapter_16_Right_FSRA_16(1:end-1) < detectionThreshold) & (data.Data.FSR_adapter_16_Right_FSRA_16(2:end) > detectionThreshold))';
RHS_times = data.Time.FSR_adapter_16_Right_FSRA_16((data.Data.FSR_adapter_16_Right_FSRA_16(1:end-1) < -0.01) & (data.Data.FSR_adapter_16_Right_FSRA_16(2:end) > -0.01))';

% LTO_times = data.Time.FSR_adapter_15_Left_FSRD_15((data.Data.FSR_adapter_15_Left_FSRD_15(1:end-1) > detectionThreshold) & (data.Data.FSR_adapter_15_Left_FSRD_15(2:end) < detectionThreshold))';
LTO_times = data.Time.FSR_adapter_15_Left_FSRD_15((data.Data.FSR_adapter_15_Left_FSRD_15(1:end-1) > 5) & (data.Data.FSR_adapter_15_Left_FSRD_15(2:end) < 5))';


% Remove gait events that are before marking start time
LHS_times = LHS_times(LHS_times>=markingTimeRange(1) & LHS_times<=markingTimeRange(2));
RTO_times = RTO_times(RTO_times>=markingTimeRange(1) & RTO_times<=markingTimeRange(2));
RHS_times = RHS_times(RHS_times>=markingTimeRange(1) & RHS_times<=markingTimeRange(2));
LTO_times = LTO_times(LTO_times>=markingTimeRange(1) & LTO_times<=markingTimeRange(2));

% Remove gait events that are too close in time to each other
LHS_times(find(diff(LHS_times)<0.75)+1) = [];
RTO_times(find(diff(RTO_times)<0.75)+1) = [];
RHS_times(find(diff(RHS_times)<0.75)+1) = [];
LTO_times(find(diff(LTO_times)<0.75)+1) = [];

% Setup output variable
gaitEventMat = [];
gaitEventMat_count = 0;

% State variables
LHS_count = 1;
RTO_count = 1;
RHS_count = 1;
LTO_count = 1;

% Determine start event
[~,startEventInd] = min([LHS_times(1),RTO_times(1),RHS_times(1),LTO_times(1)]);
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
            if abs(currGaitEvents(referenceIndexFun(i))-min(gaitEventMat(gaitEventMat_count,:))) < 1.5 && isnan(gaitEventMat(gaitEventMat_count,storingInd))
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

% Add 
gaitEventTable = addGSLTinfo(gaitEventTable,start_trim_time);

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