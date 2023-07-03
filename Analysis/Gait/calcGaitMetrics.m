function [out_table,cadence,gait_speed] = calcGaitMetrics(filename,varargin)
% Calculates spatial and temporal gait metrics from Delsys or Xsens data
%
% INPUTS:   gait_events         [=] Matrix of gait events. Must be 4
%                                   columns with the labels: LHS, RTO, RHS,
%                                   and LTO as column names in any order.
%
%           xsens_data          [=] Xsens table of data from csv file. 
%
%           xsens_filename      [=] Xsens csv file to read in.
% 
%           delsys_data         [=] Delsys struct of data.
%
%           delsys_filename     [=] Delsys .mat file to read in if delsys
%                                   data is not passed in.
%
%           level_type          [=] "none"/"single". Changes the way gait
%                                   events are automatically detected and
%                                   how gait metrics are calculated.
%                                   Default is "single".
%
%           flip_data           [=] true/false. Used to flip the 
%                                   coordinates to accurately calculate 
%                                   gait metrics. Should be set to true if 
%                                   on the xsens recording the mannequin 
%                                   body is walking into the screen from 
%                                   the origin. If walking to the left or 
%                                   right, should be set to false.
%           
%           med_state           [=] String input to denote medication state
%                                   of data. Can be on/off/low/na. Default
%                                   is na.
%
%           stim_state          [=] String input to denote stimulation
%                                   state. Can be clinical/aDBS_XXX/off/na.
%                                   If the option is aDBS_XXX, the XXX can
%                                   be changed to denote more information
%                                   about the aDBS settings. Default is na.
%
% OUTPUTS:  out_table           [=] Table of gait events by gait cycle.
%           cadence             [=] Cadence of patient during walking.
%           gait_speed          [=] Gait speed of patient during walking.
%
% EXAMPLE:
%           calcGaitMetrics('xsens_filename',filename2,'level_type','single','sample_tol',5,'preview',true);
%
% Author:   Kenneth Louie
% Date:     06/27/2023

% Example: 
% out_table = calcGaitMetrics('xsens_filename',xsens_filename,'level_type','single','flip_data',true,'baseline_nSamples',4500,'sample_tol',3,'preview',true);
% filename = '/Volumes/dwang3_shared/Patient Data/RC+S Data/RCS07/Gait/Pre-Processed Data/Xsens /Xsens Raw/MVNX/RCS07_Gait_12_17_19_walking_back_and_forth_after_task.csv'
% filename2 = '/Volumes/dwang3_shared/Patient Data/RC+S Data/RCS07/Gait/Pre-Processed Data/Xsens /Xsens Raw/MVNX/RCS07_Gait_12_17_19_walking_back_and_forth_before_task.csv'
%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:nargin/2
    switch varargin{i*2-1}
        case 'gait_events'
            gait_events = varargin{i*2};
        case 'xsens_data'
            xsens_data = varargin{i*2};
        case 'xsens_filename'
            xsens_filename = varargin{i*2};
        case 'delsys_data'
            delsys_data = varargin{i*2};
        case 'delsys_filename'
            delsys_filename = varargin{i*2};
        case 'level_type'
            level_type = varargin{i*2};
        case 'flip_data'
            % Used to flip the coordinate.
            flip_data = varargin{i*2};
        case 'med_state'
            med_state = varargin{i*2};
        case 'stim_state'
            stim_state = varargin{i*2};
        case 'subjectID'
            subjectID = varargin{i*2};
    end
    
end

%% Load data
if ~isempty(filename)
    load(filename);
    xsens_data = aligned_data.Xsens;
    delsys_data = aligned_data.Delsys;
    gait_events = aligned_data.gait_events;

    if isfield(aligned_data,'med_condition')
        med_state = aligned_data.med_condition;
    else
        med_state = 'na';
    end
    
    if isfield(aligned_data,'stim_condition')
        stim_state = aligned_data.stim_condition;
    else
        stim_state = 'na';
    end
end

if (~exist('xsens_data','var') || isempty(xsens_data)) && exist('xsens_filename','var')
    try
        xsens_data = readtable(xsens_filename);
    catch
        error('Xsens file does not exist, data is empty, or wrong filename.');
    end
end

if (~exist('delsys_data','var') || isempty(delsys_data)) && exist('delsys_filename','var')
    try
        temp_data = load(delsys_filename);
        delsys_data = temp_data.out_struct;
    catch
        error('Delsys file does not exist, data is empty, or wrong filename.');
    end
end

%% Set warnings if xsens dataset is missing
if ~exist('xsens_data','var') && exist('delsys_data','var')
    warning(sprintf('Gait Metric Calculation. \n No Xsens data detected. \n Delsys data detected \n Only temporal gait metrics will be calculated.'));
end

%% Default values
if ~exist('level_type','var') || isempty(level_type)
    level_type = 'single';
end

if ~exist('flip_data','var') || isempty(flip_data)
    flip_data = false;
end

if ~exist('med_state','var') || isempty(med_state)
    med_state = 'na';
end

if ~exist('stim_state','var') || isempty(stim_state)
    stim_state = 'na';
end

if ~exist('subjectID','var') || isempty(subjectID)
    subjectID = 'nan';
end

%% Calculate gait matrics and build output table
% out_table = array2table(nan(n_gait_events*2,11));
% out_table.Properties.VariableNames = {'Gait_Cycle','Side',...
%     'Step_Time','Step_Length','Step_Width',...
%     'Stride_Time','Stride_Length',...
%     'Swing_Time','Stance_Time',...
%     'Single_Limb_Support','Double_Limb_Support'};
% out_table.Gait_Cycle = repelem(1:n_gait_events,2)';
% out_table.Side = repmat(['R';'L'],n_gait_events,1);

if exist('xsens_data','var') && ~isempty(xsens_data)
    gait_metrics_table = getStepMetrics(gait_events,xsens_data,level_type);
elseif ~exist('xsens_data','var') && exist('delsys_data','var') && ~isempty(delsys_data)
    gait_metrics_table = getStepMetrics(gait_events,delsys_data,level_type);
end

SubjectID = repmat({subjectID},height(gait_metrics_table),1);
MedState = repmat({med_state},height(gait_metrics_table),1);
StimState = repmat({stim_state},height(gait_metrics_table),1);

gait_metrics_table = addvars(gait_metrics_table,SubjectID,MedState,StimState,'Before','Side');

cadence = [];
gait_speed = [];

end

function gait_metrics_table = getStepMetrics(gait_events,foot_data,level_type)
% [LHS,RTO,RHS,LTO]
% step length, width, time
% stride length, time
% swing time
% stance time
% double support time

% Order gait events
gait_events_turns_removed = removeGaitCyclesTurns(gait_events);
gait_events_ordered = sortGaitEvents(gait_events_turns_removed,'LHS');

if istable(foot_data)
    gc_start_search = find(gait_events_ordered.LHS>=foot_data.Time(1),1,'first');
    gc_end_search = find(gait_events_ordered.LHS<=foot_data.Time(end),1,'last');
    gait_events_ordered_trim = gait_events_ordered(gc_start_search:gc_end_search,:);
elseif isstruct(foot_data)
end

n_gait_events = size(gait_events_ordered_trim,1);
step_length = nan(n_gait_events*2,1);
step_time = nan(n_gait_events*2,1);
step_width = nan(n_gait_events*2,1);
stride_length = nan(n_gait_events*2,1);
stride_time = nan(n_gait_events*2,1);
swing_time = nan(n_gait_events*2,1);
stance_time = nan(n_gait_events*2,1);
dst = nan(n_gait_events*2,1);

if strcmp(level_type,'none')
    for i = 1:n_gait_events
        % i*2 = left metrics
        % i*2-1 = right metrics
        
        % out of bounds
        if i+1 > n_gait_events
            step_time(i*2) = nan;
            swing_time(i*2) = nan;
            sls(i*2) = nan;
            stance_time(i*2-1) = nan;
        else
            step_time(i*2) = (gait_events(i+1,1)-gait_events(i,3))/time_correction;
            swing_time(i*2) = (gait_events(i+1,1)-gait_events(i,4))/time_correction;
            sls(i*2) = (gait_events(i+1,1)-gait_events(i,4))/time_correction;
            stance_time(i*2-1) = (gait_events(i+1,2)-gait_events(i,3))/time_correction;
            
            if isnan(gait_events(i+1,1))
                step_length(i*2) = nan;
                step_width(i*2) = nan;
            else
                step_length(i*2) = abs(foot_data(gait_events_spatial(i+1,1),2)-foot_data(gait_events_spatial(i+1,1),5));
                step_width(i*2) = abs(foot_data(gait_events_spatial(i+1,1),1)-foot_data(gait_events_spatial(i+1,1),4));
            end
        end
        
        % Can't index with nan
        if isnan(gait_events(i,3))
            step_length(i*2-1) = nan;
            step_width(i*2-1) = nan;
        else
            step_length(i*2-1) = abs(foot_data(gait_events_spatial(i,3),5)-foot_data(gait_events_spatial(i,3),2));  % Right
            step_width(i*2-1) = abs(foot_data(gait_events_spatial(i,3),4)-foot_data(gait_events_spatial(i,3),1));   % Right
        end
        
        % Only looks at current gait cycle and doesn't matter if nan
        stance_time(i*2) = (gait_events(i,4)-gait_events(i,1))/time_correction;
        dls(i*2) = (gait_events(i,4)-gait_events(i,3))/time_correction;
        step_time(i*2-1) = (gait_events(i,3)-gait_events(i,1))/time_correction;
        swing_time(i*2-1) = (gait_events(i,3)-gait_events(i,2))/time_correction;
        sls(i*2-1) = (gait_events(i,3)-gait_events(i,2))/time_correction;
        dls(i*2-1) = (gait_events(i,2)-gait_events(i,1))/time_correction;
    end
end

if strcmp(level_type,'single')
    % [LHS,RTO,RHS,LTO]

    % Left events
    step_length(1:2:(n_gait_events-1)*2) = nan;   %TODO
    step_time(1:2:(n_gait_events-1)*2) = gait_events_ordered_trim.LHS(2:end)-gait_events_ordered_trim.RHS(1:end-1);     % Left step time = time from right heel strike to left heel-strike
    step_width(1:2:(n_gait_events-1)*2) = nan;    %TODO
    stride_length(1:2:(n_gait_events-1)*2) = nan; %TODO
    stride_time(1:2:(n_gait_events-1)*2) = gait_events_ordered_trim.LHS(2:end)-gait_events_ordered_trim.LHS(1:end-1);   % Left stride time = time from left heel strike to next left heel-strike
    swing_time(1:2:(n_gait_events-1)*2) = gait_events_ordered_trim.LHS(2:end)-gait_events_ordered_trim.LTO(1:end-1);    % Left swing time = Left toe-off to left heel-strike
    stance_time(1:2:n_gait_events*2) = gait_events_ordered_trim.LTO-gait_events_ordered_trim.LHS;                   % Left stance time = left heel-strike to left toe-off

    % Right events
    step_length(2:2:(n_gait_events-1)*2) = nan;     %TODO
    step_time(2:2:n_gait_events*2) = gait_events_ordered_trim.RHS-gait_events_ordered_trim.LHS;
    step_width(2:2:(n_gait_events-1)*2)= nan;       %TODO
    stride_length(2:2:(n_gait_events-1)*2) = nan;   %TODO
    stride_time(2:2:(n_gait_events-1)*2) = gait_events_ordered_trim.RHS(2:end)-gait_events_ordered_trim.RHS(1:end-1);
    swing_time(2:2:n_gait_events*2) = gait_events_ordered_trim.RHS-gait_events_ordered_trim.RTO;
    stance_time(2:2:(n_gait_events-1)*2) = gait_events_ordered_trim.RTO(2:end)-gait_events_ordered_trim.RHS(1:end-1);

    % Double support time
    dst(1:2:n_gait_events*2) = gait_events_ordered_trim.RTO-gait_events_ordered_trim.LHS;
    dst(2:2:n_gait_events*2) = gait_events_ordered_trim.LTO-gait_events_ordered_trim.RHS;
end

side = repmat({'L';'R'},n_gait_events,1);

% Create table
gait_metrics_table = table(side,step_length,step_time,step_width,stride_length,stride_time,swing_time,stance_time,dst,'VariableNames',{'Side','StepLength','StepTime','StepWidth','StrideLength','StrideTime','SwingTime','StanceTime','DoubleSupportTime'});
end

function gait_events_turns_removed = removeGaitCyclesTurns(gait_events)
gait_events_turns_removed = gait_events;
end