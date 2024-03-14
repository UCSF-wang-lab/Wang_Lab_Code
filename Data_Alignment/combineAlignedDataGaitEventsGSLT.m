function combineAlignedDataGaitEventsGSLT(filename_aligned_data,filename_gait_events)
% Adds appropriately time-aligned gait events ('gait_event' field) to the
% 'aligned_data' structure. This field is also populated with information
% about which of the gait events correspond to successfully adapted steps 
% (hits) and which to unsuccessfully adapted steps (misses).
%
% Inputs:   filename_aligned_data   [=] File of the extracted target information
%           filename_gait_events    [=] Save name and path of the .mat file
%
% Outputs:  None
%
% Author:   Kenneth Louie
% Contributor: Eleni Patelaki
% Date:     3/13/24

% load Aligned data mat file
if ~exist('filename_aligned_data','var')
    [filename_aligned_data,path] = uigetfile('*.mat');
    load(fullfile(path,filename_aligned_data));
else
    load(filename_aligned_data);
    [path,temp1,temp2] = fileparts(filename_aligned_data);
    filename_aligned_data = [temp1,temp2];
end


% load Gait event file
if ~exist('filename_gait_events','var')
    [filename_gait_events,path2] = uigetfile('*.csv');
    if filename_gait_events == 0
        [filename_gait_events,path2] = uigetfile('*.txt');
    end
    fid = fopen(fullfile(path2,filename_gait_events));
else
    fid = fopen(filename_gait_events);
    [path2,filename_gait_events] = fileparts(filename_gait_events);
end

count = 0;
gait_event_line = true;
while gait_event_line
    curr_line = fgetl(fid);
    count = count + 1;
    if contains(curr_line,{'LHS','RTO','RHS','LTO'})
        gait_event_line = false;
    end
end
fclose(fid);

opts = detectImportOptions(fullfile(path2,filename_gait_events));
opts.VariableNamesLine = count;
opts.DataLines(1) = count+1;

gait_events = readtable(fullfile(path2,filename_gait_events),opts);

gait_events_aligned = gait_events;
for i = 1:length(gait_events.Properties.VariableNames)-2
    if ~isnan(aligned_data.alignment_times(5))
        gait_events_aligned.(gait_events.Properties.VariableNames{i}) = ...
            gait_events_aligned.(gait_events.Properties.VariableNames{i}) - aligned_data.alignment_times(5) + aligned_data.pre_align_time;
    else
        gait_events_aligned.(gait_events.Properties.VariableNames{i}) = ...
            gait_events_aligned.(gait_events.Properties.VariableNames{i}) - aligned_data.pre_align_time;
    end
end

aligned_data.gait_events = gait_events_aligned;

savename = strrep(filename_aligned_data,'.mat','_w_Gait_Events_Eleni_Adaptive.mat');
% savename = strrep(filename_aligned_data,'.mat','_w_Gait_Events_Hamid.mat');
% savename = strrep(filename_aligned_data,'.mat','_w_Gait_Events_Ken.mat');
save(fullfile(path,savename),'aligned_data');
end