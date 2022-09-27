function combineAlignedDataGaitEvents()
% load Aligned data mat file
[filename_aligned_data,path] = uigetfile('*.mat');
load(fullfile(path,filename_aligned_data));

% load Gait event file
[filename_gait_events,path2] = uigetfile('*.csv');
if filename_gait_events == 0
    [filename_gait_events,path2] = uigetfile('*.txt');
end

fid = fopen(fullfile(path2,filename_gait_events));
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
for i = 1:length(gait_events.Properties.VariableNames)
    gait_events_aligned.(gait_events.Properties.VariableNames{i}) = ...
        gait_events_aligned.(gait_events.Properties.VariableNames{i}) - aligned_data.alignment_times(5) + aligned_data.pre_align_time;
end

aligned_data.gait_events = gait_events_aligned;

% savename = strrep(filename_aligned_data,'.mat','_w_Gait_Events.mat');
savename = strrep(filename_aligned_data,'.mat','_w_Gait_Events_Julia.mat');
save(fullfile(path,savename),'aligned_data');
end