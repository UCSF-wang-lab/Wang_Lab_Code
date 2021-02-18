function combineAlignedDataGaitTimings()

[aligned_data_filename,aligned_data_path] = uigetfile('*.mat');
load(fullfile(aligned_data_path,aligned_data_filename));

[gait_marker_filename,gait_marker_path] = uigetfile('*.txt');
gait_events = readtable(fullfile(gait_marker_path,gait_marker_filename));

temp = table2array(gait_events);
adjusted_gait_markers = temp - aligned_data.alignment_times(5) + aligned_data.pre_align_time;
adjusted_gait_markers = array2table(adjusted_gait_markers,'VariableNames',gait_events.Properties.VariableNames);

aligned_data.gait_events = adjusted_gait_markers;

save(fullfile(aligned_data_path,aligned_data_filename),'aligned_data');
end