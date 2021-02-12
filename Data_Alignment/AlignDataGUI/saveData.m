function saveData(src,event)
aligned_data = src.Parent.Parent.UserData.aligned_data;
aligned_data.alignment_times = src.Parent.Parent.UserData.alignment_times;
aligned_data.alignment_subvariable_source = src.Parent.Parent.UserData.alignment_source;
aligned_data.alignment_source = {'Left IPG LFP','Right IPG LFP','Left IPG Accel','Right IPG Accel','Delsys','Xsens','Force Plate'};
aligned_data.files_used = src.Parent.Parent.UserData.file_names;
aligned_data.pre_align_time = src.Parent.Parent.UserData.pre_alignment_time;

[file, path] = uiputfile('*.mat');

save(fullfile(path,file),'aligned_data');

addEvent(sprintf('Aligned data saved to: %s',fullfile(path,file)));
end