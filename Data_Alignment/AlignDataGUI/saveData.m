function saveData(src,event)
aligned_data = src.Parent.Parent.UserData.aligned_data;
aligned_data.alignment_times = src.Parent.Parent.UserData.alignment_times;
aligned_data.alignment_subvariable_source = src.Parent.Parent.UserData.alignment_source;
aligned_data.alignment_source = {'Left IPG LFP','Right IPG LFP','Left IPG Accel','Right IPG Accel','Delsys','Xsens','Force Plate','Teensey'};
aligned_data.files_used = src.Parent.Parent.UserData.file_names;
aligned_data.pre_align_time = src.Parent.Parent.UserData.pre_alignment_time;
aligned_data.DeviceSettings = src.Parent.Parent.UserData.DeviceSettings;
aligned_data.stim_condition = src.Parent.Parent.UserData.stim_condition;
aligned_data.med_condition = src.Parent.Parent.UserData.med_condition;
aligned_data.trial_num = src.Parent.Parent.UserData.trial_num;

[file, path] = uiputfile('*.mat');
save_filename = fullfile(path,file);

save(save_filename,'aligned_data');

addEvent(sprintf('Aligned data saved to: %s',save_filename));
end