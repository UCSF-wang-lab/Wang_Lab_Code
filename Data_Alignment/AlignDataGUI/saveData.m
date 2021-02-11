function saveData(src,event)
aligned_data = src.Parent.Parent.UserData.aligned_data;
aligned_data.alignment_times = src.Parent.Parent.UserData.alignment_times;
aligned_data.alignment_subvariable_source = src.Parent.Parent.UserData.alignment_source;
aligned_data.alignment_source = {'Left IPG LFP','Right IPG LFP','Left IPG Accel','Right IPG Accel','Delsys','Xsens','Force Plate'};

[file, path] = uiputfile('*.mat');

save(fullfile(path,file),'aligned_data');
end