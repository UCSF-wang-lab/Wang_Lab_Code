function loadData(src,~)
% General function to load in vicon or delsys data for APA marking.
% This function is only called as a callback to a load data button from the
% "markAPA" gui.
%
% Author:   Kenneth Louie
% Date:     02/10/2021

src.Parent.UserData.Marked_events = nan(1,5);
src.Parent.UserData.plot_markers = nan(1,5);

[aligned_data_filename,path] = uigetfile('.mat');
aligned_data = load(fullfile(path,aligned_data_filename));
src.Parent.UserData.Delsys_data = aligned_data.aligned_data.Delsys;
src.Parent.UserData.Vicon_data = aligned_data.aligned_data.FP;
src.Parent.UserData.Vicon_time = aligned_data.aligned_data.FP.Time;
src.Parent.UserData.logger.String = {sprintf('Loaded alignment file: %s',aligned_data_filename)};
src.Parent.UserData.logger.Value = length(src.Parent.UserData.logger.String);
src.Parent.UserData.save_path = path;
src.Parent.UserData.data_file = aligned_data_filename;
plotLoadedData(src.Parent);

end