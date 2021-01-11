function loadData(src,~)
% General function to load in vicon or delsys data for APA marking.
% This function is only called as a callback to a load data button from the
% "markAPA" gui.
%
% Author:   Kenneth Louie
% Date:     12/28/20

src.Parent.UserData.Marked_events = nan(1,5);
src.Parent.UserData.plot_markers = nan(1,5);

if contains(src.String,'Vicon')
    [vicon_filename,path] = uigetfile('*ForcePlateData.csv');
    src.Parent.UserData.save_path = path;
    [~,src.Parent.UserData.Vicon_filename,~] = fileparts(vicon_filename);
    
    src.Parent.UserData.Vicon_data = readtable(fullfile(path,vicon_filename));
    if isempty(src.Parent.UserData.logger.String)
        src.Parent.UserData.logger.String = {sprintf('Loaded Vicon file: %s',vicon_filename)};
        src.Parent.UserData.logger.Value = length(src.Parent.UserData.logger.String);
    else
        src.Parent.UserData.logger.String = [src.Parent.UserData.logger.String;sprintf('Loaded Vicon file: %s',vicon_filename)];
        src.Parent.UserData.logger.Value = length(src.Parent.UserData.logger.String);
    end
    src.Parent.UserData.Vicon_time = [];
    plotLoadedData(src.Parent);
else
    [delsys_filename,path] = uigetfile('.mat');
    src.Parent.UserData.Delsys_data = load(fullfile(path,delsys_filename));
    [~,src.Parent.UserData.Delsys_filename,~] = fileparts(delsys_filename);
    
    if isempty(src.Parent.UserData.logger.String)
        src.Parent.UserData.logger.String = {sprintf('Loaded Delsys file: %s',delsys_filename)};
        src.Parent.UserData.logger.Value = length(src.Parent.UserData.logger.String);
    else
        src.Parent.UserData.logger.String = [src.Parent.UserData.logger.String;sprintf('Loaded Delsys file: %s',delsys_filename)];
        src.Parent.UserData.logger.Value = length(src.Parent.UserData.logger.String);
    end
    plotLoadedData(src.Parent);
end

end