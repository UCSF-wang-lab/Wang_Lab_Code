function plotCoP(src,obj)
% Author:   Kenneth Louie
% Date:     12/29/20

if isnan(src.Parent.UserData.Marked_events(2)) || isnan(src.Parent.UserData.Marked_events(3))
    src.Parent.UserData.logger.String = [src.Parent.UserData.logger.String;...
        'Not all makers set.'];
    src.Parent.UserData.logger.Value = length(src.Parent.UserData.logger.String);
    return;
end

vicon_time = (0:height(src.Parent.UserData.Vicon_data)-1)./1000;

% Find approximate index
[~,startInd] = min(abs(vicon_time-src.Parent.UserData.Marked_events(2)));
[~,endInd] = min(abs(vicon_time-src.Parent.UserData.Marked_events(3)));
src.Parent.UserData.CoP_startInd = startInd;
src.Parent.UserData.CoP_endInd = endInd;

% Plot
cla(src.Parent.UserData.plot_axes(2));
hold(src.Parent.UserData.plot_axes(2),'on');
plot(src.Parent.UserData.plot_axes(2),...
    src.Parent.UserData.Vicon_data.Cx_mm(startInd:endInd),...
    src.Parent.UserData.Vicon_data.Cy_mm(startInd:endInd),...
    '-k','linewidth',1.5);
scatter(src.Parent.UserData.plot_axes(2),...
    src.Parent.UserData.Vicon_data.Cx_mm(startInd),...
    src.Parent.UserData.Vicon_data.Cy_mm(startInd),...
    60,'ok','filled');

src.Parent.UserData.logger.String = [src.Parent.UserData.logger.String;...
    'Plotted CoP data'];
src.Parent.UserData.logger.Value = length(src.Parent.UserData.logger.String);

xlabel(src.Parent.UserData.plot_axes(2),'Mediolateral (mm)','FontWeight','bold');
ylabel(src.Parent.UserData.plot_axes(2),'Anteroposterior (mm)','FontWeight','bold');
end