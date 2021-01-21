function markLandmark(src,obj)
%
%
% Author:   Kenneth Louie
% Date:     12/29/20

if isnan(src.Parent.UserData.Marked_events(2)) || isnan(src.Parent.UserData.Marked_events(3))
    src.Parent.UserData.logger.String = [src.Parent.UserData.logger.String;...
        'Plot CoP data before marking landmarks.'];
    src.Parent.UserData.logger.Value = length(src.Parent.UserData.logger.String);
    return;
end

if ~isfield(src.Parent.UserData,'dataCursorObj')
    dcm = datacursormode(src.Parent);
    src.Parent.UserData.dataCursorObj = dcm;
else
    dcm = src.Parent.UserData.dataCursorObj;
end

% Remove all previous marks
dcm.removeAllDataCursors;

% Annotate current marked areas
startInd = src.Parent.UserData.CoP_startInd;
if ~isnan(src.Parent.UserData.Marked_events(4))
    L1_ind = src.Parent.UserData.Marked_events(4);
    if isnan(src.Parent.UserData.plot_markers(4))
        h = scatter(src.Parent.UserData.plot_axes(2),...
            src.Parent.UserData.Vicon_data.Cx_mm(startInd+L1_ind),...
            src.Parent.UserData.Vicon_data.Cy_mm(startInd+L1_ind),...
            60,'m','filled','DisplayName','L1');
        src.Parent.UserData.plot_markers(4) = h;
    else
        delete(src.Parent.UserData.plot_markers(4));
        h = scatter(src.Parent.UserData.plot_axes(2),...
            src.Parent.UserData.Vicon_data.Cx_mm(startInd+L1_ind),...
            src.Parent.UserData.Vicon_data.Cy_mm(startInd+L1_ind),...
            60,'m','filled','DisplayName','L1');
        src.Parent.UserData.plot_markers(4) = h;
    end
end

if ~isnan(src.Parent.UserData.Marked_events(5))
    L2_ind = src.Parent.UserData.Marked_events(5);
    if isnan(src.Parent.UserData.plot_markers(5))
        h = scatter(src.Parent.UserData.plot_axes(2),...
            src.Parent.UserData.Vicon_data.Cx_mm(startInd+L2_ind),...
            src.Parent.UserData.Vicon_data.Cy_mm(startInd+L2_ind),...
            60,'c','filled','DisplayName','L2');
        src.Parent.UserData.plot_markers(5) = h;
    else
        delete(src.Parent.UserData.plot_markers(5));
        h = scatter(src.Parent.UserData.plot_axes(2),...
            src.Parent.UserData.Vicon_data.Cx_mm(startInd+L2_ind),...
            src.Parent.UserData.Vicon_data.Cy_mm(startInd+L2_ind),...
            60,'c','filled','DisplayName','L2');
        src.Parent.UserData.plot_markers(5) = h;
    end
end

plotLegend = [];
if ~isnan(src.Parent.UserData.plot_markers(4))
    plotLegend = [plotLegend,src.Parent.UserData.plot_markers(4)];
end

if ~isnan(src.Parent.UserData.plot_markers(5))
    plotLegend = [plotLegend,src.Parent.UserData.plot_markers(5)];
end

if ~isempty(plotLegend)
    legend(src.Parent.UserData.plot_axes(2),plotLegend);
end

% Add callback
% (need to update ever button call because it will remember the first
% button pressed otherwise.)
dcm.UpdateFcn = {@setPoint,src};

% Switch ability to select data points
if strcmp(dcm.Enable,'on')
    dcm.Enable = 'off';
else
    dcm.Enable = 'on';
end
end