function markMovement(src,obj)
%
%
% Author:   Kenneth Louie
% Date:     12/28/20

if ~isfield(src.Parent.UserData,'dataCursorObj')
    dcm = datacursormode(src.Parent);
    src.Parent.UserData.dataCursorObj = dcm;
else
    dcm = src.Parent.UserData.dataCursorObj;
end

% Remove all previous marks
dcm.removeAllDataCursors;

% Annotate current marked areas
if ~isnan(src.Parent.UserData.Marked_events(1))
    if isnan(src.Parent.UserData.plot_markers(1))
        h = xline(src.Parent.UserData.plot_axes(1),src.Parent.UserData.Marked_events(1),'--k','DisplayName','Cue');
        src.Parent.UserData.plot_markers(1) = h;
    else
        delete(src.Parent.UserData.plot_markers(1));
        h = xline(src.Parent.UserData.plot_axes(1),src.Parent.UserData.Marked_events(1),'--k','DisplayName','Cue');
        src.Parent.UserData.plot_markers(1) = h;
    end
end

if ~isnan(src.Parent.UserData.Marked_events(2))
    if isnan(src.Parent.UserData.plot_markers(2))
        h = xline(src.Parent.UserData.plot_axes(1),src.Parent.UserData.Marked_events(2),'--g','DisplayName','Movement onset');
        src.Parent.UserData.plot_markers(2) = h;
    else
        delete(src.Parent.UserData.plot_markers(2));
        h = xline(src.Parent.UserData.plot_axes(1),src.Parent.UserData.Marked_events(2),'--g','DisplayName','Movement onset');
        src.Parent.UserData.plot_markers(2) = h;
    end
end

if ~isnan(src.Parent.UserData.Marked_events(3))
    if isnan(src.Parent.UserData.plot_markers(3))
        h = xline(src.Parent.UserData.plot_axes(1),src.Parent.UserData.Marked_events(3),'--r','DisplayName','Stance leg toe-off');
        src.Parent.UserData.plot_markers(3) = h;
    else
        delete(src.Parent.UserData.plot_markers(3));
        h = xline(src.Parent.UserData.plot_axes(1),src.Parent.UserData.Marked_events(3),'--r','DisplayName','Stance leg toe-off');
        src.Parent.UserData.plot_markers(3) = h;
    end
end

plotLegend = [];
if ~isnan(src.Parent.UserData.plot_markers(1))
    plotLegend = [plotLegend,src.Parent.UserData.plot_markers(1)];
end

if ~isnan(src.Parent.UserData.plot_markers(2))
    plotLegend = [plotLegend,src.Parent.UserData.plot_markers(2)];
end

if ~isnan(src.Parent.UserData.plot_markers(3))
    plotLegend = [plotLegend,src.Parent.UserData.plot_markers(3)];
end

if ~isempty(plotLegend)
    legend(src.Parent.UserData.plot_axes(1),plotLegend);
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