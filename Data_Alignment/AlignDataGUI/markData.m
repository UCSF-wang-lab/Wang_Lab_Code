function markData(src,obj)
%
%
% Author:   Kenneth Louie
% Date:     12/28/20

if ~isfield(src.Parent.Parent.UserData,'dataCursorObj')
    dcm = datacursormode(src.Parent.Parent);
    src.Parent.Parent.UserData.dataCursorObj = dcm;
else
    dcm = src.Parent.Parent.UserData.dataCursorObj;
end

% Remove all previous marks
dcm.removeAllDataCursors;

% Annotate current marked times
plotAlignmentTimes();

% Add callback
% (need to update ever button call because it will remember the first
% button pressed otherwise.)
dcm.UpdateFcn = @setPoint;

% Switch ability to select data points
if strcmp(dcm.Enable,'on')
    dcm.Enable = 'off';
else
    dcm.Enable = 'on';
end
end