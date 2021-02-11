function setZoom(src,~)
% Callback function when pressing a toggle button from the "AlignData" gui.
% When toggled, will enable or disable zoom for all plots.
%
% Author:   Kenneth Louie
% Date:     01-22-2021

if ~isfield(src.Parent.Parent.UserData,'zoomObj')
    z = zoom(src.Parent.Parent);
    src.Parent.Parent.UserData.zoomObj = z;
else
    z = src.Parent.Parent.UserData.zoomObj;
end

if strcmp(src.Parent.Parent.UserData.zoomObj.Enable,'on')
    src.Parent.Parent.UserData.zoomObj.Enable = 'off';
    for i = 1:length(src.Parent.Parent.UserData.plot_axes)
        setAllowAxesZoom(z,src.Parent.Parent.UserData.plot_axes(i),false);
    end
else
    src.Parent.Parent.UserData.zoomObj.Enable = 'on';
    for i = 1:length(src.Parent.Parent.UserData.plot_axes)
        setAllowAxesZoom(z,src.Parent.Parent.UserData.plot_axes(i),true);
    end
end

end