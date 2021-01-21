function setZoom(src,~)
% Callback function when pressing a toggle button from the "markAPA" gui.
% When toggled, will enable or disable zoom for the first plot axes only.
% The first plot on the GUI plots the raw force plate data from Vicon and
% TTL pulse collected by Delsys.
%
% Author:   Kenneth Louie
% Date:     12/28/20

if ~isfield(src.Parent.UserData,'zoomObj')
    z = zoom(src.Parent);
    src.Parent.UserData.zoomObj = z;
else
    z = src.Parent.UserData.zoomObj;
end

if strcmp(src.Parent.UserData.zoomObj.Enable,'on')
    src.Parent.UserData.zoomObj.Enable = 'off';
else
    src.Parent.UserData.zoomObj.Enable = 'on';
    setAllowAxesZoom(z,src.Parent.UserData.plot_axes(2),false);
    setAllowAxesZoom(z,src.Parent.UserData.plot_axes(3),false);
    
end

end