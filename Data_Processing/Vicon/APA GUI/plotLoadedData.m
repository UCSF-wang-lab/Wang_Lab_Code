function plotLoadedData(main_window_handle)
% Function that is called once Vicon or Delsys data is loaded. Takes
% advantage of dual y axes. The left y axis is for raw Vicon force plate
% data. The right y axis is for the raw TTL pulse collected from Delsys.
% Note, the ground reaction force (Fz) is baseline subtracted so the left y
% axis range does not blow up.
%
% Author:   Kenneth Louie
% Date:     12/28/20

% Clear any previous plots
cla(main_window_handle.UserData.plot_axes(1));
set(main_window_handle.UserData.plot_axes(1),'ytick',[],'yticklabels',[],'xcolor','k','ycolor','k','linewidth',1.5);
cla(main_window_handle.UserData.plot_axes(2));
set(main_window_handle.UserData.plot_axes(2),'xcolor','k','ycolor','k','linewidth',1.5);

colors = colororder;

% Plot Vicon
if ~isempty(main_window_handle.UserData.Vicon_data)
    if isempty(main_window_handle.UserData.Vicon_time)
        vicon_time = (0:height(main_window_handle.UserData.Vicon_data)-1)./1000;
        main_window_handle.UserData.Vicon_time = vicon_time;
    else
        vicon_time = main_window_handle.UserData.Vicon_time;
    end

    baseline_Fz_N = mean(main_window_handle.UserData.Vicon_data.Fz_N(1:10000));
    
    myYLim = [-50,...
        max([max(main_window_handle.UserData.Vicon_data.Fx_N),...
        max(main_window_handle.UserData.Vicon_data.Fy_N),...
        max(main_window_handle.UserData.Vicon_data.Fz_N-baseline_Fz_N)])*1.05];
    
    yyaxis(main_window_handle.UserData.plot_axes(1),'left');
    cla(main_window_handle.UserData.plot_axes(1));
    hold(main_window_handle.UserData.plot_axes(1),'on');
    
    plot(main_window_handle.UserData.plot_axes(1),...
        vicon_time,...
        main_window_handle.UserData.Vicon_data.Fx_N,'Color',colors(1,:),...
        'LineStyle','-');
    plot(main_window_handle.UserData.plot_axes(1),...
        vicon_time,...
        main_window_handle.UserData.Vicon_data.Fy_N,'Color',colors(2,:),...
        'LineStyle','-');
    plot(main_window_handle.UserData.plot_axes(1),...
        vicon_time,...
        main_window_handle.UserData.Vicon_data.Fz_N-baseline_Fz_N,'Color',colors(3,:),...
        'LineStyle','-');
    ylim(main_window_handle.UserData.plot_axes(1),myYLim);
end

% Plot Delsys
if ~isempty(main_window_handle.UserData.Delsys_data)
    TTL_ind = find(contains(main_window_handle.UserData.Delsys_data.out_struct.Chan_names,'TTL'));
    chan_name = main_window_handle.UserData.Delsys_data.out_struct.Chan_names{TTL_ind};
    
    yyaxis(main_window_handle.UserData.plot_axes(1),'right');
    cla(main_window_handle.UserData.plot_axes(1));
    
    plot(main_window_handle.UserData.plot_axes(1),...
        main_window_handle.UserData.Delsys_data.out_struct.Time.(chan_name),...
        main_window_handle.UserData.Delsys_data.out_struct.Data.(chan_name),...
        'Color',colors(4,:));
end
xlabel(main_window_handle.UserData.plot_axes(1),'Time (s)','FontWeight','bold');

% Remove y ticks
yyaxis(main_window_handle.UserData.plot_axes(1),'left');
set(main_window_handle.UserData.plot_axes(1),'YTick',[],'YTickLabels',[],'ycolor','k');
yyaxis(main_window_handle.UserData.plot_axes(1),'right');
set(main_window_handle.UserData.plot_axes(1),'YTick',[],'YTickLabels',[],'ycolor','k');


end