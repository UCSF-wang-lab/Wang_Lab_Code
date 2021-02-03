function createPlotPanel(main_window)
%% Plots
plot_panel = uipanel('BackgroundColor','white',...
    'BorderType','none',...
    'Units','normalized',...
    'Position',[0.215 0.2 0.67 0.79]);

h(1) = axes('Parent',plot_panel,'Units','normalized','Position',[0.05,0.545,0.44,0.42]);
h(2) = axes('Parent',plot_panel,'Units','normalized','Position',[0.54,0.545,0.44,0.42]);
h(3) = axes('Parent',plot_panel,'Units','normalized','Position',[0.05,0.05,0.44,0.42]);
h(4) = axes('Parent',plot_panel,'Units','normalized','Position',[0.54,0.05,0.44,0.42]);

% Annotation


% Debugging
for i = 1:length(h)
    % Get plot axes limits
    x = xlim(h(i));
    y = ylim(h(i));
    
    % Show plot number
    text(h(i),diff(x)/2,diff(y)/2,num2str(i),'FontSize',80,'HorizontalAlignment','center');
    
    % Debugging
%     ylabel(h(i),'Voltage'); 
%     xlabel(h(i),'Time');
%     yticks(h(i),[0.001:0.005]);
%     ylim(h(i),[0.001,0.05]);
%     title(h(i),'test');
end

% Set user data to plot data to these axes
main_window.UserData.plot_axes = h;
end