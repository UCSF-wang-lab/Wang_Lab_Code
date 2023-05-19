function changeRCSTimeUnit(src,event)

for i = 1:length(src.Parent.Parent.UserData.plot_axes)
    if contains(src.Parent.Parent.UserData.plot_axes(i).Title.String,{'Left','Right'})
        plotData(src.Parent.Parent.UserData.plot_options.sources(i,2));
    end
end


end