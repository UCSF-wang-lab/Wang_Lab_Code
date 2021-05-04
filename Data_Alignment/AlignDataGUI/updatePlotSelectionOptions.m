function updatePlotSelectionOptions(src,event)
main = findobj('Tag','AlignData_main_window');

% If src and event is empty, it was called from load data event
if isempty(src) &&  isempty(event)
    str = getLoadedData(main.UserData);
    for i = 1:size(main.UserData.plot_options.sources,1)
        if length(main.UserData.plot_options.sources(i,1).String) == 1
            main.UserData.plot_options.sources(i,1).Enable = 'on';
            current_data_source = [];
        else
            current_data_source = main.UserData.plot_options.sources(i,1).String{main.UserData.plot_options.sources(i,1).Value};
        end
        main.UserData.plot_options.sources(i,1).String = [' ';str];
        
        if ~isempty(current_data_source)
            adjusted_value = find(cellfun(@(x) strcmp(current_data_source,x),main.UserData.plot_options.sources(i,1).String));
            main.UserData.plot_options.sources(i,1).Value = adjusted_value;
        end
    end
    return; % Exit function early
end

% Function was called after changing selection on plot options
selection = src.String{src.Value};
if contains(selection,'LFP')
    if contains(selection,'Left')
        str = getSubVariables(main,'LFP','Left');
    elseif contains(selection,'Right')
        str = getSubVariables(main,'LFP','Right');
    end
elseif contains(selection,'Accel')
    if contains(selection,'Left')
        str = getSubVariables(main,'Accel','Left');
    elseif contains(selection,'Right')
        str = getSubVariables(main,'Accel','Right');
    end
elseif strcmp(selection,'Delsys')
    str = getSubVariables(main,'Delsys');
elseif strcmp(selection,'Xsens')
    str = getSubVariables(main,'Xsens');
elseif strcmp(selection,'Force plate')
    str = getSubVariables(main,'FP');
elseif strcmp(selection,'Teensey')
    str = getSubVariables(main,'Teensey');
else
    str = [];
end

if contains(src.Parent.Title,'1')
    main.UserData.plot_options.sources(1,2).String = [' ';str];
    main.UserData.plot_options.sources(1,2).Enable = 'on';
    main.UserData.plot_options.sources(1,2).Value = 1;
elseif contains(src.Parent.Title,'2')
    main.UserData.plot_options.sources(2,2).String = [' ';str];
    main.UserData.plot_options.sources(2,2).Enable = 'on';
    main.UserData.plot_options.sources(2,2).Value = 1;
elseif contains(src.Parent.Title,'3')
    main.UserData.plot_options.sources(3,2).String = [' ';str];
    main.UserData.plot_options.sources(3,2).Enable = 'on';
    main.UserData.plot_options.sources(3,2).Value = 1;
elseif contains(src.Parent.Title,'4')
    main.UserData.plot_options.sources(4,2).String = [' ';str];
    main.UserData.plot_options.sources(4,2).Enable = 'on';
    main.UserData.plot_options.sources(4,2).Value = 1;
end

end

function names = getLoadedData(data_struct)
names = {};

% LFP data
if ~isempty(data_struct.LFP_data)
    if isfield(data_struct.LFP_data,'Left')
        names = [names;'Left LFP'];
    end
    if isfield(data_struct.LFP_data,'Right')
        names = [names;'Right LFP'];
    end
end

% Acceleration data
if ~isempty(data_struct.Accel_data)
    if isfield(data_struct.Accel_data,'Left')
        names = [names;'Left Accel'];
    end
    
    if isfield(data_struct.Accel_data,'Right')
        names = [names;'Right Accel'];
    end
end

% Delsys data
if ~isempty(data_struct.Delsys_data)
    names = [names;'Delsys'];
end

% Xsens data
if ~isempty(data_struct.Xsens_data)
    names = [names;'Xsens'];
end

% Force plate data
if ~isempty(data_struct.FP_data)
    names = [names;'Force plate'];
end

% Teensey data
if ~isempty(data_struct.Teensey_data)
    names = [names;'Teensey'];
end

end

function sub_variables = getSubVariables(varargin)
sub_variables = [];

% LFP and Accel have a third input
switch varargin{2}
    case 'LFP'
        all_names = varargin{1}.UserData.([varargin{2},'_data']).(varargin{3}).timeDomainDataTable.Properties.VariableNames;
        matches = cellfun(@(x) contains(x,'key'),all_names);
        sub_variables = all_names(matches);
    case 'Accel'
        all_names = varargin{1}.UserData.([varargin{2},'_data']).(varargin{3}).accelDataTable.Properties.VariableNames;
        matches = cellfun(@(x) contains(x,'Samples'),all_names);
        sub_variables = all_names(matches);
    case 'Delsys'
        sub_variables = varargin{1}.UserData.Delsys_data.out_struct.Chan_names;
    case 'Xsens'
        sub_variables = varargin{1}.UserData.Xsens_data.Properties.VariableNames;
    case 'FP'
        all_names = varargin{1}.UserData.FP_data.Properties.VariableNames;
        matches = cellfun(@(x) contains(x,'_'),all_names);
        sub_variables = all_names(matches);
    case 'Teensey'
        sub_variables = varargin{1}.UserData.Teensey_data.Properties.VariableNames;
end

sub_variables = sub_variables(:);
end