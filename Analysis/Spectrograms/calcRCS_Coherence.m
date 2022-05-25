function calcRCS_Coherence(aligned_data,varargin)
%% RCSCoherogram
% Calculate the coherogram for the entire file. The coherogram will be
% between all possible recorded areas. This function can only take in one
% file at a time. 
%
% INPUTS:  Required
%               aligned_data        [=] Struct containing all the aligned
%                                       data from the trial of interests. 
%
%          Optional
%               saveDir             [=] Directory to save generated plots.
%                                       Default is an empty variable
%
%               savePlot            [=] Boolean option to save the 
%                                       resulting plot. Default is false.
%
%   Example call:
%           load(<filename>)
%           A = calcRCS_STFT(aligned_data,[],1,0.9,[]);
%           GaitEventSpectrogram(aligned_data,[],0);
%
% Date:     05/25/2022
% Author:   Kenneth H. Louie (kenneth.louie@ucsf.edu)
% Project:  MJFF aDBS Gait

%% Option variables
for i = 1:2:nargin-2
    switch varargin{i}
        case 'saveDir'
            saveDir = varargin{i+1};
        case 'savePlot'
            savePlot = varargin{i+1};
    end
end

% Set default options if not passed in by user
if ~exist('saveDir','var') || isempty(saveDir)
    saveDir = [];
end

if ~exist('savePlot','var') || isempty(savePlot)
    savePlot = 0;
end

%% Run 
elements_to_compare = [];
if isfield(aligned_data,'left_LFP_table')
    elements_to_compare = [elements_to_compare,1:4];
end

if isfield(aligned_data,'right_LFP_table')
    elements_to_compare = [elements_to_compare,5:8];
end

fig_vec = [];
pairs = nchoosek(elements_to_compare,2);
for i = 1:size(pairs,1)
    
    % Determine the indexing
    switch pairs(i,1)
        case {1,5}
            contact1 = 'key0';
        case {2,6}
            contact1 = 'key1';
        case {3,7}
            contact1 = 'key2';
        case {4,8}
            contact1 = 'key3';
    end
    
    switch pairs(i,2)
        case {1,5}
            contact2 = 'key0';
        case {2,6}
            contact2 = 'key1';
        case {3,7}
            contact2 = 'key2';
        case {4,8}
            contact2 = 'key3';
    end
    
    if pairs(i,1) <= 4
        side1 = 'left_LFP_table';
    else
        side1 = 'right_LFP_table';
    end
    
    if pairs(i,2) <= 4
        side2 = 'left_LFP_table';
    else
        side2 = 'right_LFP_table';
    end
    
    % Determine index limits
    n_val1 = height(aligned_data.(side1));
    n_val2 = height(aligned_data.(side2));
    [n_val,min_ind] = min([n_val1,n_val2]);
    
    if min_ind == 1
        if contains(side1,'left')
            time_vec = aligned_data.left_taxis;
        else
            time_vec = aligned_data.right_taxis;
        end
    else
        if contains(side2,'left')
            time_vec = aligned_data.left_taxis;
        else
            time_vec = aligned_data.right_taxis;
        end
    end
    
    % Get sample rate
    sr = aligned_data.DeviceSettings.Left.timeDomainSettings.samplingRate(end);
    
    % Calculate wavelet coherence and plot
    [x,~,y] = wcoherence(aligned_data.(side1).(contact1)(1:n_val),aligned_data.(side2).(contact2)(1:n_val),sr,'VoicesPerOctave',10);
    fig_vec(end+1) = figure();
    h = pcolor(time_vec,log2(y),x);
    h.EdgeColor = 'none';
    curr_y_ticks = yticks;
    yticklabels(cellstr(num2str(2.^curr_y_ticks')));
    title(createFigureTitle(side1,side2,contact1,contact2));
    ylabel('Frequency (Hz)');
    xlabel('Time (sec)');
end

% Format and save figures
if savePlot
    figure_format(6,6,12,[],'painters');
    
    % check if saving folders exist
    if ~isfolder(fullfile(saveDir,'Coherence'))
        mkdir(fullfile(saveDir,'Coherence'));
    end
    
    if ~isfolder(fullfile(saveDir,'Coherence',[aligned_data.stim_condition,'_STIM']))
        mkdir(fullfile(saveDir,'Coherence',[aligned_data.stim_condition,'_STIM']))
    end
    
    folders_to_check = {'FIG_files','PDF_files','TIFF_files'};
    extension = {'.fig','.pdf','.tiff'};
    for n = 1:length(folders_to_check)
        if ~isfolder(fullfile(saveDir,'Coherence',[aligned_data.stim_condition,'_STIM'],'CWT',folders_to_check{n}))
            mkdir(fullfile(saveDir,'Coherence',[aligned_data.stim_condition,'_STIM'],'CWT',folders_to_check{n}));
        end
    end
    
    for j = 1:length(fig_vec)
        curr_axes = gca(fig_vec(j));
        save_name = curr_axes.Parent.Children(1).Title.String;
        savefig(fig_vec(i),fullfile(saveDir,'Coherence',[aligned_data.stim_condition,'_STIM'],'CWT',folders_to_check{1},strrep(save_name,' ','_')));
        
        for m = 2:length(folders_to_check)
            print(fig_vec(j),[fullfile(saveDir,'Coherence',[aligned_data.stim_condition,'_STIM'],'CWT',folders_to_check{m},strrep(save_name,' ','_')),extension{m}],'-r300',['-d',extension{m}(2:end)]);
        end
    end
end
end

function out_str = createFigureTitle(side1,side2,contact1,contact2)
out_str = [];
if contains(side1,'left')
    out_str = [out_str,'Left '];
    switch contact1
        case 'key0'
            out_str = [out_str, '+2-0'];
        case 'key1'
            out_str = [out_str, '+3-1'];
        case 'key2'
            out_str = [out_str, '+9-8'];
        case 'key3'
            out_str = [out_str, '+11-10'];
    end
else
    out_str = [out_str,'Right '];
    switch contact1
        case 'key0'
            out_str = [out_str, '+2-0'];
        case 'key1'
            out_str = [out_str, '+3-1'];
        case 'key2'
            out_str = [out_str, '+9-8'];
        case 'key3'
            out_str = [out_str, '+11-10'];
    end
end

out_str = [out_str,' to '];

if contains(side2,'left')
    out_str = [out_str,'Left '];
    switch contact2
        case 'key0'
            out_str = [out_str, '+2-0'];
        case 'key1'
            out_str = [out_str, '+3-1'];
        case 'key2'
            out_str = [out_str, '+9-8'];
        case 'key3'
            out_str = [out_str, '+11-10'];
    end
else
    out_str = [out_str,'Right '];
    switch contact2
        case 'key0'
            out_str = [out_str, '+2-0'];
        case 'key1'
            out_str = [out_str, '+3-1'];
        case 'key2'
            out_str = [out_str, '+9-8'];
        case 'key3'
            out_str = [out_str, '+11-10'];
    end
end

end