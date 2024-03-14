function diffCoherenceGSLT(varargin)
%% diffCoherenceGSLT
% This function plots the difference coherogram between 2 conditions.
%
% INPUTS:  Optional
%
%               dataPath[=]     The path to the parent folder containing the 
%                               grand coherence data of both condition 1
%                               and 2. If not provided as an argument, the
%                               user will need to provide it through the UI
%                               prompt.
%
%                  cond1[=]     String describing the condition 1, 
%                               e.g. OFFmed_ONdbs_LeftMaladaptive. 
%
%                  cond2[=]     String describing the condition 2, 
%                               e.g. OFFmed_ONdbs_LeftAdaptive. 
%
%               cohPairs[=]     Pairs to compute coherence between. The
%                               list of pairs will be used across all
%                               files. Default is to do all possible
%                               unilateral pairs.
%               
%                   keys[=]     Cell array of keys to analyze. Default is
%                               to analyze all keys.
%
%               savePlot[=]     Boolean option to save the resulting merged plots.
%                               Default is false.       
%
%
%   Example call:
%           keys = {'key0','key1','key2','key3'};
%           diffCoherenceGSLT('cond1','OFFmed_ONdbs_LeftMaladaptive','cond2','OFFmed_ONdbs_LeftAdaptive','condDiff','Left Ada-Mal','keys',keys,'savePlot',1)
%
% Date:     10/29/2023
% Author: Eleni Patelaki (eleni.patelaki@ucsf.edu)

for i = 1:2:nargin   
    switch varargin{i}
        case 'dataPath'
            dataPath = varargin{i+1};
        case 'cond1'
            cond1 = varargin{i+1};
        case 'cond2'
            cond2 = varargin{i+1};
        case 'condDiff'
            condDiff = varargin{i+1};
        case 'cohPairs'
            cohPairs = varargin{i+1};
        case 'keys'
            keys = varargin{i+1};
        case 'savePlot'
            savePlot = varargin{i+1};
    end
end
% Set default options if not passed in by user
if ~exist('keys','var') || isempty(keys)
    keys = {'key0','key1','key2','key3'};
end

if ~exist('cohPairs','var') || isempty(cohPairs)
    cohPairs = nchoosek(keys,2);
end

if ~exist('cond1','var') || isempty(cond1)
    cond1 = 'OFFmed_ONdbs_LeftMaldaptive';
end

if ~exist('cond2','var') || isempty(cond2)
    cond2 = 'OFFmed_ONdbs_LeftAptive';
end

if ~exist('condDiff','var') || isempty(condDiff)
    condDiff = 'Left Ada-Mal';
end

if ~exist('savePlot','var') || isempty(savePlot)
    savePlot = 0;   % Does not save plot by default.
else
    correct_path = false;
    while ~correct_path
        saveplot_Path = uigetdir();
        if ismac
            split_parent_dir = strsplit(saveplot_Path,'/');
        elseif ispc
            split_parent_dir = strsplit(saveplot_Path,'\');
        else
            error('Platform not supported');
        end
            
        if strcmp(split_parent_dir{end},'Data')
            correct_path = true;
        else
            warning('Please select the folder called "Data"');
        end
    end
    
    plot_save_path = fullfile(saveplot_Path,'Figures','GSLT_coherence_figs_adaptive','Merged',condDiff);
    if ~exist(plot_save_path, 'dir')
       mkdir(plot_save_path);
    end
end

% Fetch & load the grand average data for conditions 1 & 2
if ~exist('dataPath','var') || isempty(dataPath)
    correct_path = false;
    while ~correct_path
        dataPath = uigetdir();
        if ismac
            split_parent_dir = strsplit(dataPath,'/');
        elseif ispc
            split_parent_dir = strsplit(dataPath,'\');
        else
            error('Platform not supported');
        end
            
        if strcmp(split_parent_dir{end},'GSLT_coherence_data_adaptive')
            correct_path = true;
        else
            warning('Please select the folder called "GSLT_coherence_data_adaptive"');
        end
    end
end

GA1 = load(fullfile(dataPath,cond1,'grandAvgCoh.mat'));
GA1 = GA1.grandAvgCoh;
GA2 = load(fullfile(dataPath,cond2,'grandAvgCoh.mat'));
GA2 = GA2.grandAvgCoh;

% Load the frequency vector
fc = load(fullfile(dataPath,cond1,'fc_coh_part1.mat'));
fc = fc.fc;

%% Plotting
%Left
for i = 1:length(GA1.Left)
    freq_vec = squeeze(fc.Left{1}(:,:,1));

    figure;

    ax = pcolor(1:100,log2(freq_vec),GA2.Left{i}-GA1.Left{i});
    ticks = logspace(log10(2.5),log10(50),5);
    ax.Parent.YTick = log2(ticks);
    ax.Parent.YTickLabel = sprintf('%1.2f\n',round(ticks,2));
    ylim([log2(2.5),log2(50)]);
    xticks([1,10,20,30,40,50,60,70,80,90,100]);
    xticklabels({'0','10','20','30','40','50','60','70','80','90','100'});
    shading interp;
    colormap jet;
    ylabel('Frequency (Hz)');
    title({'Left';sprintf('%s to %s Grand Average Difference Coherence',cohPairs{i,1},cohPairs{i,2})});

    % Save plots
    if savePlot    
        saveas(gcf,fullfile(plot_save_path,strcat('Left_',cohPairs{i,1},'_',cohPairs{i,2})),'fig');
        saveas(gcf,fullfile(plot_save_path,strcat('Left_',cohPairs{i,1},'_',cohPairs{i,2})),'tiff');
    end
end

%Right
for i = 1:length(GA1.Right)
    freq_vec = squeeze(fc.Right{1}(:,:,1));

    figure;

    ax = pcolor(1:100,log2(freq_vec),GA2.Right{i}-GA1.Right{i});
    ticks = logspace(log10(2.5),log10(50),5);
    ax.Parent.YTick = log2(ticks);
    ax.Parent.YTickLabel = sprintf('%1.2f\n',round(ticks,2));
    ylim([log2(2.5),log2(50)]);
    xticks([1,10,20,30,40,50,60,70,80,90,100]);
    xticklabels({'0','10','20','30','40','50','60','70','80','90','100'});
    shading interp;
    colormap jet;
    ylabel('Frequency (Hz)');
    title({'Right';sprintf('%s to %s Grand Average Difference Coherence',cohPairs{i,1},cohPairs{i,2})});

    % Save plots
    if savePlot    
        saveas(gcf,fullfile(plot_save_path,strcat('Right_',cohPairs{i,1},'_',cohPairs{i,2})),'fig');
        saveas(gcf,fullfile(plot_save_path,strcat('Right_',cohPairs{i,1},'_',cohPairs{i,2})),'tiff');
    end
end

end