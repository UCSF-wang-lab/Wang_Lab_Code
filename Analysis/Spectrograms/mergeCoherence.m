function mergeCoherence(varargin)
%% MergeCoherence
% When multiple aligned datasets exist within the same experimental session, 
% this function merges the multiple occuring coherence decompositions.
%
% INPUTS:  Optional
%               dataPath    [=] The path to the folder containing the 
%                               coherence decompositions.
%
%               cohPairs    [=] Pairs to compute coherence between. The
%                               list of pairs will be used across all
%                               files. Default is to do all possible
%                               unilateral pairs.
%               
%               keys        [=] Cell array of keys to analyze. Default is
%                               to analyze all keys.
%
%               savePlot    [=] Boolean option to save the resulting merged plots.
%                               Default is false.       
%
%               condition   [=] String describing the condition analyzed, 
%                               e.g. OFFmed_ONdbs. It creates appropriate subfolders
%                               to save the plots in
%
%   Example call:
%           keys = {'key0','key1','key2','key3'};
%           mergeCoherence('keys',keys,'savePlot',1,'condition','OFFmed_ONdbs')
%
% Date:     10/24/2023
% Author: Eleni Patelaki (eleni.patelaki@ucsf.edu)
% Project:  MJFF aDBS Gait

%% Option variables
for i = 1:2:nargin
    switch varargin{i}
        case 'dataPath'
            dataPath = varargin{i+1};
        case 'cohPairs'
            cohPairs = varargin{i+1};
        case 'keys'
            keys = varargin{i+1};
        case 'savePlot'
            savePlot = varargin{i+1};
        case 'condition'
            condition = varargin{i+1};
    end
end

% Set default options if not passed in by user
if ~exist('keys','var') || isempty(keys)
    keys = {'key0','key1','key2','key3'};
end

if ~exist('cohPairs','var') || isempty(cohPairs)
    cohPairs = nchoosek(keys,2);
end

if ~exist('condition','var')
    condition = '';
end

if ~exist('savePlot','var') || isempty(savePlot)
    savePlot = 0;   % Does not save plot by default.
else
    correct_path = false;
    while ~correct_path
        parent_dir = uigetdir();
        if ismac
            split_parent_dir = strsplit(parent_dir,'/');
        elseif ispc
            split_parent_dir = strsplit(parent_dir,'\');
        else
            error('Platform not supported');
        end
            
        if strcmp(split_parent_dir{end},'Data')
            correct_path = true;
        else
            warning('Please select the folder called "Data"');
        end
    end
    
    plot_save_path = fullfile(parent_dir,'Figures','GSLT_coherence_figs','Merged',condition);
    if ~exist(plot_save_path, 'dir')
       mkdir(plot_save_path);
    end
end

% Fetch all coherence datafiles
if ~exist('dataPath','var') || isempty(dataPath)
    [cohdatafiles,dataPath] = uigetfile('gc*.mat','Select file(s)', 'MultiSelect', 'on');
else
    cohdatafiles = dir(fullfile(dataPath, 'gc*.mat'));
end

% Fetch all frequency & normalization datafiles
freqfiles = dir(fullfile(dataPath, 'fc*.mat'));
normfiles = dir(fullfile(dataPath, 'normalization*.mat'));

% check that the number of the frequency & normalization datafiles
% is equal to that of the coherence datafiles
assert((numel(cohdatafiles)==numel(freqfiles))& (numel(freqfiles)==numel(normfiles)));

% Initialize the aggregated coherence, normalization structs and gait
% cycle number structures
allgc.Left = cell(1,size(cohPairs,1));
allgc.Right = cell(1,size(cohPairs,1));
allfc.Left = cell(1,size(cohPairs,1));
allfc.Right = cell(1,size(cohPairs,1));
allnorm.Left = cell(1,size(cohPairs,1));
allnorm.Right = cell(1,size(cohPairs,1));
allgcNums.Left = cell(1,size(cohPairs,1));
allgcNums.Right = cell(1,size(cohPairs,1));

% Merge the coherence & normalization datafiles
for i = 1:numel(cohdatafiles)
    if isstruct(cohdatafiles)
        load(fullfile(dataPath,getfield(cohdatafiles(i),'name')));
    else
        load(fullfile(dataPath,cohdatafiles{i}));
    end
    load(fullfile(dataPath,getfield(freqfiles(i),'name')));
    load(fullfile(dataPath,getfield(normfiles(i),'name')));

    % Ger the number of gait cycles in each part
    gcNums.Left = cellfun(@(x) size(x,3),gc.Left,'UniformOutput',false);
    gcNums.Right = cellfun(@(x) size(x,3),gc.Right,'UniformOutput',false);

    allgc.Left = cellfun(@(x,y) cat(3,x,y),allgc.Left,gc.Left,'UniformOutput',false);
    allgc.Right = cellfun(@(x,y) cat(3,x,y),allgc.Right,gc.Right,'UniformOutput',false);
    allfc.Left = cellfun(@(x,y) cat(3,x,y),allfc.Left,fc.Left,'UniformOutput',false);
    allfc.Right = cellfun(@(x,y) cat(3,x,y),allfc.Right,fc.Right,'UniformOutput',false);
    allnorm.Left = cellfun(@(x,y) cat(3,x,y),allnorm.Left,normalization.Left,'UniformOutput',false);
    allnorm.Right = cellfun(@(x,y) cat(3,x,y),allnorm.Right,normalization.Right,'UniformOutput',false);
    allgcNums.Left = cellfun(@(x,y) cat(2,x,y),allgcNums.Left,gcNums.Left,'UniformOutput',false);
    allgcNums.Right = cellfun(@(x,y) cat(2,x,y),allgcNums.Right,gcNums.Right,'UniformOutput',false);

    clear gc gcNums normalization 
end


% Calculate the weighted average of the normalizations matrices &
% Normalize
for i = 1:size(cohPairs)
    for j = 1:numel(cohdatafiles)
        if j==1
            sumNormLeft =  allgcNums.Left{i}(j)*allnorm.Left{i}(:,:,j);
            sumNormRight =  allgcNums.Right{i}(j)*allnorm.Right{i}(:,:,j);
        else
            sumNormLeft =  sumNormLeft + allgcNums.Left{i}(j)*allnorm.Left{i}(:,:,j);
            sumNormRight =  sumNormRight + allgcNums.Right{i}(j)*allnorm.Right{i}(:,:,j);
        end
    end
    avgnormLeft = sumNormLeft/sum(allgcNums.Left{i},'all');
    avgnormRight = sumNormRight/sum(allgcNums.Right{i},'all');

    % Apply the normalization (z-scoring)
    normCoh.Left{i} = (allgc.Left{i}-avgnormLeft(:,1))./avgnormLeft(:,2);
    normCoh.Right{i} = (allgc.Right{i}-avgnormRight(:,1))./avgnormRight(:,2);

    clear avgnormLeft avgnormRight sumNormLeft sumNormRight
end

% Calculate grand averages
for i = 1:length(normCoh.Left)
    grandAvgCoh.Left{i} = mean(normCoh.Left{i},3);
end

for i = 1:length(normCoh.Right)
    grandAvgCoh.Right{i} = mean(normCoh.Right{i},3);
end

%% Plotting
for i = 1:length(grandAvgCoh.Left)
    freq_vec = squeeze(allfc.Left{1}(:,:,1));

    figure;

    ax = pcolor(1:100,log2(freq_vec),grandAvgCoh.Left{i});
    ticks = logspace(log10(2.5),log10(50),5);
    ax.Parent.YTick = log2(ticks);
    ax.Parent.YTickLabel = sprintf('%1.2f\n',round(ticks,2));
    ylim([log2(2.5),log2(50)]);
    xticks([1,10,20,30,40,50,60,70,80,90,100]);
    xticklabels({'0','10','20','30','40','50','60','70','80','90','100'});
    shading interp;
    colormap jet;
    ylabel('Frequency (Hz)');
    title({'Left';sprintf('%s to %s Grand Average Coherence',cohPairs{i,1},cohPairs{i,2})});

    % Save plots
    if savePlot    
        saveas(gcf,fullfile(plot_save_path,strcat('Left_',cohPairs{i,1},'_',cohPairs{i,2})),'fig');
        saveas(gcf,fullfile(plot_save_path,strcat('Left_',cohPairs{i,1},'_',cohPairs{i,2})),'tiff');
    end
end

for i = 1:length(grandAvgCoh.Right)
    freq_vec = squeeze(allfc.Left{1}(:,:,1));

    figure;
    ax = pcolor(1:100,log2(freq_vec),grandAvgCoh.Right{i});
    ticks = logspace(log10(2.5),log10(50),5);
    ax.Parent.YTick = log2(ticks);
    ax.Parent.YTickLabel = sprintf('%1.2f\n',round(ticks,2));
    ylim([log2(2.5),log2(50)]);
    xticks([1,10,20,30,40,50,60,70,80,90,100]);
    xticklabels({'0','10','20','30','40','50','60','70','80','90','100'});
    shading interp;
    colormap jet;
    ylabel('Frequency (Hz)');
    title({'Right';sprintf('%s to %s Grand Average Coherence',cohPairs{i,1},cohPairs{i,2})});

    % Save plots
    if savePlot
        saveas(gcf,fullfile(plot_save_path,strcat('Right_',cohPairs{i,1},'_',cohPairs{i,2})),'fig');
        saveas(gcf,fullfile(plot_save_path,strcat('Right_',cohPairs{i,1},'_',cohPairs{i,2})),'tiff');
    end
end
end

