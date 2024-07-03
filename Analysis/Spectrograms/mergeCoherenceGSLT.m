function mergeCoherenceGSLT(varargin)
% When multiple aligned datasets exist within the same experimental session, 
% this function merges the multiple occuring coherence decompositions.
%
% INPUTS:  Optional
%               dataPath[=]     The path to the folder containing the 
%                               coherence decompositions.
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
%               saveData[=]     Boolean option to save the coherence data.
%                               Default is false.
%
%               condition[=]    String describing the condition analyzed, 
%                               e.g. OFFmed_ONdbs. It creates appropriate subfolders
%                               to save the plots in.
%
%               unilat[=]       0 if bilat (default), -1 if left
%                               unilateral, 1 if right unilateral
%
%
% Example call:
%           keys = {'key0','key1','key2','key3'};
%           mergeCoherenceGSLT('keys',keys,'savePlot',1,'condition','OFFmed_ONdbs_RightAdaptive')
%
% Date:     10/24/2023
% Author:   Eleni Patelaki (eleni.patelaki@ucsf.edu)

%% Parse optional arguments
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
        case 'saveData'
            saveData = varargin{i+1};
        case 'condition'
            condition = varargin{i+1};
        case 'unilat'
            unilat = varargin{i+1};
    end
end

%% Set to default values if not passed in by user
if ~exist('keys','var') || isempty(keys)
    keys = {'key0','key1','key2','key3'};
end

if ~exist('cohPairs','var') || isempty(cohPairs)
    cohPairs = nchoosek(keys,2);
end

if ~exist('condition','var') || isempty(condition)
    condition = '';
end

if ~exist('unilat','var') || isempty(unilat)
    unilat = 0;
end

if ~exist('saveData','var') || isempty(saveData)
    saveData = 0;
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
    
    plot_save_path = fullfile(parent_dir,'Figures','GSLT','Coherence_figs_adaptive_masked','Merged',condition);
    if ~exist(plot_save_path, 'dir')
       mkdir(plot_save_path);
    end
end

%% Fetch necessary files
% All coherence datafiles
if ~exist('dataPath','var') || isempty(dataPath)
    [cohdatafiles,dataPath] = uigetfile('step_coh*.mat','Select file(s)', 'MultiSelect', 'on');
else
    cohdatafiles = dir(fullfile(dataPath, 'step_coh*.mat'));
end
if ~iscell(cohdatafiles)==1
    cohdatafiles = {cohdatafiles};
end

% All frequency & normalization datafiles that correspond to the selected coherence data files
cohdatafiles_unnormDur = cellfun(@(x) strrep(x,'step','step_unnormDur'),cohdatafiles,'UniformOutput',false);
freqfnames = cellfun(@(x) strrep(x,'step','fc'),cohdatafiles,'UniformOutput',false);
normfnames = cellfun(@(x) strrep(x,'step_coh','normalization'),cohdatafiles,'UniformOutput',false);
lenNormfnames = cellfun(@(x) strrep(x,'step_coh','lenNorm'),cohdatafiles,'UniformOutput',false);
tofnames = cellfun(@(x) strrep(x,'step_coh','toeoff'),cohdatafiles,'UniformOutput',false);
blockNumfnames = cellfun(@(x) strrep(x,'step_coh','blockNum'),cohdatafiles,'UniformOutput',false);
tofnames_unnormDur = cellfun(@(x) strrep(x,'step_coh','toeoff_unnormDur'),cohdatafiles,'UniformOutput',false);

%% Ensure that the number of the frequency & normalization datafiles is equal to that of the coherence datafiles
assert((numel(cohdatafiles)==numel(freqfnames)) & (numel(freqfnames)==numel(normfnames)) & (numel(freqfnames)==numel(tofnames)));

%% Initialize the aggregated coherence, normalization structs and gait cycle number structures
allsteps.Left = cell(1,size(cohPairs,1));
allsteps.Right = cell(1,size(cohPairs,1));
allstepsUD.Left = cell(1,size(cohPairs,1));
allstepsUD.Right = cell(1,size(cohPairs,1));
allfc.Left = cell(1,size(cohPairs,1));
allfc.Right = cell(1,size(cohPairs,1));
allnorm.Left = cell(1,size(cohPairs,1));
allnorm.Right = cell(1,size(cohPairs,1));
allLenNorm.Left = cell(1,size(cohPairs,1));
allLenNorm.Right = cell(1,size(cohPairs,1));
allstepNums.Left = cell(1,size(cohPairs,1));
allstepNums.Right = cell(1,size(cohPairs,1));
alltoeOffs.Left = cell(1,size(cohPairs,1));
alltoeOffs.Right = cell(1,size(cohPairs,1));
alltoeOffsUD.Left = cell(1,size(cohPairs,1));
alltoeOffsUD.Right = cell(1,size(cohPairs,1));
allblockNums.Left = cell(1,size(cohPairs,1));
allblockNums.Right = cell(1,size(cohPairs,1));


%% Load & merge the coherence & normalization datafiles
for i = 1:numel(cohdatafiles)
    if isstruct(cohdatafiles)
        load(fullfile(dataPath,getfield(cohdatafiles(i),'name')));
        load(fullfile(dataPath,getfield(cohdatafiles_unnormDur(i),'name')));
    else
        load(fullfile(dataPath,cohdatafiles{i}));
        load(fullfile(dataPath,cohdatafiles_unnormDur{i}));
    end
    load(fullfile(dataPath,freqfnames{i}));
    load(fullfile(dataPath,normfnames{i}));
    load(fullfile(dataPath,lenNormfnames{i}));
    load(fullfile(dataPath,tofnames{i}));
    load(fullfile(dataPath,blockNumfnames{i}));
    load(fullfile(dataPath,tofnames_unnormDur{i}));

    if sum(step.Left{1},'all')~=0 % && sum(step.Right{1},'all')~=0
        % Get the number of gait cycles in each part
        stepNums.Left = cellfun(@(x) size(x,3),step.Left,'UniformOutput',false);
        stepNums.Right = cellfun(@(x) size(x,3),step.Right,'UniformOutput',false);
    
        allsteps.Left = cellfun(@(x,y) cat(3,x,y),allsteps.Left,step.Left,'UniformOutput',false);
        allsteps.Right = cellfun(@(x,y) cat(3,x,y),allsteps.Right,step.Right,'UniformOutput',false);
        allstepsUD.Left = cellfun(@(x,y) [x,y],allstepsUD.Left,stepUnnormDur.Left,'UniformOutput',false);
        allstepsUD.Right = cellfun(@(x,y) [x,y],allstepsUD.Right,stepUnnormDur.Right,'UniformOutput',false);
        allfc.Left = cellfun(@(x,y) cat(3,x,y),allfc.Left,fc.Left,'UniformOutput',false);
        allfc.Right = cellfun(@(x,y) cat(3,x,y),allfc.Right,fc.Right,'UniformOutput',false);
        allnorm.Left = cellfun(@(x,y) cat(3,x,y),allnorm.Left,normalization.Left,'UniformOutput',false);
        allnorm.Right = cellfun(@(x,y) cat(3,x,y),allnorm.Right,normalization.Right,'UniformOutput',false);
        allLenNorm.Left = cellfun(@(x,y) cat(1,x,y),allLenNorm.Left,lenNorm.Left,'UniformOutput',false);
        allLenNorm.Right = cellfun(@(x,y) cat(1,x,y),allLenNorm.Right,lenNorm.Right,'UniformOutput',false);
        allstepNums.Left = cellfun(@(x,y) cat(2,x,y),allstepNums.Left,stepNums.Left,'UniformOutput',false);
        allstepNums.Right = cellfun(@(x,y) cat(2,x,y),allstepNums.Right,stepNums.Right,'UniformOutput',false);
        alltoeOffs.Left = cellfun(@(x,y) cat(2,x,y),alltoeOffs.Left,toeOff.Left,'UniformOutput',false);
        alltoeOffs.Right = cellfun(@(x,y) cat(2,x,y),alltoeOffs.Right,toeOff.Right,'UniformOutput',false);
        alltoeOffsUD.Left = cellfun(@(x,y) cat(2,x,y),alltoeOffsUD.Left,toeOffUnnormDur.Left,'UniformOutput',false);
        alltoeOffsUD.Right = cellfun(@(x,y) cat(2,x,y),alltoeOffsUD.Right,toeOffUnnormDur.Right,'UniformOutput',false);
        allblockNums.Left = cellfun(@(x,y) cat(2,x,y),allblockNums.Left,blockNum.Left,'UniformOutput',false);
        allblockNums.Right = cellfun(@(x,y) cat(2,x,y),allblockNums.Right,blockNum.Right,'UniformOutput',false);

        clear step stepUnnormDur stepNums normalization toeOff toeOffUnnormDur lenNorm
    end
end


%% Calculate the weighted average of the normalizations matrices & normalize
for i = 1:size(cohPairs,1)
    for j = 1:numel(cohdatafiles)
        if j==1
            if unilat==0 || unilat==-1
                sumNormLeft =  allLenNorm.Left{i}(j)*allnorm.Left{i}(:,:,j);
            end
            if unilat==0 || unilat==1
                sumNormRight =  allLenNorm.Right{i}(j)*allnorm.Right{i}(:,:,j);
            end
        else
            if unilat==0 || unilat==-1
                sumNormLeft =  sumNormLeft + allLenNorm.Left{i}(j)*allnorm.Left{i}(:,:,j);
            end
            if unilat==0 || unilat==1
                sumNormRight =  sumNormRight + allLenNorm.Right{i}(j)*allnorm.Right{i}(:,:,j);
            end
        end
    end
    
    if unilat==0 || unilat==-1
        avgnormLeft = sumNormLeft/sum(allLenNorm.Left{i},'all');

        % Apply the normalization (z-scoring)
        normCoh.Left{i} = (allsteps.Left{i}-avgnormLeft(:,1))./avgnormLeft(:,2);
        normCohUD.Left{i} = cellfun(@(x) (x-avgnormLeft(:,1))./avgnormLeft(:,2),allstepsUD.Left{i},'UniformOutput',false);

        % Store the average normalization matrix
        avgnorm.Left{i} = avgnormLeft;
    end
    if unilat==0 || unilat==1
        avgnormRight = sumNormRight/sum(allLenNorm.Right{i},'all');

        % Apply the normalization (z-scoring)
        normCoh.Right{i} = (allsteps.Right{i}-avgnormRight(:,1))./avgnormRight(:,2);
        normCohUD.Right{i} = cellfun(@(x) (x-avgnormRight(:,1))./avgnormRight(:,2),allstepsUD.Right{i},'UniformOutput',false);

        % Store the average normalization matrix
        avgnorm.Right{i} = avgnormRight;
    end
    
    clear avgnormLeft avgnormRight sumNormLeft sumNormRight
end

%% Identify the statistically significant effects using FDR correction
% Left
if unilat==0 || unilat==-1
    h_all_left = cell(1,length(normCoh.Left));
    p_all_left = cell(1,length(normCoh.Left));
    for i = 1:length(normCoh.Left)
        data_snip = normCoh.Left{i};
        h_all_left{i} = NaN(size(data_snip,1),size(data_snip,2));
        p_all_left{i} = NaN(size(data_snip,1),size(data_snip,2));
        for f = 1:size(data_snip,1)
            for t = 1:size(data_snip,2)
                [h,p] = ttest(data_snip(f,t,:));
                if h==0
                    h_new =0.5;
                elseif h==1
                    h_new = 1;
                end
                h_all_left{i}(f,t) = h_new;
                p_all_left{i}(f,t) = p;

                clear h h_new p
            end
        end
    end
end

% Right
if unilat==0 || unilat==1
    h_all_right = cell(1,length(normCoh.Right));
    p_all_right = cell(1,length(normCoh.Right));
    for i = 1:length(normCoh.Right)
        data_snip = normCoh.Right{i};
        h_all_right{i} = NaN(size(data_snip,1),size(data_snip,2));
        p_all_right{i} = NaN(size(data_snip,1),size(data_snip,2));
        for f = 1:size(data_snip,1)
            for t = 1:size(data_snip,2)
                [h,p] = ttest(data_snip(f,t,:));
                if h==0
                    h_new =0.5;
                elseif h==1
                    h_new = 1;
                end
                h_all_right{i}(f,t) = h_new;
                p_all_right{i}(f,t) = p;

                clear h h_new p
            end
        end
    end
end

% Calculate grand averages
if unilat==0 || unilat==-1
    for i = 1:length(normCoh.Left)
        grandAvgCoh.Left{i} = mean(normCoh.Left{i},3);
        avgToeOff.Left{i} = mean(alltoeOffs.Left{i},'omitnan');
    end
end

if unilat==0 || unilat==1
    for i = 1:length(normCoh.Right)
        grandAvgCoh.Right{i} = mean(normCoh.Right{i},3);
        avgToeOff.Right{i} = mean(alltoeOffs.Right{i},'omitnan');
    end
end

% Save the grand average data (optional)
if saveData
    save(fullfile(dataPath,'grandAvgCoh.mat'),'grandAvgCoh','-v7.3');
    save(fullfile(dataPath,'normCoh.mat'),'normCoh','-v7.3');
    save(fullfile(dataPath,'normCohUD.mat'),'normCohUD','-v7.3');
    save(fullfile(dataPath,'allfc.mat'),'allfc','-v7.3');
    save(fullfile(dataPath,'alltoeOffs.mat'),'alltoeOffs','-v7.3');
    save(fullfile(dataPath,'allblockNums.mat'),'allblockNums','-v7.3');
    save(fullfile(dataPath,'alltoeOffsUD.mat'),'alltoeOffsUD','-v7.3');
    save(fullfile(dataPath,'avgToeOff.mat'),'avgToeOff','-v7.3');
    save(fullfile(dataPath,'avgnorm.mat'),'avgnorm','-v7.3');
end

%% Plotting
if unilat==0 || unilat==-1
    for i = 1:length(grandAvgCoh.Left)
        freq_vec = squeeze(allfc.Left{1}(:,:,1));

        figure;
        h = imagesc(1:100,log2(freq_vec),grandAvgCoh.Left{i}, 'Interpolation', 'bilinear');
        ax = gca;
        ticks = logspace(log10(2.5),log10(50),5);
        ax.YTick = log2(ticks);
        ax.YTickLabel = sprintf('%1.2f\n',round(ticks,2));
        xline(avgToeOff.Left{i},'LineStyle','--','LineWidth',1.5);
        ylim([log2(2.5),log2(50)]);
        xticks([1,10,20,30,40,50,60,70,80,90,100]);
        xticklabels({'0','10','20','30','40','50','60','70','80','90','100'});
        colormap jet;
        
        % Apply transparency mask to non-significant p-values
        alphaMask = h_all_left{i};
        set(h, 'AlphaData', alphaMask);
        set(gca,'YDir','normal') 

        ylabel('Frequency (Hz)');
        title({'Left';sprintf('%s to %s Grand Average Coherence',cohPairs{i,1},cohPairs{i,2})});

        % Save plots
        if savePlot    
            saveas(gcf,fullfile(plot_save_path,strcat('Left_',cohPairs{i,1},'_',cohPairs{i,2})),'fig');
            saveas(gcf,fullfile(plot_save_path,strcat('Left_',cohPairs{i,1},'_',cohPairs{i,2})),'tiff');
        end
    end
end

if unilat==0 || unilat==1
    for i = 1:length(grandAvgCoh.Right)
        freq_vec = squeeze(allfc.Left{1}(:,:,1));

        figure;
        h = imagesc(1:100,log2(freq_vec),grandAvgCoh.Right{i}, 'Interpolation', 'bilinear');
        ax = gca;
        ticks = logspace(log10(2.5),log10(50),5);
        ax.YTick = log2(ticks);
        ax.YTickLabel = sprintf('%1.2f\n',round(ticks,2));
        xline(avgToeOff.Right{i},'LineStyle','--','LineWidth',1.5);
        ylim([log2(2.5),log2(50)]);
        xticks([1,10,20,30,40,50,60,70,80,90,100]);
        xticklabels({'0','10','20','30','40','50','60','70','80','90','100'});
        colormap jet;

        % Apply transparency mask to non-significant p-values
        alphaMask = h_all_right{i};
        set(h, 'AlphaData', alphaMask); 
        set(gca,'YDir','normal') 

        ylabel('Frequency (Hz)');
        title({'Right';sprintf('%s to %s Grand Average Coherence',cohPairs{i,1},cohPairs{i,2})});

        % Save plots
        if savePlot
            saveas(gcf,fullfile(plot_save_path,strcat('Right_',cohPairs{i,1},'_',cohPairs{i,2})),'fig');
            saveas(gcf,fullfile(plot_save_path,strcat('Right_',cohPairs{i,1},'_',cohPairs{i,2})),'tiff');
        end
    end
end

end

