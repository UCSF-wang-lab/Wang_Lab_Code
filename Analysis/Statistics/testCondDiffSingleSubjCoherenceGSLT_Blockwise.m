function testCondDiffSingleSubjCoherenceGSLT_Blockwise(varargin)
% Tests for differences in coherence between two task conditions (e.g. hits
% and misses), for a single subject
%
% INPUTS:  Optional
%               dataPath1[=]    The path to the folder containing the
%                               coherence data of the 1st condition.
%
%               dataPath2[=]    The path to the folder containing the
%                               coherence data of the 2nd condition.
%
%                cohPairs[=]    Pairs to compute coherence between. The
%                               list of pairs will be used across all
%                               files. Default is to do all possible
%                               unilateral pairs.
%               
%                    keys[=]    Cell array of keys to analyze. Default is
%                               to analyze all keys.
%
%                savePlot[=]    Boolean option to save the resulting merged plots.
%                               Default is false.
%
%                saveData[=]    Boolean option to save the coherence data.
%                               Default is false.
%
%              comparison[=]    String describing the comparison tested, 
%                               e.g. adaptVSmaladapt_right. It creates appropriate subfolders
%                               to save the plots in.
%
%                blocknum[=]    Block number to be analyzed
%
%                  unilat[=]    0 if bilat (default), -1 if left
%                               unilateral, 1 if right unilateral
%
%               stat_test[=]    it can be either 'perm' (permutation tests)
%                               or ttests ('ttest')
%           
%               mult_corr[=]    1 applies FDR correction, 0 for uncorrected
%               
%
%   Example call:
%           keys = {'key0','key1','key2','key3'};
%           testCondDiffSingleSubjCoherenceGSLT('keys',keys,'saveData',0,'savePlot',0,'comparison','adaptVSmaladapt_right','unilat',0);
%
% Date:     12/22/2023
% Author:   Eleni Patelaki (eleni.patelaki@ucsf.edu)

%% Parse optional arguments
for i = 1:2:nargin
    switch varargin{i}
        case 'dataPath1'
            dataPath1 = varargin{i+1};
        case 'dataPath2'
            dataPath2 = varargin{i+1};
        case 'keys'
            keys = varargin{i+1};
        case 'cohPairs'
            cohPairs = varargin{i+1};
        case 'savePlot'
            savePlot = varargin{i+1};
        case 'saveData'
            saveData = varargin{i+1};
        case 'comparison'
            comparison = varargin{i+1};
        case 'blocknum'
            blocknum = varargin{i+1};
        case 'unilat'
            unilat = varargin{i+1};
        case 'stat_test'
            stat_test = varargin{i+1};
        case 'mult_corr'
            mult_corr = varargin{i+1};
    end
end

%% Set to default values if not passed in by user
if ~exist('keys','var') || isempty(keys)
    keys = {'key0','key1','key2','key3'};
end

if ~exist('cohPairs','var') || isempty(cohPairs)
    cohPairs = nchoosek(keys,2);
end

if ~exist('comparison','var') || isempty(comparison)
    comparison = '';
end

if ~exist('blocknum','var') || isempty(blocknum)
    blocknum = '';
end

if ~exist('unilat','var') || isempty(unilat)
    unilat = 0;
end

if ~exist('saveData','var') || isempty(saveData)
    saveData = 0;
end

if ~exist('stat_test','var') || isempty(stat_test)
    stat_test = 'perm';
end

if ~exist('mult_corr','var') || isempty(mult_corr)
    mult_corr = 0;
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
    
    plot_save_path = fullfile(parent_dir,'Figures','GSLT','Coherence_figs_adaptive_condDiffs','Merged',comparison);
    if ~exist(plot_save_path, 'dir')
       mkdir(plot_save_path);
    end
end

%% Load the data
% Fetch & load coherence data for the first condition
if ~exist('dataPath1','var') || isempty(dataPath1)
    [~,dataPath1] = uigetfile('normCoh.mat','Select the data file of the 1st condition.');
end
data_cond1 = load(fullfile(dataPath1, 'normCoh.mat'));
fc_cond1 = load(fullfile(dataPath1, 'allfc.mat'));
avgTO{1} = load(fullfile(dataPath1, 'avgToeOff.mat'));
blocknum_cond1 = load(fullfile(dataPath1, 'allblockNums.mat'));

% Fetch & load coherence data for the second condition
if ~exist('dataPath2','var') || isempty(dataPath2)
    [~,dataPath2] = uigetfile('normCoh.mat','Select the data file of the 2nd condition.');
end
data_cond2 = load(fullfile(dataPath2, 'normCoh.mat'));
fc_cond2 = load(fullfile(dataPath2, 'allfc.mat'));
avgTO{2} = load(fullfile(dataPath2, 'avgToeOff.mat'));
blocknum_cond2 = load(fullfile(dataPath2, 'allblockNums.mat'));


%% Identify the statistically significant effects using permutation tests (FDR correction)
% Left
if unilat==0 || unilat==-1
    h_all_left = cell(length(data_cond1.normCoh.Left),max(blocknum_cond1.allblockNums.Left{1}));
    p_all_left = cell(length(data_cond1.normCoh.Left),max(blocknum_cond1.allblockNums.Left{1}));
    normCoh_block.Left = cell(length(data_cond1.normCoh.Left),max(blocknum_cond1.allblockNums.Left{1}),2);
    for i = 1:length(data_cond1.normCoh.Left)
        data_snip_cond1 = data_cond1.normCoh.Left{i};
        data_snip_cond2 = data_cond2.normCoh.Left{i};
        blocknum_cohPair_cond1 = blocknum_cond1.allblockNums.Left{i};
        blocknum_cohPair_cond2 = blocknum_cond2.allblockNums.Left{i};
        assert(size(data_snip_cond1,1)==size(data_snip_cond2,1));
        assert(size(data_snip_cond1,2)==size(data_snip_cond2,2));
        
        for j = 1:max(blocknum_cond1.allblockNums.Left{1})
            data_snip_cond1_block = data_snip_cond1(:,:,blocknum_cohPair_cond1==j); 
            data_snip_cond2_block = data_snip_cond2(:,:,blocknum_cohPair_cond2==j);

            h_all_left{i,j} = NaN(size(data_snip_cond1_block,1),size(data_snip_cond1_block,2));
            p_all_left{i,j} = NaN(size(data_snip_cond1_block,1),size(data_snip_cond1_block,2));
            for f = 1:size(data_snip_cond1_block,1)
                for t = 1:size(data_snip_cond1_block,2)
                    if strcmpi(stat_test,'perm')
                        n_perm = 500;
                        [p,~,~] = permutationTest(squeeze(data_snip_cond1_block(f,t,:)),squeeze(data_snip_cond2_block(f,t,:)),n_perm,'exact',1);
                        h=p<0.05;
                    elseif strcmpi(stat_test,'ttest')
                        [h,p] = ttest2(squeeze(data_snip_cond1_block(f,t,:)),squeeze(data_snip_cond2_block(f,t,:)));
                    end
    
                    if h==0
                        h =0.5;
                    end
    
                    p_all_left{i,j}(f,t) = p;
                    h_all_left{i,j}(f,t) = h;
    
                    clear h p
                end
    
                if mult_corr
                    % FDR
                    FDR = mafdr(p_all_left{i,j}(f,:));
        
                    for t=1:length(FDR)
                        if FDR(t)<0.05
                            h_fdr_new =0.5;
                        else
                            h_fdr_new = 1;
                        end
                        h_all_left{i,j}(f,t) = h_fdr_new;
    
                        clear h_fdr_new
                    end
        
                    clear FDR
                end
            end
            normCoh_block.Left{i,j,1} = data_snip_cond1_block;
            normCoh_block.Left{i,j,2} = data_snip_cond2_block;

            clear data_snip_cond1_block data_snip_cond2_block
        end
        clear data_snip_cond1 data_snip_cond2 blocknum_cohPair_cond1 blocknum_cohPair_cond2
    end
end

% Right
if unilat==0 || unilat==1
    h_all_right = cell(length(data_cond1.normCoh.Right),max(blocknum_cond1.allblockNums.Right{1}));
    p_all_right = cell(length(data_cond1.normCoh.Right),max(blocknum_cond1.allblockNums.Right{1}));
    normCoh_block.Right = cell(length(data_cond1.normCoh.Right),max(blocknum_cond1.allblockNums.Right{1}),2);
    for i = 1:length(data_cond1.normCoh.Right)
        data_snip_cond1 = data_cond1.normCoh.Right{i};
        data_snip_cond2 = data_cond2.normCoh.Right{i};
        blocknum_cohPair_cond1 = blocknum_cond1.allblockNums.Right{i};
        blocknum_cohPair_cond2 = blocknum_cond2.allblockNums.Right{i};
        assert(size(data_snip_cond1,1)==size(data_snip_cond2,1));
        assert(size(data_snip_cond1,2)==size(data_snip_cond2,2));

        for j = 1:max(blocknum_cond1.allblockNums.Right{1})
            data_snip_cond1_block = data_snip_cond1(:,:,blocknum_cohPair_cond1==j);
            data_snip_cond2_block = data_snip_cond2(:,:,blocknum_cohPair_cond2==j);

            h_all_right{i,j} = NaN(size(data_snip_cond1_block,1),size(data_snip_cond1_block,2));
            p_all_right{i,j} = NaN(size(data_snip_cond1_block,1),size(data_snip_cond1_block,2));
            for f = 1:size(data_snip_cond1,1)
                for t = 1:size(data_snip_cond1,2)
                    if strcmpi(stat_test,'perm')
                        n_perm = 500;
                        [p,~,~] = permutationTest(squeeze(data_snip_cond1_block(f,t,:)),squeeze(data_snip_cond2_block(f,t,:)),n_perm,'exact',1);
                        h=p<0.05;
                    elseif strcmpi(stat_test,'ttest')
                        [h,p] = ttest2(squeeze(data_snip_cond1_block(f,t,:)),squeeze(data_snip_cond2_block(f,t,:)));
                    end
    
                    if h==0
                        h =0.5;
                    end
    
                    p_all_right{i,j}(f,t) = p;
                    h_all_right{i,j}(f,t) = h;
    
                    clear h p
                end
    
                if mult_corr
                    % FDR
                    FDR = mafdr(p_all_right{i,j}(f,:));
        
                    for t=1:length(FDR)
                        if FDR(t)<0.05
                            h_fdr_new =0.5;
                        else
                            h_fdr_new = 1;
                        end
                        h_all_right{i,j}(f,t) = h_fdr_new;
    
                        clear h_fdr_new
                    end
        
                    clear FDR
                end
            end
            normCoh_block.Right{i,j,1} = data_snip_cond1_block;
            normCoh_block.Right{i,j,2} = data_snip_cond2_block;

            clear data_snip_cond1_block data_snip_cond2_block
        end
        clear data_snip_cond1 data_snip_cond2 blocknum_cohPair_cond1 blocknum_cohPair_cond2
    end
end

% Calculate grand averages
if unilat==0 || unilat==-1
    for i = 1:length(data_cond1.normCoh.Left)
        for j = 1:max(blocknum_cond1.allblockNums.Left{1})
            grandAvgCoh{1}.Left{i,j} = mean(normCoh_block.Left{i,j,1},3);
            grandAvgCoh{2}.Left{i,j} = mean(normCoh_block.Left{i,j,2},3);
        end
    end
end

if unilat==0 || unilat==1
    for i = 1:length(data_cond1.normCoh.Right)
        for j = 1:max(blocknum_cond1.allblockNums.Right{1})
            grandAvgCoh{1}.Right{i,j} = mean(normCoh_block.Right{i,j,1},3);
            grandAvgCoh{2}.Right{i,j} = mean(normCoh_block.Right{i,j,2},3);
        end
    end
end

%% Plotting
if unilat==0 || unilat==-1
    for i = 1:size(grandAvgCoh{1}.Left,1)
        for j = 1:max(blocknum_cond1.allblockNums.Left{1})
            for cond=1:length(grandAvgCoh)
                assert(sum(fc_cond1.allfc.Left{1}(:,:,1)-fc_cond2.allfc.Left{1}(:,:,1),'all')==0);
                freq_vec = squeeze(fc_cond1.allfc.Left{1}(:,:,1));
        
                figure;
                h = imagesc(1:100,log2(freq_vec),grandAvgCoh{cond}.Left{i,j}, 'Interpolation', 'bilinear');
                ax = gca;
                % ticks = logspace(log10(2.5),log10(50),5);
                ticks = logspace(log10(2.5),log10(70),5);
                ax.YTick = log2(ticks);
                ax.YTickLabel = sprintf('%1.2f\n',round(ticks,2));
                xline(avgTO{cond}.avgToeOff.Left{1},'LineStyle','--','LineWidth',1.5);
                % ylim([log2(2.5),log2(50)]);
                ylim([log2(2.5),log2(70)]);
                xticks([1,10,20,30,40,50,60,70,80,90,100]);
                xticklabels({'0','10','20','30','40','50','60','70','80','90','100'});
                colormap jet;
                
                % Apply transparency mask to non-significant p-values
                alphaMask = h_all_left{i,j};
                set(h, 'AlphaData', alphaMask);
                set(gca,'YDir','normal');
        
                ylabel('Frequency (Hz)');
                title({'Left';sprintf('%s to %s Grand Average Coherence',cohPairs{i,1},cohPairs{i,2})});
        
                % Save plots
                if savePlot    
                    if ~exist(fullfile(plot_save_path,strcat('Blockwise_',stat_test),strcat('Block ',num2str(j))),'dir')
                        mkdir(fullfile(plot_save_path,strcat('Blockwise_',stat_test),strcat('Block ',num2str(j))));
                    end

                    saveas(gcf,fullfile(plot_save_path,strcat('Blockwise_',stat_test),strcat('Block ',num2str(j)),strcat('Left_',cohPairs{i,1},'_',cohPairs{i,2},'_cond',num2str(cond))),'fig');
                    saveas(gcf,fullfile(plot_save_path,strcat('Blockwise_',stat_test),strcat('Block ',num2str(j)),strcat('Left_',cohPairs{i,1},'_',cohPairs{i,2},'_cond',num2str(cond))),'tiff');
                end
            end
        end
        close all
    end
end

if unilat==0 || unilat==1
    for i = 1:size(grandAvgCoh{1}.Right,1)
        for j = 1:max(blocknum_cond1.allblockNums.Right{1})
            for cond=1:length(grandAvgCoh)
                assert(sum(fc_cond1.allfc.Right{1}(:,:,1)-fc_cond2.allfc.Right{1}(:,:,1),'all')==0);
                freq_vec = squeeze(fc_cond1.allfc.Right{1}(:,:,1));
        
                figure;
                h = imagesc(1:100,log2(freq_vec),grandAvgCoh{cond}.Right{i,j}, 'Interpolation', 'bilinear');
                ax = gca;
                % ticks = logspace(log10(2.5),log10(50),5);
                ticks = logspace(log10(2.5),log10(70),5);
                ax.YTick = log2(ticks);
                ax.YTickLabel = sprintf('%1.2f\n',round(ticks,2));
                xline(avgTO{cond}.avgToeOff.Right{1},'LineStyle','--','LineWidth',1.5);
                % ylim([log2(2.5),log2(50)]);
                ylim([log2(2.5),log2(70)]);
                xticks([1,10,20,30,40,50,60,70,80,90,100]);
                xticklabels({'0','10','20','30','40','50','60','70','80','90','100'});
                colormap jet;
        
                % Apply transparency mask to non-significant p-values
                alphaMask = h_all_right{i,j};
                set(h, 'AlphaData', alphaMask);
                set(gca,'YDir','normal') 
        
                ylabel('Frequency (Hz)');
                title({'Right';sprintf('%s to %s Grand Average Coherence',cohPairs{i,1},cohPairs{i,2})});
        
                % Save plots
                if savePlot
                    if ~exist(fullfile(plot_save_path,strcat('Blockwise_',stat_test),strcat('Block ',num2str(j))),'dir')
                        mkdir(fullfile(plot_save_path,strcat('Blockwise_',stat_test),strcat('Block ',num2str(j))));
                    end

                    saveas(gcf,fullfile(plot_save_path,strcat('Blockwise_',stat_test),strcat('Block ',num2str(j)),strcat('Right_',cohPairs{i,1},'_',cohPairs{i,2},'_cond',num2str(cond))),'fig');
                    saveas(gcf,fullfile(plot_save_path,strcat('Blockwise_',stat_test),strcat('Block ',num2str(j)),strcat('Right_',cohPairs{i,1},'_',cohPairs{i,2},'_cond',num2str(cond))),'tiff');
                end
            end
        end
        close all
    end
end

end

