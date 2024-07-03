function testCondDiffSingleSubjCoherenceGSLT_StepLenControl(varargin)
% Tests for differences in coherence between two task conditions (e.g. hits
% and misses), for a single subject
%
% INPUTS:  Optional
%               dataPath[=]    The path to the folder containing the
%                               coherence data of the 1st condition.
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
%                  unilat[=]    0 if bilat (default), -1 if left
%                               unilateral, 1 if right unilateral
%
%               stat_test[=]    corr, for now
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
        case 'dataPath'
            dataPath = varargin{i+1};
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

if ~exist('unilat','var') || isempty(unilat)
    unilat = 0;
end

if ~exist('saveData','var') || isempty(saveData)
    saveData = 0;
end

if ~exist('stat_test','var') || isempty(stat_test)
    stat_test = 'corr';
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
    
    plot_save_path = fullfile(parent_dir,'Figures','GSLT','Coherence_figs_adaptive_condDiffs','Merged_StepLenControl',comparison);
    if ~exist(plot_save_path, 'dir')
       mkdir(plot_save_path);
    end
end

%% Load the data
% Fetch & load coherence data for the condition to be analyzed
if ~exist('dataPath','var') || isempty(dataPath)
    [~,dataPath] = uigetfile('normCoh.mat','Select the data file of the 1st condition.');
end
data_cond = load(fullfile(dataPath, 'normCoh.mat'));
fc_cond = load(fullfile(dataPath, 'allfc.mat'));
avgTO = load(fullfile(dataPath, 'avgToeOff.mat'));
stepMod = load(fullfile(dataPath, 'allStepMod.mat'));

%% Identify the statistically significant effects using permutation tests (FDR correction)
% Left
if unilat==0 || unilat==-1
    h_all_left = cell(1,length(data_cond.normCoh.Left));
    p_all_left = cell(1,length(data_cond.normCoh.Left));
    for i = 1:length(data_cond.normCoh.Left)
        data_snip_cond = data_cond.normCoh.Left{i};
        h_all_left{i} = NaN(size(data_snip_cond,1),size(data_snip_cond,2));
        p_all_left{i} = NaN(size(data_snip_cond,1),size(data_snip_cond,2));
        for f = 1:size(data_snip_cond,1)
            for t = 1:size(data_snip_cond,2)
                if strcmpi(stat_test,'corr')
                    [rho,p] = corr(squeeze(data_snip_cond(f,t,:)),stepMod.allStepMod.Left{i}');
                    h=p<0.05;
                else
                    % TODO: Add alternative stat. testing options here
                end

                if h==0
                    h =0.5;
                end
                
                rho_all_left{i}(f,t) = rho;
                p_all_left{i}(f,t) = p;
                h_all_left{i}(f,t) = h;

                clear h p rho
            end

            if mult_corr
                % FDR
                FDR = mafdr(p_all_left{i}(f,:));
    
                for t=1:length(FDR)
                    if FDR(t)<0.05
                        h_fdr_new =0.5;
                    else
                        h_fdr_new = 1;
                    end
                    h_all_left{i}(f,t) = h_fdr_new;

                    clear h_fdr_new
                end
    
                clear FDR
            end
        end
        clear data_snip_cond data_snip_cond2
    end
end

% Right
if unilat==0 || unilat==1
    h_all_right = cell(1,length(data_cond.normCoh.Right));
    p_all_right = cell(1,length(data_cond.normCoh.Right));
    for i = 1:length(data_cond.normCoh.Right)
        data_snip_cond = data_cond.normCoh.Right{i};
        h_all_right{i} = NaN(size(data_snip_cond,1),size(data_snip_cond,2));
        p_all_right{i} = NaN(size(data_snip_cond,1),size(data_snip_cond,2));
        for f = 1:size(data_snip_cond,1)
            for t = 1:size(data_snip_cond,2)
                if strcmpi(stat_test,'corr')
                    [rho,p] = corr(squeeze(data_snip_cond(f,t,:)),stepMod.allStepMod.Right{i}');
                    h=p<0.05;
                else
                     % TODO: Add alternative stat. testing options here
                end

                if h==0
                    h =0.5;
                end

                rho_all_left{i}(f,t) = rho;
                p_all_right{i}(f,t) = p;
                h_all_right{i}(f,t) = h;

                clear h p rho
            end

            if mult_corr
                % FDR
                FDR = mafdr(p_all_right{i}(f,:));
    
                for t=1:length(FDR)
                    if FDR(t)<0.05
                        h_fdr_new =0.5;
                    else
                        h_fdr_new = 1;
                    end
                    h_all_right{i}(f,t) = h_fdr_new;

                    clear h_fdr_new
                end
    
                clear FDR
            end
        end
        clear data_snip_cond data_snip_cond2
    end
end

% Calculate grand averages
if unilat==0 || unilat==-1
    for i = 1:length(data_cond.normCoh.Left)
        grandAvgCoh{1}.Left{i} = mean(data_cond.normCoh.Left{i},3);
    end
end

if unilat==0 || unilat==1
    for i = 1:length(data_cond.normCoh.Right)
        grandAvgCoh{1}.Right{i} = mean(data_cond.normCoh.Right{i},3);
    end
end

%% Plotting
if unilat==0 || unilat==-1
    for i = 1:length(grandAvgCoh{1}.Left)
        for cond=1:length(grandAvgCoh)
            freq_vec = squeeze(fc_cond.allfc.Left{1}(:,:,1));
    
            figure;
            h = imagesc(1:100,log2(freq_vec),grandAvgCoh{cond}.Left{i}, 'Interpolation', 'bilinear');
            ax = gca;
            % ticks = logspace(log10(2.5),log10(50),5);
            ticks = logspace(log10(2.5),log10(70),5);
            ax.YTick = log2(ticks);
            ax.YTickLabel = sprintf('%1.2f\n',round(ticks,2));
            xline(avgTO.avgToeOff.Left{1},'LineStyle','--','LineWidth',1.5);
            % ylim([log2(2.5),log2(50)]);
            ylim([log2(2.5),log2(70)]);
            xticks([1,10,20,30,40,50,60,70,80,90,100]);
            xticklabels({'0','10','20','30','40','50','60','70','80','90','100'});
            colormap jet;
            
            % Apply transparency mask to non-significant p-values
            alphaMask = h_all_left{i};
            set(h, 'AlphaData', alphaMask);
            set(gca,'YDir','normal');
    
            ylabel('Frequency (Hz)');
            title({'Left';sprintf('%s to %s Grand Average Coherence',cohPairs{i,1},cohPairs{i,2})});

            figure;
            h = imagesc(1:100,log2(freq_vec),rho_all_left{i}, 'Interpolation', 'bilinear');
            ax = gca;
            % ticks = logspace(log10(2.5),log10(50),5);
            ticks = logspace(log10(2.5),log10(70),5);
            ax.YTick = log2(ticks);
            ax.YTickLabel = sprintf('%1.2f\n',round(ticks,2));
            xline(avgTO.avgToeOff.Left{1},'LineStyle','--','LineWidth',1.5);
            % ylim([log2(2.5),log2(50)]);
            ylim([log2(2.5),log2(70)]);
            xticks([1,10,20,30,40,50,60,70,80,90,100]);
            xticklabels({'0','10','20','30','40','50','60','70','80','90','100'});
            colormap jet;
            
            % Apply transparency mask to non-significant p-values
            alphaMask = h_all_left{i};
            set(h, 'AlphaData', alphaMask);
            set(gca,'YDir','normal');
    
            ylabel('Frequency (Hz)');
            title({'Left';sprintf('%s to %s Grand Average Coherence',cohPairs{i,1},cohPairs{i,2})});
    
            % Save plots
            if savePlot    
                saveas(gcf,fullfile(plot_save_path,strcat('Left_',cohPairs{i,1},'_',cohPairs{i,2},'_cond',num2str(cond))),'fig');
                saveas(gcf,fullfile(plot_save_path,strcat('Left_',cohPairs{i,1},'_',cohPairs{i,2},'_cond',num2str(cond))),'tiff');
            end
        end
    end
end

if unilat==0 || unilat==1
    for i = 1:length(grandAvgCoh{1}.Right)
        for cond=1:length(grandAvgCoh)
            freq_vec = squeeze(fc_cond.allfc.Right{1}(:,:,1));
    
            figure;
            h = imagesc(1:100,log2(freq_vec),grandAvgCoh{cond}.Right{i}, 'Interpolation', 'bilinear');
            ax = gca;
            % ticks = logspace(log10(2.5),log10(50),5);
            ticks = logspace(log10(2.5),log10(70),5);
            ax.YTick = log2(ticks);
            ax.YTickLabel = sprintf('%1.2f\n',round(ticks,2));
            xline(avgTO.avgToeOff.Right{1},'LineStyle','--','LineWidth',1.5);
            % ylim([log2(2.5),log2(50)]);
            ylim([log2(2.5),log2(70)]);
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
                saveas(gcf,fullfile(plot_save_path,strcat('Right_',cohPairs{i,1},'_',cohPairs{i,2},'_cond',num2str(cond))),'fig');
                saveas(gcf,fullfile(plot_save_path,strcat('Right_',cohPairs{i,1},'_',cohPairs{i,2},'_cond',num2str(cond))),'tiff');
            end
        end
    end
end

end

