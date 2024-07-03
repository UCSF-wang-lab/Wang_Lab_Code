function testCondDiffMultSubjCoherenceGSLT_StepLenControl(subjList,cond,corrType,varargin)
% Tests for differences in coherence between two task conditions (e.g. hits
% and misses), for multiple subjects
%
% INPUTS:  Required
%                subjList[=]    Vector containing the patients numbers to
%                               be analyzed, e.g. [2,3,5] for gRCS02,
%                               grRCS03 and gRCS05.
%
%                   cond1[=]    Condition 1 string. Must be the same as the
%                               corresponding data folder name, 
%                               e.g. ONmed_ONdbs_RightAdaptive
%
%                corrType[=]    String describing the correlation tested, 
%                               e.g. ONmed_AdaptRight_StepLenControl. It creates appropriate subfolders
%                               to save the plots in.
%
%           Optional
%
%                    keys[=]    Cell array of keys to analyze. Default is
%                               to analyze all keys.
%           
%                cohPairs[=]    Pairs of keys to compute coherence between.
%                               The list of pairs will be used across all
%                               files. Default is to do all possible
%                               unilateral pairs.
%
%                savePlot[=]    Boolean option to save the resulting merged plots.
%                               Default is false.
%
%                saveData[=]    Boolean option to save the coherence data.
%                               Default is false.
%
%                  unilat[=]    0 if bilat (default), -1 if left
%                               unilateral, 1 if right unilateral. Default
%                               is 0.
%
%               mult_corr[=]    0 for no correction, 1 for FDR correction
%
%               stat_test[=]    'corr' for correlation (only option for now)
%
%                dbs_cond[=]    'dbsopt' for data recorded during the dbs-optimized phase
%                               'pre-prog' for data recorded during the pre-programming phase
%
%   Example call:
%           keys = {'key0','key1','key2','key3'};
%           testCondDiffMultSubjCoherenceGSLT([2,3,4],'ONmed_ONdbs_LeftAdaptive','ONmed_ONdbs_LeftMaladaptive','ONmed_ONdbs_AdaptvsMaladapt_Left','saveData',0,'savePlot',1,'unilat',1,'stat_test','ttest','dbs_cond','dbsopt');
%
% Date:     3/24/2024
% Author:   Eleni Patelaki (eleni.patelaki@ucsf.edu)

%% Parse optional arguments
for i = 1:2:nargin-4
    switch varargin{i}
        case 'keys'
            keys = varargin{i+1};
        case 'cohPairs'
            cohPairs = varargin{i+1};
        case 'savePlot'
            savePlot = varargin{i+1};
        case 'saveData'
            saveData = varargin{i+1};
        case 'unilat'
            unilat = varargin{i+1};
        case 'mult_corr'
            mult_corr = varargin{i+1};
        case 'stat_test'
            stat_test = varargin{i+1};
        case 'dbs_cond'
            dbs_cond = varargin{i+1};
    end
end

%% Set to default values if not passed in by user
if ~exist('keys','var') || isempty(keys)
    keys = {'key0','key1','key2','key3'};
end

if ~exist('cohPairs','var') || isempty(cohPairs)
    cohPairs = nchoosek(keys,2);
end

if ~exist('unilat','var') || isempty(unilat)
    unilat = 0;
end

if ~exist('saveData','var') || isempty(saveData)
    saveData = 0;
else
    data_save_path = fullfile('X:\Patient Data\RC+S Data\gait_RCS aggregate data\Eleni\GSLT_coherence\Data',[corrType,'_',dbs_cond,'_',stat_test,'_subs_3_4_5_HigherFreq_upto70_v1']);
    if ~exist(data_save_path, 'dir')
       mkdir(data_save_path);
    end
end

if ~exist('mult_corr','var') || isempty(mult_corr)
    mult_corr = 0;
end

if ~exist('stat_test','var') || isempty(stat_test)
    stat_test = 'ttest';
end

if ~exist('dbs_cond','var') || isempty(dbs_cond)
    dbs_cond = 'dbsopt';
end

if ~exist('savePlot','var') || isempty(savePlot)
    savePlot = 0;   % Does not save plot by default.
else
    plot_save_path = fullfile('X:\Patient Data\RC+S Data\gait_RCS aggregate data\Eleni\GSLT_coherence\Figures',[corrType,'_',dbs_cond,'_',stat_test,'_subs_2_3_4_RHO_HigherFreq_upto70_v1']);
    if ~exist(plot_save_path, 'dir')
       mkdir(plot_save_path);
    end
end

%% Load the data iteratively for each patient
basePath = 'X:\Patient Data\RC+S Data';
allsteps{1}.Left = cell(1,size(cohPairs,1));
allsteps{1}.Right = cell(1,size(cohPairs,1));
allTO{1}.Left = cell(1,size(cohPairs,1));
allTO{1}.Right = cell(1,size(cohPairs,1));
allStepMod{1}.Left = cell(1,size(cohPairs,1));
allStepMod{1}.Right = cell(1,size(cohPairs,1));
allfc.Left = cell(1,size(cohPairs,1));
allfc.Right = cell(1,size(cohPairs,1));
subjVar{1} = [];

for sub = 1:length(subjList)
    dirList = dir(fullfile(basePath,sprintf('gait_RCS_%02d',subjList(sub))));
    fileIdx = find(contains(lower({dirList.name}),dbs_cond));

    if length(fileIdx)==0
        error('No appropriate folder was found. Check')
    elseif length(fileIdx)==1
        fileIdx_final = fileIdx;
    else
        warning('More than 1 folder matches for subject gait_RCS_%02d'.',subjList(sub));
        if subjList(sub)==2
            cnt = 0;
            for folderMatch = 1:length(fileIdx)
                if contains(lower(dirList(fileIdx(folderMatch)).name),'side')
                    fileIdx_final = fileIdx(folderMatch);

                    cnt = cnt + 1;
                    if cnt == 2
                        error('Found 2 matches.');
                    end
                end
            end
        end
    end
    assert(length(fileIdx_final)==1);

    % Load data for condition 1
    data_path = fullfile(basePath,sprintf('gait_RCS_%02d',subjList(sub)),dirList(fileIdx_final).name,'Data\Analysis Data\GSLT\Coherence_data_adaptive',cond);
    data = load(fullfile(data_path, 'normCoh.mat'));
    fc = load(fullfile(data_path, 'allfc.mat'));
    toeOff = load(fullfile(data_path, 'alltoeOffs.mat'));
    stepMod = load(fullfile(data_path, 'allStepMod.mat'));

    % Left
    if unilat==0 || unilat==-1
        allsteps{1}.Left = cellfun(@(x,y) cat(3,x,y),allsteps{1}.Left,data.normCoh.Left,'UniformOutput',false);
        allTO{1}.Left = cellfun(@(x,y) cat(2,x,y),allTO{1}.Left,toeOff.alltoeOffs.Left,'UniformOutput',false);
        allStepMod{1}.Left = cellfun(@(x,y) cat(2,x,y),allStepMod{1}.Left,stepMod.allStepMod.Left,'UniformOutput',false);
        subjVar{1} = [subjVar{1}; sub*ones(size(data.normCoh.Left{1},3),1)];

        allfc.Left = fc.allfc.Left;
    end

    % Right
    if unilat==0 || unilat==1
        allsteps{1}.Right = cellfun(@(x,y) cat(3,x,y),allsteps{1}.Right,data.normCoh.Right,'UniformOutput',false);
        allTO{1}.Right = cellfun(@(x,y) cat(2,x,y),allTO{1}.Right,toeOff.alltoeOffs.Right,'UniformOutput',false);
        allStepMod{1}.Right = cellfun(@(x,y) cat(2,x,y),allStepMod{1}.Right,stepMod.allStepMod.Right,'UniformOutput',false);
        if unilat==1
            subjVar{1} = [subjVar{1}; sub*ones(size(data.normCoh.Right{1},3),1)];
        end

        allfc.Right = fc.allfc.Right;
    end

end

%% Identify the statistically significant effects using permutation tests (FDR correction)
% Left
if unilat==0 || unilat==-1
    h_all.Left = cell(1,length(allsteps{1}.Left));
    p_all.Left = cell(1,length(allsteps{1}.Left));
    rho_all.Left = cell(1,length(allsteps{1}.Left));
    subjVar_corr = subjVar{1};
    stepMod_corr = allStepMod{1}.Left;
    for i = 1:length(allsteps{1}.Left)
        data_snip = allsteps{1}.Left{i};
        h_all.Left{i} = NaN(size(data_snip,1),size(data_snip,2));
        p_all.Left{i} = NaN(size(data_snip,1),size(data_snip,2));
        rho_all.Left{i} = NaN(size(data_snip,1),size(data_snip,2));
        for f = 1:size(data_snip,1)
            for t = 1:size(data_snip,2)
                if strcmpi(stat_test,'corr')
                    [rho,p] = partialcorr(squeeze(data_snip(f,t,:)),stepMod_corr{i}',subjVar_corr);
                else 
                    % Add more stat testing options here
                end

                p_all.Left{i}(f,t) = p;
                rho_all.Left{i}(f,t) = rho;

                if p<0.05
                    h_new =1;
                else
                    h_new = 0.5;
                end
                h_all.Left{i}(f,t) = h_new;

                clear h_new p lme coh_stacked rho
            end

            if mult_corr
                % FDR
                FDR = mafdr(p_all.Left{i}(f,:));
    
                for t=1:length(FDR)
                    if FDR(t)<0.05
                        h_fdr_new =0.5;
                    else
                        h_fdr_new = 1;
                    end
                    h_all.Left{i}(f,t) = h_fdr_new;

                    clear h_fdr_new
                end
    
                clear FDR
            end
        end
        clear data_snip
    end
end

% Right
if unilat==0 || unilat==1
    h_all.Right = cell(1,length(allsteps{1}.Right));
    p_all.Right = cell(1,length(allsteps{1}.Right));
    rho_all.Right = cell(1,length(allsteps{1}.Right));
    subjVar_corr = subjVar{1};
    stepMod_corr = allStepMod{1}.Right;
    for i = 1:length(allsteps{1}.Right)
        data_snip = allsteps{1}.Right{i};
        h_all.Right{i} = NaN(size(data_snip,1),size(data_snip,2));
        p_all.Right{i} = NaN(size(data_snip,1),size(data_snip,2));
        rho_all.Right{i} = NaN(size(data_snip,1),size(data_snip,2));
        for f = 1:size(data_snip,1)
            for t = 1:size(data_snip,2)
                if strcmpi(stat_test,'corr')
                   [rho,p] = partialcorr(squeeze(data_snip(f,t,:)),stepMod_corr{i}',subjVar_corr);
                else 
                     % Add more stat testing options here
                end

                p_all.Right{i}(f,t) = p;
                rho_all.Right{i}(f,t) = rho;

                if p<0.05
                    h_new =1;
                else
                    h_new = 0.5;
                end

                h_all.Right{i}(f,t) = h_new;

                clear h_new p lme coh_stacked rho
            end

            if mult_corr
                % FDR
                FDR = mafdr(p_all.Right{i}(f,:));
    
                for t=1:length(FDR)
                    if FDR(t)<0.05
                        h_fdr_new =0.5;
                    else
                        h_fdr_new = 1;
                    end
                    h_all.Right{i}(f,t) = h_fdr_new;

                    clear h_fdr_new
                end
    
                clear FDR
            end
        end
        clear data_snip data_snip_cond2
    end
end

% Calculate grand averages
if unilat==0 || unilat==-1
    for i = 1:length(data.normCoh.Left)
        grandAvgCoh{1}.Left{i} = mean(allsteps{1}.Left{i},3);
    end
end

if unilat==0 || unilat==1
    for i = 1:length(data.normCoh.Right)
        grandAvgCoh{1}.Right{i} = mean(allsteps{1}.Right{i},3);
    end
end

% Save the significance mask
if saveData
    save(fullfile(data_save_path,'significance_mask.mat'),'h_all','-v7.3');
end

%% Plotting
if unilat==0 || unilat==-1
    for i = 1:length(grandAvgCoh{1}.Left)
        for cond=1:length(grandAvgCoh)
            freq_vec = squeeze(allfc.Left{1}(:,:,1));
    
            figure;
            % h = imagesc(1:100,log2(freq_vec),grandAvgCoh{cond}.Left{i}, 'Interpolation', 'bilinear');
            h = imagesc(1:100,log2(freq_vec),rho_all.Left{i}, 'Interpolation', 'bilinear');
            ax = gca;
            % ticks = logspace(log10(2.5),log10(50),5);
            ticks = logspace(log10(2.5),log10(70),5);
            ax.YTick = log2(ticks);
            ax.YTickLabel = sprintf('%1.2f\n',round(ticks,2));
            xline(mean(allTO{cond}.Left{1}),'LineStyle','--','LineWidth',1.5);
            % ylim([log2(2.5),log2(50)]);
            ylim([log2(2.5),log2(70)]);
            xticks([1,10,20,30,40,50,60,70,80,90,100]);
            xticklabels({'0','10','20','30','40','50','60','70','80','90','100'});
            colormap jet;
            
            % Apply transparency mask to non-significant p-values
            alphaMask = h_all.Left{i};
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
            freq_vec = squeeze(allfc.Right{1}(:,:,1));
    
            figure;
            % h = imagesc(1:100,log2(freq_vec),grandAvgCoh{cond}.Right{i}, 'Interpolation', 'bilinear');
            h = imagesc(1:100,log2(freq_vec),rho_all.Right{i}, 'Interpolation', 'bilinear');
            ax = gca;
            % ticks = logspace(log10(2.5),log10(50),5);
            ticks = logspace(log10(2.5),log10(70),5);
            ax.YTick = log2(ticks);
            ax.YTickLabel = sprintf('%1.2f\n',round(ticks,2));
            xline(mean(allTO{cond}.Right{1}),'LineStyle','--','LineWidth',1.5);
            % ylim([log2(2.5),log2(50)]);
            ylim([log2(2.5),log2(70)]);
            % % xticks([1,10,20,30,40,50,60,70,80,90,100]);
            % % xticklabels({'0','10','20','30','40','50','60','70','80','90','100'});
            xticks([1,20,40,60,80,100]);
            xticklabels({'0','10','20','30','40','50'});
            colormap jet;
    
            % Apply transparency mask to non-significant p-values
            alphaMask = h_all.Right{i};
            set(h, 'AlphaData', alphaMask);
            set(gca,'YDir','normal') 
    
            ylabel('Frequency (Hz)');
            xlabel('% Gait Cycle');
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

