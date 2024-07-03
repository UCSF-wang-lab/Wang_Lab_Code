function testCondDiffMultSubjSpectrumGSLT_RvS_PostHoc(subjList,cond1,cond2,comparison,cluster_frange,cluster_trange, varargin)
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
%                   cond2[=]    Condition 2 string. Must be the same as the
%                               corresponding data folder name, 
%                               e.g. ONmed_ONdbs_RightMaladaptive
%
%              comparison[=]    String describing the comparison tested, 
%                               e.g. ONmed_adaptVSmaladapt_right. It creates appropriate subfolders
%                               to save the plots in.
%           
%          cluster_frange[=]    Frequency range of the analyzed cluster
%                               e.g. [5,15] (in Hz)
%
%          cluster_trange[=]    Time range of the analyzed cluster
%                               e.g. [0,20] in (% gait cycle)
%
%           Optional
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
%                  unilat[=]    0 if bilat (default), -1 if left
%                               unilateral, 1 if right unilateral. Default
%                               is 0.
%
%               stat_test[=]    'lme' for linear mixed-effects
%                               models,'ttest'for ttests
%
%                dbs_cond[=]    'dbsopt' for data recorded during the dbs-optimized phase
%                               'pre-prog' for data recorded during the pre-programming phase
%
%   Example call:
%           keys = {'key0','key1','key2','key3'};
%           testCondDiffMultSubjCoherenceGSLT('keys',keys,'saveData',0,'savePlot',0,'comparison','adaptVSmaladapt_right','unilat',0);
%
% Date:     7/1/2024
% Author:   Eleni Patelaki (eleni.patelaki@ucsf.edu)

%% Parse optional arguments
for i = 1:2:nargin-6
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
        case 'stat_test'
            stat_test = varargin{i+1};
        case 'dbs_cond'
            dbs_cond = varargin{i+1};
    end
end

%% Define the gait cycle phase
if cluster_trange(2)<50
    gcphase = 'preSwing';
else
    gcphase = 'midSwing';
end

%% Set to default values if not passed in by user
if ~exist('keys','var') || isempty(keys)
    keys = {'key0','key1','key2','key3'};
end

if ~exist('saveData','var') || isempty(saveData)
    saveData = 0;
else
    data_save_path = fullfile('X:\Patient Data\RC+S Data\gait_RCS aggregate data\Eleni\GSLT_coherence\Data',[comparison,'_',dbs_cond,'_',stat_test,'_',gcphase,'_subs_2_3_4_HigherFreq_upto70_midSwingClust_v1']);
    if ~exist(data_save_path, 'dir')
       mkdir(data_save_path);
    end
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
    plot_save_path = fullfile('X:\Patient Data\RC+S Data\gait_RCS aggregate data\Eleni\GSLT_spectrum\Figures',[comparison,'_',dbs_cond,'_',stat_test,'_',gcphase,'_subs_2_3_4_HigherFreq_upto70_v1']);
    if ~exist(plot_save_path, 'dir')
       mkdir(plot_save_path);
    end
end

%% Load the data iteratively for each patient
basePath = 'X:\Patient Data\RC+S Data';
allsteps{1}.Left = cell(1,size(keys,1));
allsteps{1}.Right = cell(1,size(keys,1));
allsteps{2}.Left = cell(1,size(keys,1));
allsteps{2}.Right = cell(1,size(keys,1));
allTO{1}.Left = cell(1,size(keys,1));
allTO{1}.Right = cell(1,size(keys,1));
allTO{2}.Left = cell(1,size(keys,1));
allTO{2}.Right = cell(1,size(keys,1));
allfc.Left = cell(1,size(keys,1));
allfc.Right = cell(1,size(keys,1));
subjVar{1} = [];
subjVar{2} = [];

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
    data_path1 = fullfile(basePath,sprintf('gait_RCS_%02d',subjList(sub)),dirList(fileIdx_final).name,'Data\Analysis Data\GSLT\Spectrum_data_adaptive',cond1);
    data_cond1 = load(fullfile(data_path1, 'normSpec.mat'));
    fc_cond1 = load(fullfile(data_path1, 'allfc.mat'));
    toeOff_cond1 = load(fullfile(data_path1, 'alltoeOffs.mat'));
    
    % Load data for condition 2
    data_path2 = fullfile(basePath,sprintf('gait_RCS_%02d',subjList(sub)),dirList(fileIdx_final).name,'Data\Analysis Data\GSLT\Spectrum_data_adaptive',cond2);
    data_cond2 = load(fullfile(data_path2, 'normSpec.mat'));
    fc_cond2 = load(fullfile(data_path2, 'allfc.mat'));
    toeOff_cond2 = load(fullfile(data_path2, 'alltoeOffs.mat'));

    % Left
    if unilat==0 || unilat==-1
        % Cond 1
        allsteps{1}.Left = cellfun(@(x,y) cat(3,x,y),allsteps{1}.Left,data_cond1.normSpec.Left,'UniformOutput',false);
        allTO{1}.Left = cellfun(@(x,y) cat(2,x,y),allTO{1}.Left,toeOff_cond1.alltoeOffs.Left,'UniformOutput',false);
        subjVar{1} = [subjVar{1}; sub*ones(size(data_cond1.normSpec.Left{1},3),1)];

        % Cond 2
        allsteps{2}.Left = cellfun(@(x,y) cat(3,x,y),allsteps{2}.Left,data_cond2.normSpec.Left,'UniformOutput',false);
        allTO{2}.Left = cellfun(@(x,y) cat(2,x,y),allTO{2}.Left,toeOff_cond2.alltoeOffs.Left,'UniformOutput',false);
        subjVar{2} = [subjVar{2}; sub*ones(size(data_cond2.normSpec.Left{1},3),1)];

        allfc.Left = fc_cond1.allfc.Left;
    end

    % Right
    if unilat==0 || unilat==1
        % Cond 1
        allsteps{1}.Right = cellfun(@(x,y) cat(3,x,y),allsteps{1}.Right,data_cond1.normSpec.Right,'UniformOutput',false);
        allTO{1}.Right = cellfun(@(x,y) cat(2,x,y),allTO{1}.Right,toeOff_cond1.alltoeOffs.Right,'UniformOutput',false);
        if unilat==1
            subjVar{1} = [subjVar{1}; sub*ones(size(data_cond1.normSpec.Right{1},3),1)];
        end

        % Cond 2
        allsteps{2}.Right = cellfun(@(x,y) cat(3,x,y),allsteps{2}.Right,data_cond2.normSpec.Right,'UniformOutput',false);
        allTO{2}.Right = cellfun(@(x,y) cat(2,x,y),allTO{2}.Right,toeOff_cond2.alltoeOffs.Right,'UniformOutput',false);
        if unilat==1
            subjVar{2} = [subjVar{2}; sub*ones(size(data_cond2.normSpec.Right{1},3),1)];
        end

        allfc.Right = fc_cond1.allfc.Right;
    end

end

%% Load the significance mask
[mask_fname,mask_path] = uigetfile('*.mat');
mask = load(fullfile(mask_path,mask_fname));

%% Identify the statistically significant effects using permutation tests (FDR correction)
% Left
if unilat==0 || unilat==-1
    subjVar_stacked = [subjVar{1}; subjVar{2}];
    cond_stacked = [ones(length(subjVar{1}),1); 2*ones(length(subjVar{2}),1)];
    for i = 1:length(allsteps{1}.Left)
        avgSpec_cond1 = [];
        avgSpec_cond2 = [];
        data_snip_cond1 = allsteps{1}.Left{i};
        data_snip_cond2 = allsteps{2}.Left{i};
        assert(size(data_snip_cond1,1)==size(data_snip_cond2,1));
        assert(size(data_snip_cond1,2)==size(data_snip_cond2,2));

        % Calculate the average coherence within the cluster
        % [ind_freq,ind_time] = find(mask.h_all.Left{i}==1);
        % ind_freq_1 = ind_freq(allfc.Left{1}(ind_freq)>= cluster_frange(1) & allfc.Left{1}(ind_freq)<= cluster_frange(2) & ind_time>=cluster_trange(1) & ind_time<=cluster_trange(2));
        % ind_time_1 = ind_time(allfc.Left{1}(ind_freq)>= cluster_frange(1) & allfc.Left{1}(ind_freq)<= cluster_frange(2) & ind_time>=cluster_trange(1) & ind_time<=cluster_trange(2));
        % sumCoh_cond1 = squeeze(sum(data_snip_cond1(ind_freq_1,ind_time_1,:),[1,2]));
        % sumCoh_cond2 = squeeze(sum(data_snip_cond2(ind_freq_1,ind_time_1,:),[1,2]));
        % 
        % clear ind_freq ind_time

        for itrial = 1:size(data_snip_cond1,3)
            sumSpec_cond1 = 0;
            count_cond1 = 0;
            for ifreq=1:size(data_snip_cond1,1)
                for itime=1:size(data_snip_cond1,2)
                    if mask.h_all.Left{i}(ifreq,itime)==1 && allfc.Left{1}(ifreq)>= cluster_frange(1) && allfc.Left{1}(ifreq)<= cluster_frange(2) && itime>=cluster_trange(1) && itime<=cluster_trange(2)
                        sumSpec_cond1 = sumSpec_cond1 + data_snip_cond1(ifreq,itime,itrial);
                        count_cond1 = count_cond1 + 1;
                    end
                end
            end
            avgSpec_cond1 = [avgSpec_cond1,sumSpec_cond1/count_cond1];
        end


        for itrial = 1:size(data_snip_cond2,3)
            sumSpec_cond2 = 0;
            count_cond2 = 0;
            for ifreq=1:size(data_snip_cond2,1)
                for itime=1:size(data_snip_cond2,2)
                    if mask.h_all.Left{i}(ifreq,itime)==1 && allfc.Left{1}(ifreq)>= cluster_frange(1) && allfc.Left{1}(ifreq)<= cluster_frange(2) && itime>=cluster_trange(1) && itime<=cluster_trange(2)
                        sumSpec_cond2 = sumSpec_cond2 + data_snip_cond2(ifreq,itime,itrial);
                        count_cond2 = count_cond2 + 1;
                    end
                end
            end
            avgSpec_cond2 = [avgSpec_cond2,sumSpec_cond2/count_cond2];
        end

        if strcmpi(stat_test,'ttest')
            % [p,~] = ranksum(sumCoh_cond1,sumCoh_cond2);
            [~,p] = ttest2(avgSpec_cond1,avgSpec_cond2);
        elseif strcmpi(stat_test,'lme')
            spec_stacked = [avgSpec_cond1';avgSpec_cond2'];
           
            if all(isnan(spec_stacked))
                p = NaN;
            else
                % LME
                tbl = table(cond_stacked,subjVar_stacked,spec_stacked);
                tbl.cond_stacked = nominal(tbl.cond_stacked);
                tbl.subjVar_stacked = nominal(tbl.subjVar_stacked);
    
                lme = fitlme(tbl,'spec_stacked ~ 1  + cond_stacked + (1|subjVar_stacked)');
                p = coefTest(lme);
            end

        end

        % Display the p-values
        fprintf('Left %s to %s spectrum cluster p = %f\n',cohPairs{i,1},cohPairs{i,2},p);

        clear data_snip_cond1 data_snip_cond2 ind_freq ind_time sumSpec_cond1 sumSpec_cond2 spec_stacked lme p
    end
end

% Right
if unilat==0 || unilat==1
    subjVar_stacked = [subjVar{1}; subjVar{2}];
    cond_stacked = [ones(length(subjVar{1}),1); 2*ones(length(subjVar{2}),1)];
    for i = 1:length(allsteps{1}.Right)
        avgSpec_cond1 = [];
        avgSpec_cond2 = [];
        data_snip_cond1 = allsteps{1}.Right{i};
        data_snip_cond2 = allsteps{2}.Right{i};
        assert(size(data_snip_cond1,1)==size(data_snip_cond2,1));
        assert(size(data_snip_cond1,2)==size(data_snip_cond2,2));

        % Calculate the average coherence within the cluster
        % [ind_freq,ind_time] = find(mask.h_all.Right{i}==1);
        % ind_freq_1 = ind_freq(allfc.Right{1}(ind_freq)>= cluster_frange(1) & allfc.Right{1}(ind_freq)<= cluster_frange(2) & ind_time>=cluster_trange(1) & ind_time<=cluster_trange(2));
        % ind_time_1 = ind_time(allfc.Right{1}(ind_freq)>= cluster_frange(1) & allfc.Right{1}(ind_freq)<= cluster_frange(2) & ind_time>=cluster_trange(1) & ind_time<=cluster_trange(2));
        % sumCoh_cond1 = squeeze(sum(data_snip_cond1(ind_freq_1,ind_time_1,:),[1,2]));
        % sumCoh_cond2 = squeeze(sum(data_snip_cond2(ind_freq_1,ind_time_1,:),[1,2]));
        % clear ind_freq ind_time

        for itrial = 1:size(data_snip_cond1,3)
            sumSpec_cond1 = 0;
            count_cond1 = 0;
            for ifreq=1:size(data_snip_cond1,1)
                for itime=1:size(data_snip_cond1,2)
                    if mask.h_all.Right{i}(ifreq,itime)==1 && allfc.Right{1}(ifreq)>= cluster_frange(1) && allfc.Right{1}(ifreq)<= cluster_frange(2) && itime>=cluster_trange(1) && itime<=cluster_trange(2)
                        sumSpec_cond1 = sumSpec_cond1 + data_snip_cond1(ifreq,itime,itrial);
                        count_cond1 = count_cond1 + 1;
                    end
                end
            end
            avgSpec_cond1 = [avgSpec_cond1,sumSpec_cond1/count_cond1];
        end


        for itrial = 1:size(data_snip_cond2,3)
            sumSpec_cond2 = 0;
            count_cond2 = 0;
            for ifreq=1:size(data_snip_cond2,1)
                for itime=1:size(data_snip_cond2,2)
                    if mask.h_all.Right{i}(ifreq,itime)==1 && allfc.Right{1}(ifreq)>= cluster_frange(1) && allfc.Right{1}(ifreq)<= cluster_frange(2) && itime>=cluster_trange(1) && itime<=cluster_trange(2)
                        sumSpec_cond2 = sumSpec_cond2 + data_snip_cond2(ifreq,itime,itrial);
                        count_cond2 = count_cond2 + 1;
                    end
                end
            end
            avgSpec_cond2 = [avgSpec_cond2,sumSpec_cond2/count_cond2];
        end

        if strcmpi(stat_test,'ttest')
            % [p,~] = ranksum(sumCoh_cond1,sumCoh_cond2);
            [~,p] = ttest2(avgSpec_cond1,avgSpec_cond2);
        elseif strcmpi(stat_test,'lme')
            spec_stacked = [avgSpec_cond1';avgSpec_cond2'];
           
            if all(isnan(spec_stacked))
                p = NaN;
            else
                % LME
                tbl = table(cond_stacked,subjVar_stacked,spec_stacked);
                tbl.cond_stacked = nominal(tbl.cond_stacked);
                tbl.subjVar_stacked = nominal(tbl.subjVar_stacked);
    
                lme = fitlme(tbl,'spec_stacked ~ 1  + cond_stacked + (1|subjVar_stacked)');
                p = coefTest(lme);
            end
        end

        % Display the p-values
        fprintf('Right %s to %s coherence cluster p = %f\n',cohPairs{i,1},cohPairs{i,2},p);

        clear data_snip_cond1 data_snip_cond2 ind_freq ind_time sumSpec_cond1 sumSpec_cond2 spec_stacked lme p
    end
end

% Calculate grand averages
if unilat==0 || unilat==-1
    for i = 1:length(data_cond1.normSpec.Left)
        grandAvgSpec{1}.Left{i} = mean(allsteps{1}.Left{i},3);
        grandAvgSpec{2}.Left{i} = mean(allsteps{2}.Left{i},3);
    end
end

if unilat==0 || unilat==1
    for i = 1:length(data_cond1.normSpec.Right)
        grandAvgSpec{1}.Right{i} = mean(allsteps{1}.Right{i},3);
        grandAvgSpec{2}.Right{i} = mean(allsteps{2}.Right{i},3);
    end
end

%% Plotting
if unilat==0 || unilat==-1
    for i = 1:length(grandAvgSpec{1}.Left)
        for cond=1:length(grandAvgSpec)
            freq_vec = squeeze(allfc.Left{1}(:,:,1));

            figure;
            h = imagesc(1:100,log2(freq_vec),grandAvgSpec{cond}.Left{i}, 'Interpolation', 'bilinear');
            ax = gca;
            % ticks = logspace(log10(2.5),log10(50),5);
            ticks = logspace(log10(2.5),log10(70),5);
            ax.YTick = log2(ticks);
            ax.YTickLabel = sprintf('%1.2f\n',round(ticks,2));
            xline(mean(allTO{cond}.Left{1}),'LineStyle','--','LineWidth',1.5);
            % ylim([log2(2.5),log2(50)]);
            ylim([log2(2.5),log2(70)]);
            % xticks([1,10,20,30,40,50,60,70,80,90,100]);
            % xticklabels({'0','10','20','30','40','50','60','70','80','90','100'});
            xticks([1,20,40,60,80,100]);
            xticklabels({'0','10','20','30','40','50'});
            colormap jet;

            % Apply transparency mask to non-significant p-values
            alphaMask = ones(size(mask.h_all.Left{1}))*0.5;
            for ifreq=1:size(mask.h_all.Left{i},1)
                for itime=1:size(mask.h_all.Left{1},2)
                    if mask.h_all.Left{i}(ifreq,itime)==1 && allfc.Left{1}(ifreq)>= cluster_frange(1) && allfc.Left{1}(ifreq)<= cluster_frange(2) && itime>=cluster_trange(1) && itime<=cluster_trange(2)
                        alphaMask(ifreq,itime) = 1;
                    end
                end
            end

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
    for i = 1:length(grandAvgSpec{1}.Right)
        for cond=1:length(grandAvgSpec)
            freq_vec = squeeze(allfc.Right{1}(:,:,1));

            figure;
            h = imagesc(1:100,log2(freq_vec),grandAvgSpec{cond}.Right{i}, 'Interpolation', 'bilinear');
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
            % alphaMask = mask.h_all.Right{i};
            alphaMask = ones(size(mask.h_all.Right{1}))*0.5;
            for ifreq=1:size(mask.h_all.Right{i},1)
                for itime=1:size(mask.h_all.Right{1},2)
                    if mask.h_all.Right{i}(ifreq,itime)==1 && allfc.Right{1}(ifreq)>= cluster_frange(1) && allfc.Right{1}(ifreq)<= cluster_frange(2) && itime>=cluster_trange(1) && itime<=cluster_trange(2)
                        alphaMask(ifreq,itime) = 1;
                    end
                end
            end

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

