function [p_left,p_right] = testCondDiffMultSubjToeOffsGSLT(subjList,cond1,cond2,comparison,varargin)
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
%           Optional
%
%
%                savePlot[=]    Boolean option to save the resulting merged plots.
%                               Default is false.
%
%                saveData[=]    Boolean option to save the coherence data.
%                               Default is false.
%
%                dbs_cond[=]    'dbsopt' for data recorded during the dbs-optimized phase
%                               'pre-prog' for data recorded during the pre-programming phase
%
%   Example call:
%           keys = {'key0','key1','key2','key3'};
%           testCondDiffMultSubjCoherenceGSLT('keys',keys,'saveData',0,'savePlot',0,'comparison','adaptVSmaladapt_right');
%
% Date:     3/24/2024
% Author:   Eleni Patelaki (eleni.patelaki@ucsf.edu)

%% Parse optional arguments
for i = 1:2:nargin-4
    switch varargin{i}
        case 'savePlot'
            savePlot = varargin{i+1};
        case 'saveData'
            saveData = varargin{i+1};
        case 'dbs_cond'
            dbs_cond = varargin{i+1};
    end
end

%% Set to default values if not passed in by user
if ~exist('saveData','var') || isempty(saveData)
    saveData = 0;
end

if ~exist('dbs_cond','var') || isempty(dbs_cond)
    dbs_cond = 'dbsopt';
end

if ~exist('savePlot','var') || isempty(savePlot)
    savePlot = 0;   % Does not save plot by default.
else
    plot_save_path = fullfile('X:\Patient Data\RC+S Data\gait_RCS aggregate data\Eleni\GSLT_coherence\Figures\ToeOffs');
    if ~exist(plot_save_path, 'dir')
       mkdir(plot_save_path);
    end
end

%% Load the data iteratively for each patient
basePath = 'X:\Patient Data\RC+S Data';
allTO{1}.Left = [];
allTO{1}.Right = [];
allTO{2}.Left = [];
allTO{2}.Right = [];
subjVar{1}.Right = [];
subjVar{1}.Left = [];
subjVar{2}.Right = [];
subjVar{2}.Left = [];

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
    data_path1 = fullfile(basePath,sprintf('gait_RCS_%02d',subjList(sub)),dirList(fileIdx_final).name,'Data\Analysis Data\GSLT\Coherence_data_adaptive',cond1);
    toeOff_cond1 = load(fullfile(data_path1, 'alltoeOffs.mat'));
    
    % Load data for condition 2
    data_path2 = fullfile(basePath,sprintf('gait_RCS_%02d',subjList(sub)),dirList(fileIdx_final).name,'Data\Analysis Data\GSLT\Coherence_data_adaptive',cond2);
    toeOff_cond2 = load(fullfile(data_path2, 'alltoeOffs.mat'));

    % Cond 1
    allTO{1}.Left = cat(2,allTO{1}.Left,toeOff_cond1.alltoeOffs.Left{1});
    allTO{1}.Right = cat(2,allTO{1}.Right,toeOff_cond1.alltoeOffs.Right{1});
    subjVar{1}.Right = [subjVar{1}.Right; subjList(sub)*ones(size(toeOff_cond1.alltoeOffs.Right{1},2),1)];
    subjVar{1}.Left = [subjVar{1}.Left; subjList(sub)*ones(size(toeOff_cond1.alltoeOffs.Left{1},2),1)];
    
    % Cond 2
    allTO{2}.Left = cat(2,allTO{2}.Left,toeOff_cond2.alltoeOffs.Left{2});
    allTO{2}.Right = cat(2,allTO{2}.Right,toeOff_cond2.alltoeOffs.Right{2});
    subjVar{2}.Right = [subjVar{2}.Right; subjList(sub)*ones(size(toeOff_cond2.alltoeOffs.Right{1},2),1)];
    subjVar{2}.Left = [subjVar{2}.Left; subjList(sub)*ones(size(toeOff_cond2.alltoeOffs.Left{1},2),1)];
end

%% Test for significant differences between hits and misses (LMEs, condition = fixed effect, patient = random effect)
% Left leg LME
subjVar_stacked = [subjVar{1}.Left; subjVar{2}.Left];
cond_stacked = [ones(length(subjVar{1}.Left),1); 2*ones(length(subjVar{2}.Left),1)];
toeOffs_stacked = [allTO{1}.Left'; allTO{2}.Left'];
   
tbl = table(cond_stacked,subjVar_stacked,toeOffs_stacked);
tbl.cond_stacked = nominal(tbl.cond_stacked);
tbl.subjVar_stacked = nominal(tbl.subjVar_stacked);

lme_left = fitlme(tbl,'toeOffs_stacked ~ 1  + cond_stacked + (1|subjVar_stacked)');
p_left = coefTest(lme_left);

clear lme_left subjVar_stacked cond_stacked toeOffs_stacked tbl

% Right leg LME
subjVar_stacked = [subjVar{1}.Right; subjVar{2}.Right];
cond_stacked = [ones(length(subjVar{1}.Right),1); 2*ones(length(subjVar{2}.Right),1)];
toeOffs_stacked = [allTO{1}.Right'; allTO{2}.Right'];

tbl = table(cond_stacked,subjVar_stacked,toeOffs_stacked);
tbl.cond_stacked = nominal(tbl.cond_stacked);
tbl.subjVar_stacked = nominal(tbl.subjVar_stacked);

lme_right = fitlme(tbl,'toeOffs_stacked ~ 1  + cond_stacked + (1|subjVar_stacked)');
p_right = coefTest(lme_right);

clear lme_right subjVar_stacked cond_stacked toeOffs_stacked tbl

%% Plot the toe-off diffs between hits and misses
figure;

% Scatter plots
% Set up the necessary variables
scatterFactor = 2;
colors = [0.7 0.7 0; 0 0.7 0; 0.7 0 0.7; 0.7 0 0; 0 0 0.7];

aggregTO{1} = allTO{1}.Left';
aggregTO{2} = allTO{2}.Left';
aggregTO{3} = allTO{1}.Right';
aggregTO{4} = allTO{2}.Right';

% aggregSubjVar{1} = subjVar{1}.Left';
% aggregSubjVar{2} = subjVar{2}.Left';
% aggregSubjVar{3} = subjVar{1}.Right';
% aggregSubjVar{4} = subjVar{2}.Right';
% 
% Set up the color matrices
% for box = 1:length(aggregTO)
%     for pnt = 1:length(aggregTO{box})
%         colorMat{box}(pnt,:) = colors(aggregSubjVar{box}(pnt),:);
%     end
% end
% 
% for box = 1:length(aggregTO)
%     hold on;
% 
%     hold on
%     s = scatter((box-1) + ones(size(aggregTO{box})).*(1+(rand(size(aggregTO{box}))-0.5)/scatterFactor),aggregTO{box},30,colorMat{box},'filled'); %[0,0.4,0.4]: teal
%     s.MarkerFaceAlpha = 0.5;
% 
%     clear s
% end
% hold off;

x = [];
for box = 1:length(aggregTO)
    hold on
    x = [x; (box-1) + ones(size(aggregTO{box})).*(1+(rand(size(aggregTO{box}))-0.5)/scatterFactor)];
end

y = [allTO{1}.Left'; allTO{2}.Left'; allTO{1}.Right'; allTO{2}.Right'];
% divide the toe-off occurence values by 2 to convert from % step cycle to
% % gait cycle
y = y/2;

g_subj = [subjVar{1}.Left; subjVar{2}.Left; subjVar{1}.Right; subjVar{2}.Right];
g_side = [ones(length(allTO{1}.Left),1); 1*ones(length(allTO{2}.Left),1); 2*ones(length(allTO{1}.Right),1); 2*ones(length(allTO{2}.Right),1)];
% g = {g_subj,g_side};
g = {g_subj};
hs = gscatter(x,y,g,colors,'o',5,'on',[],[]);

% Create the legend
hLegend = legend('location','northeastoutside');
for legend_item = 1:length(hLegend.String)
    hLegend.String{legend_item}= ['gRCS0',hLegend.String{legend_item}];
end

for line_iterm = 1:length(hs)
    set(hs(line_iterm), 'LineWidth',1);
end

% Title and ticks
xlim([0.5 2.5]);
ylim([0 50]);
xticks([1 2]);
xticklabels({'Hits','Misses'});
title('Toe-off Occurence');
ylabel('% Gait Cycle');
set(gca,'FontSize',15);

% Save plot
if savePlot
    saveas(gcf,fullfile(plot_save_path,[comparison,'_',dbs_cond,'.fig']));
    saveas(gcf,fullfile(plot_save_path,[comparison,'_',dbs_cond,'.tiff']));
end

end

