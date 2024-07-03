function plotAggregateHitRates(subjList,condStr,dbs_cond,varargin)
% Calculates and plots GSLT performance across blocks
%
% INPUTS:    
%       Required:
%                           subjList[=]     Vector containing the patients numbers to
%                                           be analyzed, e.g. [2,3,5] for gRCS02,
%                                           grRCS03 and gRCS05.
%
%                            condStr[=]     String indicating the condition to be analyzed. 
%                                           e.g. ONmed_ONdbs 
%
%                           dbs_cond[=]     'dbsopt' for data recorded during the dbs-optimized phase
%                                           'pre-prog' for data recorded during the pre-programming phase
%
%       Optional:
%                          
%
% OUTPUTS:   NONE
%
% Example call: plotAggregateHitRates([2,3,4],'ONmed_ONdbs','dbsopt');
%
% Author:   Eleni Patelaki
% Date:     3/27/2024

%% Parse optional arguments
for i = 1:2:nargin-3
    switch varargin{i}
        case 'savePlot'
            savePlot = varargin{i+1};
    end
end

%% Set to default values if not passed in by user
if ~exist('savePlot','var') || isempty(savePlot)
    savePlot = 0;   % Does not save plot by default.
else
    plot_save_path = fullfile('X:\Patient Data\RC+S Data\gait_RCS aggregate data\Eleni\GSLT_Behavioral_performance\Figures');
    if ~exist(plot_save_path, 'dir')
       mkdir(plot_save_path);
    end
end

%% Load hit rate data from all the patients
basePath = 'X:\Patient Data\RC+S Data';
for sub = 1:length(subjList)
    dirList = dir(fullfile(basePath,sprintf('gait_RCS_%02d',subjList(sub))));
    fileIdx = find(contains(lower({dirList.name}),dbs_cond));

    if isempty(fileIdx)
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
    data_path = fullfile(basePath,sprintf('gait_RCS_%02d',subjList(sub)),dirList(fileIdx_final).name,'Data\Analysis Data\GSLT\Behavioral_performance');
    currHR = load(fullfile(data_path, strcat('hitrate_all_',condStr)));
    aggregHR(sub,:) = currHR.hitrate_all;
end

%% Plot the hit rates from all the patients
colors = [0.7 0.7 0; 0 0.7 0; 0.7 0 0.7; 0.7 0 0; 0 0 0.7];

figure; 
colororder(colors);
plot((1:7)',aggregHR','LineWidth',1.5,'Marker','.','MarkerSize',20);
legend(stArray = arrayfun(@(x) strcat('gRCS0',num2str(x)), subjList, 'UniformOutput', false));

% Title and ticks
xlim([0.5, 7.5]);
xticks(1:7);
xticklabels({'R','R','S','S','S','R','S'});
xlabel('Block Order');
ylim([0 1]);
ylabel('% Hits');
title('Step Adaptation Performance');
set(gca,'FontSize',15);

% Save plots
if savePlot
    saveas(gcf,fullfile(plot_save_path,[condStr,'_',dbs_cond,'.fig']));
    saveas(gcf,fullfile(plot_save_path,[condStr,'_',dbs_cond,'.tiff']));
end


end