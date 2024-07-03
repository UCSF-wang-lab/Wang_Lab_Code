function convertToMergedCoherenceStruct(unilat)

keys = {'key0','key1','key2','key3'};
cohPairs = nchoosek(keys,2);

% Get the file that needs to be converted.
% [cohdatafname,dataPath] = uigetfile('step_coh*.mat','Select file');
% load(fullfile(dataPath,cohdatafname));

% All coherence datafiles
if ~exist('dataPath','var') || isempty(dataPath)
    [cohdatafnames,dataPath] = uigetfile('step_coh*.mat','Select file(s)', 'MultiSelect', 'on');
else
    cohdatafnames = dir(fullfile(dataPath, 'step_coh*.mat'));
end
if ~iscell(cohdatafnames)==1
    cohdatafnames = {cohdatafnames};
end
tofnames = cellfun(@(x) strrep(x,'step_coh','toeoff'),cohdatafnames,'UniformOutput',false);


load(fullfile(dataPath,'avgnorm.mat'));

allsteps.Left = cell(1,size(cohPairs,1));
allsteps.Right = cell(1,size(cohPairs,1));
alltoeOffs.Left = cell(1,size(cohPairs,1));
alltoeOffs.Right = cell(1,size(cohPairs,1));

%% Load & merge the coherence & normalization datafiles
for i = 1:numel(cohdatafnames)
    if isstruct(cohdatafnames)
        load(fullfile(dataPath,getfield(cohdatafnames(i),'name')));
    else
        load(fullfile(dataPath,cohdatafnames{i}));
    end
    load(fullfile(dataPath,tofnames{i}));

    if sum(step.Left{1},'all')~=0 % && sum(step.Right{1},'all')~=0
        allsteps.Left = cellfun(@(x,y) cat(3,x,y),allsteps.Left,step.Left,'UniformOutput',false);
        allsteps.Right = cellfun(@(x,y) cat(3,x,y),allsteps.Right,step.Right,'UniformOutput',false);
        alltoeOffs.Left = cellfun(@(x,y) cat(2,x,y),alltoeOffs.Left,toeOff.Left,'UniformOutput',false);
        alltoeOffs.Right = cellfun(@(x,y) cat(2,x,y),alltoeOffs.Right,toeOff.Right,'UniformOutput',false);

        clear step stepNums normalization toeOff
    end
end

for i = 1:size(cohPairs,1)
    % Apply the normalization (z-scoring)
    if unilat==0 || unilat==-1
        normCoh.Left{i} = (allsteps.Left{i}-avgnorm.Left{i}(:,1))./avgnorm.Left{i}(:,2);
    end
    if unilat==0 || unilat==1
        normCoh.Right{i} = (allsteps.Right{i}-avgnorm.Right{i}(:,1))./avgnorm.Right{i}(:,2);
    end
end

save(fullfile(dataPath,'normCoh.mat'),'normCoh','-v7.3');


freqfnames = cellfun(@(x) strrep(x,'step','fc'),cohdatafnames,'UniformOutput',false);
% assert(length(freqfnames)==1);
load(fullfile(dataPath,freqfnames{1}));
allfc = fc;
save(fullfile(dataPath,'allfc.mat'),'allfc','-v7.3');

% Calculate the toe-off grand averages
if unilat==0 || unilat==-1
    for i = 1:length(alltoeOffs.Left)
        grandAvgCoh.Left{i} = mean(normCoh.Left{i},3);
        avgToeOff.Left{i} = mean(alltoeOffs.Left{i});
    end
end

if unilat==0 || unilat==1
    for i = 1:length(normCoh.Right)
        grandAvgCoh.Right{i} = mean(normCoh.Right{i},3);
        avgToeOff.Right{i} = mean(alltoeOffs.Right{i});
    end
end

save(fullfile(dataPath,'avgToeOff.mat'),'avgToeOff','-v7.3');
save(fullfile(dataPath,'alltoeOffs.mat'),'alltoeOffs','-v7.3');

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
        
        % % Apply transparency mask to non-significant p-values
        % alphaMask = h_all_left{i};
        % set(h, 'AlphaData', alphaMask);
        set(gca,'YDir','normal') 

        ylabel('Frequency (Hz)');
        title({'Left';sprintf('%s to %s Grand Average Coherence',cohPairs{i,1},cohPairs{i,2})});

        % % Save plots
        % if savePlot    
        %     saveas(gcf,fullfile(plot_save_path,strcat('Left_',cohPairs{i,1},'_',cohPairs{i,2})),'fig');
        %     saveas(gcf,fullfile(plot_save_path,strcat('Left_',cohPairs{i,1},'_',cohPairs{i,2})),'tiff');
        % end
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

        % % Apply transparency mask to non-significant p-values
        % alphaMask = h_all_right{i};
        % set(h, 'AlphaData', alphaMask); 
        set(gca,'YDir','normal') 

        ylabel('Frequency (Hz)');
        title({'Right';sprintf('%s to %s Grand Average Coherence',cohPairs{i,1},cohPairs{i,2})});

        % % Save plots
        % if savePlot
        %     saveas(gcf,fullfile(plot_save_path,strcat('Right_',cohPairs{i,1},'_',cohPairs{i,2})),'fig');
        %     saveas(gcf,fullfile(plot_save_path,strcat('Right_',cohPairs{i,1},'_',cohPairs{i,2})),'tiff');
        % end
    end
end

end