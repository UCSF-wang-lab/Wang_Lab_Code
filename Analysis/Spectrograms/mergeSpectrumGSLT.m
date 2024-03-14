function mergeSpectrumGSLT(varargin)
% When multiple aligned datasets exist within the same experimental session, 
% this function merges the multiple occuring spectral decompositions.
%
% INPUTS:  Optional
%                dataPath[=]    The path to the folder containing the 
%                               spectral decompositions.
%               
%                    keys[=]    Cell array of keys to analyze. Default is
%                               to analyze all keys.
%
%                savePlot[=]    Boolean option to save the resulting merged plots.
%                               Default is false.
%
%                saveData[=]    Boolean option to save the spectral data.
%                               Default is false.
%
%               condition[=]    String describing the condition analyzed, 
%                               e.g. OFFmed_ONdbs. It creates appropriate subfolders
%                               to save the plots in.
%
%                  unilat[=]    0 if bilat (default), -1 if left
%                               unilateral, 1 if right unilateral
%
%
% Example call:
%           keys = {'key0','key1','key2','key3'};
%           mergeSpectrumGSLT('keys',keys,'savePlot',1,'condition','OFFmed_ONdbs_RightAdaptive')
%
% Date:     10/24/2023
% Author:   Eleni Patelaki (eleni.patelaki@ucsf.edu)

%% Parse optional arguments
for i = 1:2:nargin
    switch varargin{i}
        case 'dataPath'
            dataPath = varargin{i+1};
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
    
    plot_save_path = fullfile(parent_dir,'Figures','GSLT_spectrum_figs_adaptive_masked','Merged',condition);
    if ~exist(plot_save_path, 'dir')
       mkdir(plot_save_path);
    end
end

%% Fetch necessary files
% All spectrum datafiles
if ~exist('dataPath','var') || isempty(dataPath)
    [specdatafiles,dataPath] = uigetfile('step*.mat','Select file(s)', 'MultiSelect', 'on');
else
    specdatafiles = dir(fullfile(dataPath, 'step*.mat'));
end

% All frequency & normalization datafiles that correspond to the selected spectrum data files
freqfnames = cellfun(@(x) strrep(x,'step','fc'),specdatafiles,'UniformOutput',false);
normfnames = cellfun(@(x) strrep(x,'step_spec','normalization'),specdatafiles,'UniformOutput',false);
tofnames = cellfun(@(x) strrep(x,'step_spec','toeoff'),specdatafiles,'UniformOutput',false);

%% Ensure that the number of the frequency & normalization datafiles is equal to that of the spectrum datafiles
assert((numel(specdatafiles)==numel(freqfnames))& (numel(freqfnames)==numel(normfnames)) & (numel(freqfnames)==numel(tofnames)));

%% Initialize the aggregated spectrum, normalization structs and gait cycle number structures
allsteps.Left = cell(1,length(keys));
allsteps.Right = cell(1,length(keys));
allfc.Left = cell(1,length(keys));
allfc.Right = cell(1,length(keys));
allnorm.Left = cell(1,length(keys));
allnorm.Right = cell(1,length(keys));
allstepNums.Left = cell(1,length(keys));
allstepNums.Right = cell(1,length(keys));
alltoeOffs.Left = cell(1,length(keys));
alltoeOffs.Right = cell(1,length(keys));

%% Load & merge the spectrum & normalization datafiles
for i = 1:numel(specdatafiles)
    if isstruct(specdatafiles)
        load(fullfile(dataPath,getfield(specdatafiles(i),'name')));
    else
        load(fullfile(dataPath,specdatafiles{i}));
    end
    load(fullfile(dataPath,freqfnames{i}));
    load(fullfile(dataPath,normfnames{i}));
    load(fullfile(dataPath,tofnames{i}));

    if sum(step.Left{1},'all')~=0 %&& sum(step.Right{1},'all')~=0
        stepNums.Left = cellfun(@(x) size(x,3),step.Left,'UniformOutput',false);
        stepNums.Right = cellfun(@(x) size(x,3),step.Right,'UniformOutput',false);
    
        allsteps.Left = cellfun(@(x,y) cat(3,x,y),allsteps.Left,step.Left,'UniformOutput',false);
        allsteps.Right = cellfun(@(x,y) cat(3,x,y),allsteps.Right,step.Right,'UniformOutput',false);
        allfc.Left = cellfun(@(x,y) cat(3,x,y),allfc.Left,fc.Left,'UniformOutput',false);
        allfc.Right = cellfun(@(x,y) cat(3,x,y),allfc.Right,fc.Right,'UniformOutput',false);
        allnorm.Left = cellfun(@(x,y) cat(3,x,y),allnorm.Left,normalization.Left,'UniformOutput',false);
        allnorm.Right = cellfun(@(x,y) cat(3,x,y),allnorm.Right,normalization.Right,'UniformOutput',false);
        allstepNums.Left = cellfun(@(x,y) cat(2,x,y),allstepNums.Left,stepNums.Left,'UniformOutput',false);
        allstepNums.Right = cellfun(@(x,y) cat(2,x,y),allstepNums.Right,stepNums.Right,'UniformOutput',false);
        alltoeOffs.Left = cellfun(@(x,y) cat(2,x,y),alltoeOffs.Left,toeOff.Left,'UniformOutput',false);
        alltoeOffs.Right = cellfun(@(x,y) cat(2,x,y),alltoeOffs.Right,toeOff.Right,'UniformOutput',false);
    
        clear step stepNums normalization toeOff
    end
end


%% Calculate the weighted average of the normalizations matrices & normalize
for i = 1:length(keys)
    for j = 1:numel(specdatafiles)
        if j==1
            sumNormLeft =  allstepNums.Left{i}(j)*allnorm.Left{i}(:,:,j);
            sumNormRight =  allstepNums.Right{i}(j)*allnorm.Right{i}(:,:,j);
        else
            sumNormLeft =  sumNormLeft + allstepNums.Left{i}(j)*allnorm.Left{i}(:,:,j);
            sumNormRight =  sumNormRight + allstepNums.Right{i}(j)*allnorm.Right{i}(:,:,j);
        end
    end
    avgnormLeft = sumNormLeft/sum(allstepNums.Left{i},'all');
    avgnormRight = sumNormRight/sum(allstepNums.Right{i},'all');

    % Apply the normalization (z-scoring)
    if unilat==0 || unilat==-1
        normSpec.Left{i} = (allsteps.Left{i}-avgnormLeft(:,1))./avgnormLeft(:,2);
    end
    if unilat==0 || unilat==1
        normSpec.Right{i} = (allsteps.Right{i}-avgnormRight(:,1))./avgnormRight(:,2);
    end
    clear avgnormLeft avgnormRight sumNormLeft sumNormRight
end

%% Identify the statistically significant effects using FDR correction
% Left
if unilat==0 || unilat==-1
    h_all_left = cell(1,length(normSpec.Left));
    p_all_left = cell(1,length(normSpec.Left));
    for i = 1:length(normSpec.Left)
        data_snip = normSpec.Left{i};
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
    h_all_right = cell(1,length(normSpec.Right));
    p_all_right = cell(1,length(normSpec.Right));
    for i = 1:length(normSpec.Right)
        data_snip = normSpec.Right{i};
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

%% Calculate grand averages
if unilat==0 || unilat==-1
    for i = 1:length(normSpec.Left)
        grandAvgSpec.Left{i} = mean(normSpec.Left{i},3);
        avgToeOff.Left{i} = mean(alltoeOffs.Left{i});
    end
end

if unilat==0 || unilat==1
    for i = 1:length(normSpec.Right)
        grandAvgSpec.Right{i} = mean(normSpec.Right{i},3);
        avgToeOff.Right{i} = mean(alltoeOffs.Right{i});
    end
end

%% Save the grand average data (optional)
if saveData
    save(fullfile(dataPath,'grandAvgSpec.mat'),'grandAvgSpec','-v7.3');
    save(fullfile(dataPath,'normSpec.mat'),'normSpec','-v7.3');
    save(fullfile(dataPath,'allfc.mat'),'allfc','-v7.3');
    save(fullfile(dataPath,'avgToeOff.mat'),'avgToeOff','-v7.3');
end

%% Plotting
% Left
if unilat==0 || unilat==-1
    for i = 1:length(grandAvgSpec.Left)
        freq_vec = squeeze(allfc.Left{1}(:,:,1));

        figure;
        h = imagesc(1:100,log2(freq_vec),grandAvgSpec.Left{i}, 'Interpolation', 'bilinear');
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
        title({'Left';sprintf('%s Grand Average Spectrum',keys{i})});

        % Save plots
        if savePlot    
            saveas(gcf,fullfile(plot_save_path,strcat('Left_',keys{i})),'fig');
            saveas(gcf,fullfile(plot_save_path,strcat('Left_',keys{i})),'tiff');
        end
    end
end

% Right
if unilat==0 || unilat==1
    for i = 1:length(grandAvgSpec.Right)
        freq_vec = squeeze(allfc.Left{1}(:,:,1));

        figure;
        h = imagesc(1:100,log2(freq_vec),grandAvgSpec.Right{i}, 'Interpolation', 'bilinear');
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
        title({'Right';sprintf('%s Grand Average Spectrum',keys{i})});

        % Save plots
        if savePlot
            saveas(gcf,fullfile(plot_save_path,strcat('Right_',keys{i})),'fig');
            saveas(gcf,fullfile(plot_save_path,strcat('Right_',keys{i})),'tiff');
        end
    end
end

end

