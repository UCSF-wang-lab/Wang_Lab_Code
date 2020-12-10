function [successPerc,nSuccess,nFail,targetTable] = calcTaskPerformance(filename,plotFlag)
%% calcTaskPerformance
% Calculates the patient's performance on the Gait Sequence Learning Task.
% Assumes the number of blocks, block order, and number of targets per
% block. Outputs the success percent by block, the number of successes and
% and fails of the left and right targets, and the original table that is
% grabs from the file. Run this file after extracting only the target data
% from the larger csv file using "processTargetStepInfo."
%
% INPUTS:   filename    [=] The name of the file containing ONLY target
%                           data. 
%           plotFlag    [=] Optional integer input to plot the
%           performance. 0 or empty = no plots, 1 = plot.
%
% OUTPUTS:  successPerc [=] The success percentage across all blocks. A
%                           3xnBlock matrix. First row is the combined left
%                           and right target success, second row is the 
%                           number of successful left target hits, and the 
%                           third is the number of successful right target 
%                           hits.
%           nSuccess    [=] The number of left and right target hits.
%                           2xnBlock matrix. First row is left, second row
%                           is right.
%           nFail       [=] Left and right target misses. 2xnBlock matrix.
%                           First row is left and second row is right.
%           targetTable [=] The raw data in table format.
%
% Author:   Kenneth Louie
% Date:     10/21/2020

% Constants
nBlocks = 7;
blockCond = [0,0,1,1,1,0,1];    % 0 = rand, 1 = sequence
blockNames = {'R1','R2','S1','S2','S3','R3','S4'};
targetsPerBlock = 120;

% Read in data
targetTable = readtable(filename);

% Preallocate output
successPerc = zeros(3,nBlocks);
nSuccess = zeros(2,nBlocks);
nFail = zeros(2,nBlocks);

% Go through each block calculating success percent, number of successes
% for left and right foot, and number of fails for left and right foot
for i = 1:nBlocks
    % Find targets in block i
    blockInds = find(targetTable.TargetNumber >= (i-1)*targetsPerBlock ...
        & targetTable.TargetNumber < i*targetsPerBlock);
    
    % Extract block data
    blockData = targetTable(blockInds,:);
    
    % Extract left and right hits
    leftHitInds = strcmp(blockData.Side,'Left');
    leftHits = sum(leftHitInds);
    rightHits = sum(~leftHitInds);
    
    % Success
    successPerc(1,i) = (leftHits+rightHits)/targetsPerBlock;
    successPerc(2,i) = leftHits/(targetsPerBlock/2);
    successPerc(3,i) = rightHits/(targetsPerBlock/2);
    
    nSuccess(1,i) = leftHits;
    nSuccess(2,i) = rightHits;
    
    % Fails
    nFail(1,i) = (targetsPerBlock/2)-leftHits;
    nFail(2,i) = (targetsPerBlock/2)-rightHits;
end

if exist('plotFlag','var') && plotFlag
    % Colors
    colormap = CBMap('Paired',9);
    colors = [colormap(2,:); ...
        colormap(1,:);...
        colormap(9,:)];
    
    figure(1); clf;
    barHand = bar(successPerc'.*100,1,'facecolor','flat');
    for k = 1:size(successPerc',2)
        barHand(k).FaceColor = colors(k,:);
    end
    xlabel('Block #'); ylabel('Success (%)');
    xticklabels(blockNames);
    legend('Combined','Left','Right','Box','off');
    
    figure(2); clf;
    bar(successPerc(1,:).*100,0.8,'FaceColor',colors(1,:));
    xlabel('Block #'); ylabel('Success (%)'); title('Combined');
    xticklabels(blockNames);
    
    figure(3); clf;
    bar(successPerc(2,:).*100,0.8,'FaceColor',colors(2,:));
    xlabel('Block #'); ylabel('Success (%)'); title('Left Targets');
    xticklabels(blockNames);
    
    figure(4); clf;
    bar(successPerc(3,:).*100,0.8,'FaceColor',colors(3,:));
    xlabel('Block #'); ylabel('Success (%)'); title('Right Targets');
    xticklabels(blockNames);
    
    figure_format(6,6,11);
end

end