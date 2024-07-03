function gaitEventTable = addGSLTinfo(gaitEventTable,start_trim_time)
% Adds the GSLT-related information to the gait event table, i.e. which
% steps are hits (values: 0.8, 1, 1.2) and which are misses (values: 0)
%
% INPUTS: Required 
%           gaitEventTable[=]  table containing the gait events as columns: RHS, RTO, LHS, LTO      
% 
%          start_trim_time[=]  time mark before which gait events are
%                              discarded since they probably happened
%                              before the task started (before the first
%                              successfully adaptive step)
%
% OUTPUTS:    
%           gaitEventTable[=]  updated version of the gait event table,
%                              containing 2 new columns: AdaptiveLeft,
%                              AdaptiveRight containing information about
%                              whether the steps of each row were
%                              hits/misses.
%
% Author:   Eleni Patelaki
% Contributor: Eleni Patelaki
% Date:     3/11/24

% Read the processed and Xsens-aligned Cirris file
[filename,filepath]=uigetfile('*processed_targets_aligned.csv');
cirris_table=readtable(fullfile(filepath,filename));

cirris_max_time = max(cirris_table.TaskTimer);

% Set to NaN all the gaitEventTable events that happened before a given
% timepoint (optional, for cases where there is noise at the beginning of
% the Cirris stream) and after the end of the last Cirris event
%RTO
rto_tmp = gaitEventTable.RTO;
rto_tmp((rto_tmp<start_trim_time)|(rto_tmp>cirris_max_time))=NaN;
gaitEventTable.RTO = rto_tmp;

%RHS
rhs_tmp = gaitEventTable.RHS;
rhs_tmp((rhs_tmp<start_trim_time)|(rhs_tmp>cirris_max_time))=NaN;
gaitEventTable.RHS = rhs_tmp;

%RTO
lto_tmp = gaitEventTable.LTO;
lto_tmp((lto_tmp<start_trim_time)|(lto_tmp>cirris_max_time))=NaN;
gaitEventTable.LTO = lto_tmp;

%RHS
lhs_tmp = gaitEventTable.LHS;
lhs_tmp((lhs_tmp<start_trim_time)|(lhs_tmp>cirris_max_time))=NaN;
gaitEventTable.LHS = lhs_tmp;


% Get rid of rows with all nan
removeInd = [];
for i = 1:height(gaitEventTable)
    if sum(isnan(gaitEventTable{i,:})) == 4
        removeInd = [removeInd;i];
    end
end
gaitEventTable(removeInd,:) = [];

% Extract the vectors of each gait event separately
rto = gaitEventTable.RTO;
rhs = gaitEventTable.RHS;
lto = gaitEventTable.LTO;
lhs = gaitEventTable.LHS;

% Add two news columns (one for right, one for left) initilaized with 0
% (i.e maladaptive). All steps are considered maladaptive unless they can
% be found in the cirris file. If adaptive, the columns contain the step
% modification ratio
gaitEventTable.AdaptiveRight = zeros(height(gaitEventTable),1);
gaitEventTable.AdaptiveLeft = zeros(height(gaitEventTable),1);
gaitEventTable.BlockNumRight = zeros(height(gaitEventTable),1);
gaitEventTable.BlockNumLeft = zeros(height(gaitEventTable),1);

count1 = 0;
count11 = 0;
count12 = 0;
count2 = 0;
count3 = 0;
count4 = 0;
timediff_1 = [];
timediff_2 = [];

for i = 1:height(cirris_table)
    curr_line = cirris_table(i,:);
    side = curr_line.Side{1};
    time = curr_line.TaskTimer;
    step_mod = curr_line.StepModifier;
    target_num = curr_line.TargetNumber;
    block_num  = floor(target_num/120)+1;

    if strcmp(side,'Right') 

        % Define the tolernace level for the time difference between 
        % the Cirris and the Xsens events--Right foot
        tol=0.3;

        % idx_next_rto = find(time<rto,1);
        % idx_prev_rto = find(time>rto,'last');
        % idx_prev_rhs = find(time>rto,'last');
        % if sum(time>=rhs & (time<=rto))>=1
        % if sum(time>=rhs & (rto-time<=tol))>=1
        %     if sum(time>=rhs & (time<=rhs+tol))==1
        %         idx= find(time>=rhs & (time<=rhs+tol));
        %         gaitEventTable.AdaptiveRight(idx)= step_mod;
        %         gaitEventTable.BlockNumRight(idx) = block_num;
        % 
        %         num_matches_right = num_matches_right + 1;
        %     else
        %         error('More than one matches--check!');
        %     end
        % end

        idx_next_rto = find(time<=rto & abs(time-rto)<median(diff(rto)/2,'omitnan'),1);
        idx_next_lto = find(time<=lto & abs(time-lto)<median(diff(lto)/2,'omitnan'),1);
        idx_next_lhs = find(time<=lhs & abs(time-lhs)<median(diff(lhs)/2,'omitnan'),1);
        % idx_next_rhs = find(abs(time-rhs)<tol,1,'last');
        idx_next_rhs = find(abs(time-rhs)<tol);

        % idx_next_rhs = find(time>=rhs & abs(time-rhs)<tol);
        % idx_next_rhs = find(time>=rhs & abs(time-rhs)<median(diff(rhs)/2,'omitnan'),1,'last');

        if ~isempty(idx_next_rhs) 
            idx = idx_next_rhs;
            gaitEventTable.AdaptiveRight(idx)= step_mod;
            gaitEventTable.BlockNumRight(idx) = block_num;
            count1 = count1 + 1;

            if length(idx_next_rhs)>=2
                sprintf('Multiple indices found. RH index %d',idx_next_rhs)
            end
            if time<rhs(idx_next_rhs)
                count11 = count11 + 1;
                timediff_1 = [timediff_1, rhs(idx_next_rhs)-time];
            else
                count12 = count12 + 1;
                timediff_2 = [timediff_2, time-rhs(idx_next_rhs)];
            end
        elseif ~isempty(idx_next_lto)
            idx = idx_next_lto;
            gaitEventTable.AdaptiveRight(idx)= step_mod;
            gaitEventTable.BlockNumRight(idx) = block_num;
            count3 = count3 + 1;
        elseif ~isempty(idx_next_lhs)
            idx = idx_next_lhs;
            gaitEventTable.AdaptiveRight(idx)= step_mod;
            gaitEventTable.BlockNumRight(idx) = block_num;
            count4 = count4 + 1;
        elseif ~isempty(idx_next_rto) && idx_next_rto>=2
            idx = idx_next_rto-1;
            gaitEventTable.AdaptiveRight(idx)= step_mod;
            gaitEventTable.BlockNumRight(idx) = block_num;
            count2 = count2 + 1;
        end

        % if ~isempty(idx_next_rto) && ~isempty(idx_next_rhs) && (idx_next_rto == idx_next_rhs + 1)
        %     idx = idx_next_rhs;
        %     gaitEventTable.AdaptiveRight(idx)= step_mod;
        %     gaitEventTable.BlockNumRight(idx) = block_num;
        % elseif ~isempty(idx_next_lhs) && ~isempty(idx_next_rhs) && (idx_next_lhs == idx_next_rhs + 1)
        %     idx = idx_next_rhs;
        %     gaitEventTable.AdaptiveRight(idx)= step_mod;
        %     gaitEventTable.BlockNumRight(idx) = block_num;
        % elseif ~isempty(idx_next_rto) && isempty(idx_next_rhs) 
        %     idx = idx_next_rto-1;
        %     gaitEventTable.AdaptiveRight(idx)= step_mod;
        %     gaitEventTable.BlockNumRight(idx) = block_num;
        % elseif isempty(idx_next_rto) && ~isempty(idx_next_rhs) 
        %     idx = idx_next_rhs;
        %     gaitEventTable.AdaptiveRight(idx)= step_mod;
        %     gaitEventTable.BlockNumRight(idx) = block_num;
        % end
       
    elseif strcmp(side,'Left')

        % Define the tolernace level for the time difference between 
        % the Cirris and the Xsens events--Left foot
        % tol=median(diff(lhs),'omitnan');
        tol=0.3;

        % 
        % if sum(time>=lhs & (time<=lhs+tol))>=1
        %     if sum(time>=lhs & (time<=lhs+tol))==1
        %         idx= find(time>=lhs & (time<=lhs+tol));
        %         gaitEventTable.AdaptiveLeft(idx)= step_mod;
        %         gaitEventTable.BlockNumLeft(idx) = block_num;
        % 
        %         num_matches_left = num_matches_left + 1;
        %     else
        %         error('More than one matches--check!');
        %     end
        % end


        idx_next_lto = find(time<=lto & abs(time-lto)<median(diff(lto)/2,'omitnan'),1);
        idx_next_rto = find(time<=rto & abs(time-rto)<median(diff(rto)/2,'omitnan'),1);
        idx_next_rhs = find(time<=rhs & abs(time-rhs)<median(diff(rhs)/2,'omitnan'),1);
        idx_next_lhs = find(abs(time-lhs)<tol);
        % idx_next_lhs = find(time>=lhs & abs(time-lhs)<tol);
        % idx_next_lhs = find(time>=lhs & abs(time-lhs)<median(diff(lhs)/2,'omitnan'),1,'last');
       
        if ~isempty(idx_next_lhs) 
            idx = idx_next_lhs;
            gaitEventTable.AdaptiveLeft(idx)= step_mod;
            gaitEventTable.BlockNumLeft(idx) = block_num;
        elseif ~isempty(idx_next_rto) && idx_next_rto>=2
            idx = idx_next_rto-1;
            gaitEventTable.AdaptiveLeft(idx)= step_mod;
            gaitEventTable.BlockNumLeft(idx) = block_num;
        elseif ~isempty(idx_next_rhs) && idx_next_rhs>=2
            idx = idx_next_rhs-1;
            gaitEventTable.AdaptiveLeft(idx)= step_mod;
            gaitEventTable.BlockNumLeft(idx) = block_num;
        elseif ~isempty(idx_next_lto) && idx_next_lto>=2
            idx = idx_next_lto-1;
            gaitEventTable.AdaptiveLeft(idx)= step_mod;
            gaitEventTable.BlockNumLeft(idx) = block_num;
        end

        % if ~isempty(idx_next_lto) && ~isempty(idx_next_lhs) && (idx_next_lto == idx_next_lhs + 1)
        %     idx = idx_next_lhs;
        %     gaitEventTable.AdaptiveLeft(idx)= step_mod;
        %     gaitEventTable.BlockNumLeft(idx) = block_num;
        % elseif ~isempty(idx_next_lto) && isempty(idx_next_lhs) 
        %     idx = idx_next_lto-1;
        %     gaitEventTable.AdaptiveLeft(idx)= step_mod;
        %     gaitEventTable.BlockNumLeft(idx) = block_num;
        % elseif isempty(idx_next_lto) && ~isempty(idx_next_lhs) 
        %     idx = idx_next_lhs;
        %     gaitEventTable.AdaptiveLeft(idx)= step_mod;
        %     gaitEventTable.BlockNumLeft(idx) = block_num;
        % end
    end
    clear idx
end

%% Fill the zero-gaps in the block number columns 
max_block_num = max([gaitEventTable.BlockNumRight; gaitEventTable.BlockNumLeft]);
for b = 1:max_block_num
    min_ind_right = find(gaitEventTable.BlockNumRight==b, 1);
    max_ind_right = find(gaitEventTable.BlockNumRight==b, 1, 'last');

    min_ind_left = find(gaitEventTable.BlockNumLeft==b, 1);
    max_ind_left = find(gaitEventTable.BlockNumLeft==b, 1, 'last');

    % Fill in the right-foot block numbers
    % Block start
    if min_ind_right<=min_ind_left
        min_ind_right_new = min_ind_right;
    else 
        if rhs(min_ind_left)>lhs(min_ind_left)
            min_ind_right_new = min_ind_left;
        else
            min_ind_right_new = min_ind_left+1;
        end
    end
    % Block end
    if max_ind_right>=max_ind_left
        max_ind_right_new = max_ind_right;
    else 
        if rhs(max_ind_left)<lhs(max_ind_left)
            max_ind_right_new = max_ind_left;
        else
            max_ind_right_new = max_ind_left-1;
        end
    end

    % Update the gaitEventTable for the right foot
    gaitEventTable.BlockNumRight(min_ind_right_new:max_ind_right_new) = b;

    % Fill in the left-foot block numbers
    % Block start
    if min_ind_left<=min_ind_right
        min_ind_left_new = min_ind_left;
    else 
        if lhs(min_ind_right)>rhs(min_ind_right)
            min_ind_left_new = min_ind_right;
        else
            min_ind_left_new = min_ind_right+1;
        end
    end
    % Block end
    if max_ind_left>=max_ind_right
        max_ind_left_new = max_ind_left;
    else 
        if rhs(max_ind_right)<lhs(max_ind_right)
            max_ind_left_new = max_ind_right;
        else
            max_ind_left_new = max_ind_right-1;
        end
    end

    % Update the gaitEventTable for the right foot
    gaitEventTable.BlockNumLeft(min_ind_left_new:max_ind_left_new) = b;
end

end
