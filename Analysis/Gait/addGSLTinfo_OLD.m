function gaitEventTable = addGSLTinfo_OLD(gaitEventTable,start_trim_time)
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

% Define the tolernace level for the time difference between the Cirris and
% the Xsens events
tol=0.3;
for i = 1:height(cirris_table)
    curr_line = cirris_table(i,:);
    side = curr_line.Side{1};
    time = curr_line.TaskTimer;
    step_mod = curr_line.StepModifier;
    target_num = curr_line.TargetNumber;
    block_num  = floor(target_num/120)+1;

    if strcmp(side,'Right')
        % if sum((time>=rto & time<=rhs)>=1)  
        %     if sum((time>=rto & time<=rhs)==1)
        %         idx= find((time>=rto & time<=rhs));
        %         gaitEventTable.AdaptiveRight(idx)= step_mod;  
        %         gaitEventTable.BlockNumRight(idx) = block_num;
        % 
        %         disp(['case 1 Right, idx',num2str(idx)]);
        %     else
        %         error('More than one matches --check!');
        %     end
        % end

        if sum(time<rto & (time>=rto-tol))>=1
            if sum(time<rto & (time>=rto-tol))==1
                % idx= find(abs(time-rto)<tol);
                idx= find(time<rto & (time>=rto-tol));

                gaitEventTable.AdaptiveRight(idx)= step_mod;
                gaitEventTable.BlockNumRight(idx) = block_num;

                disp(['case 2 Right, idx ',num2str(idx)]);

            else
                error('More than one matches --check!');
            end
        end 

        if sum(time>rhs & (time<=rhs+tol))>=1
            if sum(time>rhs & (time<=rhs+tol))==1
                % idx= find(abs(time-rhs)<tol);
                idx= find(time>rhs & (time<=rhs+tol));
                gaitEventTable.AdaptiveRight(idx)= step_mod;
                gaitEventTable.BlockNumRight(idx) = block_num;

                disp(['case 3 Right, idx',num2str(idx)]);
            else
                error('More than one matches--check!');
            end
        end

    elseif strcmp(side,'Left')
        % if sum((time>=lto & time<=lhs)>=1)  
        %     if sum((time>=lto & time<=lhs)==1)
        %         idx= find((time>=lto & time<=lhs));
        %         gaitEventTable.AdaptiveLeft(idx)= step_mod; 
        %         gaitEventTable.BlockNumLeft(idx) = block_num;
        % 
        %         disp(['case 1 Left, idx ',num2str(idx)]);
        %     else
        %         error('More than one matches--check!');
        %     end
        % end

        if sum(time<lto & (time>=lto-tol))>=1
            if sum(time<lto & (time>=lto-tol))==1
                % idx= find(abs(time-lto)<tol);
                idx= find(time<lto & (time>=lto-tol));
                gaitEventTable.AdaptiveLeft(idx)= step_mod;
                gaitEventTable.BlockNumLeft(idx) = block_num;

                disp(['case 2 Left, idx ',num2str(idx)]);
            else
                error('More than one matches--check!');
            end
        end

        if sum(time>lhs & (time<=lhs+tol))>=1
            if sum(time>lhs & (time<=lhs+tol))==1
                % idx= find(abs(time-lhs)<tol);
                idx= find(time>lhs & (time<=lhs+tol));
                gaitEventTable.AdaptiveLeft(idx)= step_mod;
                gaitEventTable.BlockNumLeft(idx) = block_num;

                disp(['case 3 Left, idx ',num2str(idx)]);
            else
                error('More than one matches--check!');
            end
        end
    end
    
end

%% Fill the zero-gaps in the block number columns 
max_block_num = max([gaitEventTable.BlockNumRight; gaitEventTable.BlockNumLeft]);
for b = 1:max_block_num

end

end

