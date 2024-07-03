function gaitEventTable = fixGaitEvents(gaitEventTable,gaitEventOrder)

% If there are problematic instances, i.e. gait events are not arranged in
% chronological order, fix that first.
for row=2:height(gaitEventTable)
    for col = 1:4
        other_cols = setdiff([1 2 3 4],col);
        if any(gaitEventTable{row,other_cols} - gaitEventTable{row,col}<0)
            gaitEventMat = [gaitEventTable{1:row,:}; 
                              NaN NaN NaN NaN;
                              gaitEventTable{(row+1):end,:}];
            gaitEventTable = array2table(gaitEventMat,'VariableNames',gaitEventOrder);

            gaitEventTable{row+1,col} = gaitEventTable{row,col};
            gaitEventTable{row,col} = NaN;
        end
    end
end

% Organize approprietly and get rid of unnecessary NaNs
for row=2:height(gaitEventTable)
    for col = 1:4
        if isnan(gaitEventTable{row,col})
            next_rowVal = NaN;
            next_row = row;
            while isnan(next_rowVal) && next_row+1<=height(gaitEventTable)
                next_row = next_row + 1;
                next_rowVal = gaitEventTable{next_row,col};
            end

            %Check
            if col==2 || col==3
                if next_rowVal>gaitEventTable{row,col-1} && abs(next_rowVal-gaitEventTable{row,col-1})<1 %&& nextrow_val<gaitEventTable{row,col+1}
                    gaitEventTable{row,col} = next_rowVal;
                    gaitEventTable{next_row,col} = NaN;
                end
            end

            if col==4
                if next_rowVal>gaitEventTable{row,col-1} && abs(next_rowVal-gaitEventTable{row,col-1})<1 %&& nextrow_val<gaitEventTable{row,col+1}
                    gaitEventTable{row,col} = next_rowVal;
                    gaitEventTable{next_row,col} = NaN;
                end
            end

            if col==1
                if next_rowVal>gaitEventTable{row-1,4} && abs(next_rowVal-gaitEventTable{row-1,4})<1
                    gaitEventTable{row,col} = next_rowVal;
                    gaitEventTable{next_row,col} = NaN;
                end
            end
        end
    end
end

% Remove any lines with all NaN
remove_ind = [];
for row = 1:height(gaitEventTable)
    if sum(isnan(gaitEventTable{row,1:4})) == 4 || sum(isnan(gaitEventTable{row,1:4})) == 3
        remove_ind = [remove_ind,row];
    end
end
gaitEventTable(remove_ind,:) = [];

end