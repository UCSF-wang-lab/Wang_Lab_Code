function gaitEventTableFixed = fixGaitEvents2(gaitEventTable)

% Make a copy of the gait events table
gaitEventTableFixed = gaitEventTable;

% Estimate gait cycle period
gc_period = median(diff(gaitEventTable.RHS),'omitnan');

% Arrange the events in the right chronological order
for row=1:height(gaitEventTable)
    for col = 1:3
        if ~isnan(gaitEventTable{row,col}) %&& gaitEventTable{row,col}>99
            if col==1 || col==2 || col==3
                if isnan(gaitEventTable{row,col+1}) || gaitEventTable{row,col+1}<gaitEventTable{row,col}
                    for scan_rows=1:height(gaitEventTable)
                        if ~isnan(gaitEventTable{scan_rows,col+1}) && gaitEventTable{scan_rows,col+1}>gaitEventTable{row,col} && (gaitEventTable{scan_rows,col+1}-gaitEventTable{row,col})<1
                            gaitEventTableFixed{row,col+1} = gaitEventTable{scan_rows,col+1};     
                        else
                            gaitEventTableFixed{row,col+1} = NaN;
                        end
                        break;
                    end
                end
            end
        end
    end
end
for row=1:height(gaitEventTable)-1
    if col==4
        if isnan(gaitEventTable{row+1,1}) || gaitEventTable{row+1,1}<gaitEventTable{row,4}
            for scan_rows=1:height(gaitEventTable)-1
                if ~isnan(gaitEventTable{scan_rows,1}) && gaitEventTable{scan_rows,1}>gaitEventTable{row,4} && (gaitEventTable{scan_rows,1}-gaitEventTable{row,4})<1
                    gaitEventTableFixed{row+1,1} = gaitEventTable{scan_rows,1};
                else
                    gaitEventTableFixed{row+1,1} = NaN;
                end
                break; 
            end
        end
    end
end

for row=1:height(gaitEventTable)
    for col = 1:4
        aberrant = 0;
        if ~isnan(gaitEventTable{row,col})
            rest_cols = [1,2,3,4];
            rest_cols(col) = [];

            nonan_cols_lower = find(~isnan(gaitEventTable{row,rest_cols}) & rest_cols<col);
            nonan_cols_higher = find(~isnan(gaitEventTable{row,rest_cols}) & rest_cols>col);

            if (any(gaitEventTable{row,rest_cols(nonan_cols_lower)}>=gaitEventTable{row,col})) || (any(gaitEventTable{row,rest_cols(nonan_cols_higher)}<=gaitEventTable{row,col}))
                aberrant = 1;
            end

            clear nonan_cols_lower nonan_cols_higher rest_cols
        end

        if isnan(gaitEventTable{row,col}) || aberrant
            found_idx = [];
            rest_cols = [1,2,3,4];
            rest_cols(col) = [];

            if any(~isnan(gaitEventTable{row,rest_cols})) %&& row==16 %&& all(gaitEventTable{row,rest_cols}>234)
                nonan_cols_lower = find(~isnan(gaitEventTable{row,rest_cols}) & rest_cols<col);
                nonan_cols_higher = find(~isnan(gaitEventTable{row,rest_cols}) & rest_cols>col);
                     
                if ~isempty(nonan_cols_lower) && ~isempty(nonan_cols_higher)
                    idx1 = find(abs(max(gaitEventTable{row,rest_cols(nonan_cols_higher)})-gaitEventTable{:,col})<gc_period);
                    idx2 = find(gaitEventTable{:,col}<gaitEventTable{row,rest_cols(nonan_cols_higher)});
                    idx3 = find(abs(gaitEventTable{:,col}-min(gaitEventTable{row,rest_cols(nonan_cols_lower)}))<gc_period);
                    idx4 =  find(gaitEventTable{:,col}>gaitEventTable{row,rest_cols(nonan_cols_lower)});
                    intersect_1 = intersect(idx1,idx2);
                    intersect_2 = intersect(idx3,idx4);

                    found_idx = intersect(intersect_1,intersect_2);
                elseif isempty(nonan_cols_lower) && ~isempty(nonan_cols_higher)
                    idx1 = find(abs(max(gaitEventTable{row,rest_cols(nonan_cols_higher)})-gaitEventTable{:,col})<gc_period);
                    idx2 = find(gaitEventTable{:,col}<gaitEventTable{row,rest_cols(nonan_cols_higher)});
                    found_idx = intersect(idx1,idx2);
                elseif isempty(nonan_cols_higher) && ~isempty(nonan_cols_lower)
                    idx1 = find(abs(gaitEventTable{:,col}-min(gaitEventTable{row,rest_cols(nonan_cols_lower)}))<gc_period);
                    idx2 =  find(gaitEventTable{:,col}>gaitEventTable{row,rest_cols(nonan_cols_lower)});
                    found_idx = intersect(idx1,idx2);
                end
            end

            clear nonan_cols_lower nonan_cols_higher rest_cols idx1 idx2 idx3 idx4

            if ~isempty(found_idx)
                if length(found_idx)>1
                    found_idx = found_idx(1);
                end
                if found_idx>row
                    gaitEventTableFixed{row,col} = gaitEventTable{found_idx,col};
                end
            else
                gaitEventTableFixed{row,col} = NaN;
            end       
        end
    end
end

% Remove any lines with all NaN
remove_ind = [];
for row = 1:height(gaitEventTableFixed)
    if sum(isnan(gaitEventTableFixed{row,1:4})) == 4 || sum(isnan(gaitEventTableFixed{row,1:4})) == 3
        remove_ind = [remove_ind,row];
    end
end
gaitEventTableFixed(remove_ind,:) = [];

end