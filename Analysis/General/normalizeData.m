function normalized_data = normalizeData(raw_data,normalization_type,normalization_input)
% normalization_input   [=] [mu,std] column vectors for each frequency.
%                           Number of rows must match the raw data. 

if strcmp(normalization_type,'percent_change')
    normalized_data = (raw_data-normalization_input(:,1))./normalization_input(:,1);   
elseif strcmp(normalization_type,'zscore')
    normalized_data = (raw_data-normalization_input(:,1))./normalization_input(:,2);
end

end