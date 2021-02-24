function varargout = addEmptyData(time_vec,LFP_signal,sampling_freq,type)

if ~exist('type','var')
    type = 'blank';
end

out_data = LFP_signal;
out_time = time_vec;

% Find time gaps
time_diff = diff(time_vec);
[diff_val,ind] = uniquetol(time_diff,1e-6);

% Remove time difference that is at the sampling frequency and sort
remove_inds = abs(diff_val-(1/sampling_freq)) < 1e-6;
diff_val(remove_inds) = [];
ind(remove_inds) = [];
[sorted_ind,order] = sort(ind);
sorted_diff_val = diff_val(order);

cum_new_val_added = 0;
for i = 1:length(sorted_ind)
    new_time_vals{i} = (time_vec(sorted_ind(i))+1/sampling_freq:1/sampling_freq:time_vec(sorted_ind(i)+1)-1/sampling_freq)';
    if strcmp(type,'blank')
        new_data_vals{i} = zeros(length(new_time_vals{i}),1);
    elseif strcmp(type,'interp')
        % USE INTERP1
        new_data_vals{i} = zeros(length(new_time_vals{i}),1);
    end
    
    out_time = [out_time(1:sorted_ind(i)+cum_new_val_added);new_time_vals{i};out_time(sorted_ind(i)+cum_new_val_added+1:end)];
    out_data = [out_data(1:sorted_ind(i)+cum_new_val_added);new_data_vals{i};out_data(sorted_ind(i)+cum_new_val_added+1:end)];
    
    cum_new_val_added = cum_new_val_added + length(new_time_vals{i});
end

varargout{1} = out_data;
varargout{2} = out_time;
end