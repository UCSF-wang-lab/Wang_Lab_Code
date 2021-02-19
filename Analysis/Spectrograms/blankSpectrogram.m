function blankSpectrogram(time_vec,cell_spect,spect_time,sampling_freq)

% Go through all ind and blank out spectrogram
spect_time_resolution = uniquetol(diff(spect_time),1e-6);
for i = 1:length(cell_spect)
    for j = 1:length(sorted_ind)
        if time_vec(sorted_ind(j))-spect_time(1) > 0
            
        end
    end
end

end