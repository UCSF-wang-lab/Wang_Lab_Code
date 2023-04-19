function [band_indicies,band_names] = getCanonicalFreqBandInd(freq_vec)
band_indicies = zeros(8,2);
band_names = {'Delta','Theta','Alpha','Beta','Low Beta','High Beta','Gamma','High Gamma'};

% delta
[~,band_indicies(1,1)] = min(abs(freq_vec-0.5));
[~,band_indicies(1,2)] = min(abs(freq_vec-3));

% theta
[~,band_indicies(2,1)] = min(abs(freq_vec-4));
[~,band_indicies(2,2)] = min(abs(freq_vec-7));

% alpha
[~,band_indicies(3,1)] = min(abs(freq_vec-8));
[~,band_indicies(3,2)] = min(abs(freq_vec-12));

% beta
[~,band_indicies(4,1)] = min(abs(freq_vec-13));
[~,band_indicies(4,2)] = min(abs(freq_vec-30));

% low beta
[~,band_indicies(5,1)] = min(abs(freq_vec-13));
[~,band_indicies(5,2)] = min(abs(freq_vec-20));

% high beta
[~,band_indicies(6,1)] = min(abs(freq_vec-21));
[~,band_indicies(6,2)] = min(abs(freq_vec-30));

% gamma
[~,band_indicies(7,1)] = min(abs(freq_vec-31));
[~,band_indicies(7,2)] = min(abs(freq_vec-60));

% high gamma
[~,band_indicies(8,1)] = min(abs(freq_vec-61));
[~,band_indicies(8,2)] = min(abs(freq_vec-100));

end