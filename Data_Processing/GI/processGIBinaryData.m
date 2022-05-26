function processGIBinaryData(filePath,savePath)

fid = fopen(filePath);
data = fread(fid,[3,inf],'double');

outTable = array2table(data','VariableNames',{'GI Trial Type','Inter-trial State','Computer Time'});

[~,filename] = fileparts(filePath);
writetable(outTable,fullfile(savePath,[filename,'.csv']));
save(fullfile(savePath,[filename,'.mat']),"outTable");
end