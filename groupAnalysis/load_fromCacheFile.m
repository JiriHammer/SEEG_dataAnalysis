function d_struct = load_fromCacheFile(dataStruct)
% loads results stored in cache file

% (c) Jiri, Oct20

d_struct = struct;

%% cache file
outName = [dataStruct.outName_prefix '_' dataStruct.this_neur '_' dataStruct.this_traj];
cacheFile = [dataStruct.pathBeg filesep outName filesep dataStruct.this_subj filesep 'cacheFile.mat'];
assert(exist(cacheFile,'file') == 2);
d_struct.cacheFile = cacheFile;

%% load varName -> X
disp(['loading: ' dataStruct.var2show ' from: ' cacheFile]);
load(cacheFile, dataStruct.var2show);
X = eval(dataStruct.var2show);

%% extract values: single channel decoding results
if strcmp(dataStruct.var2show, 'decoding')  % decoding{lag,ch}
    nCh = size(X,2);
    nLags = size(X,1);
    d_struct.yVals = nan(nLags,nCh);
    d_struct.yErrs = nan(nLags,nCh);
    d_struct.xVals = nan(nLags,1);
    
    for ch = 1:nCh
        for lag = 1:nLags
            d_struct.yVals(lag,ch) = mean(X{lag,ch}.cc,2);
            d_struct.yErrs(lag,ch) =  sem(X{lag,ch}.cc,2);
            d_struct.xVals(lag) = X{lag,ch}.tau;
        end
    end
end


