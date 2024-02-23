function soa_allCh = ch2roi_soa_trialsData(params)
% TBD (now=Oct23)
% returns significance of activations (soa) from trialsData 
% at selected frequency bands for all channels

% (c) Jiri, Nov22

%% list of frequency bands -> list_FB
list_FB = params.soa.selFB;

%% cache file
if ~isfield(params.soa, 'outDir')
    cacheFile = [params.storage.dir_results filesep params.storage.outName filesep params.storage.subjTag filesep 'cacheFile.mat'];
else
    cacheFile = [params.soa.outDir filesep params.storage.subjTag filesep 'cacheFile.mat'];
end
assert(exist(cacheFile,'file') == 2);

%% load structure H & selected channels
clear H selCh_H_resp;
load(cacheFile, 'H', 'selCh_H_resp');    
selCh_H = selCh_H_resp;

%% define freq band
soa_allCh = [];
for freq = 1:size(list_FB,1)
    freqBand = list_FB{freq,1};
    disp([' - loading trialsData, subj = ' params.storage.subjTag ', freq = ' freqBand ' ...']);

    %% load freq. band activations (e.g. trialsData_beta) -> trialsData
    trialsData = [];
    varName = ['trialsData_' freqBand]; 
    % try to load freq. band activation from cache file
    if ~ismember(varName, who('-file', cacheFile))   % if the variable 'varName' exists
        error(['WARNING: subject = ' params.storage.subjTag ', variable = ' varName ' was not found in cacheFile: ' cacheFile]);
    else
        clear trialsData;
        load(cacheFile, varName);
        assert(exist(varName,'var') == 1);
        eval(['trialsData = ' varName ';']); 
        clear(varName); 
    end
    assert(size(selCh_H,2) == size(trialsData.hVals,2));
    
    %% time of interest
    i_t = closestval(trialsData.xVals, params.soa.selTime(1)):closestval(trialsData.xVals, params.soa.selTime(2));
    
    %% was any significant? (over time dimension)
    soa_allCh = cat(1, soa_allCh, any(trialsData.hVals(i_t,:),1));
end

%% was any significant? (over FBs)
soa_allCh = any(soa_allCh,1);
assert(size(selCh_H,2) == size(soa_allCh,2));
