function subjList = subjList_updateAnalyzed(params)
% updates subject list based on those that were analyzed (e.g. if some paradigm parts could not be run)

% (c) Jiri, Jul22

if exist([params.storage.dir_results filesep params.storage.outName],'dir') ~= 7
    error(['output directory for results not found, dir = ' params.storage.dir_results filesep params.storage.outName]);
end

disp('Pooling analysis: updating subject list ...');
freqBands = cat(1, {'erp', [nan, nan]}, params.triggering.freqBands);   % add ERP to freq. bands
subjList = [];
n = 1;
for subj = 1:size(params.storage.subjList,1)
    subjTag = params.storage.subjList{subj,1};
    includeThisSubj = true;

    % cache file
    cacheFile = [params.storage.dir_results filesep params.storage.outName filesep subjTag filesep 'cacheFile.mat'];
    if exist(cacheFile,'file') ~= 2
        disp([' - WARNING: cache file not found: subj = ' subjTag ', cache file = ' cacheFile]);
        continue;
    end
    
    % load freq. band activations
    for freq = 1:size(freqBands,1)
        
        % is freq. band activation (e.g. trialsData_beta) in cache file?
        varName = ['trialsData_' freqBands{freq,1}];  
        if ismember(varName, who('-file', cacheFile))   % if the variable 'varName' exists

            % check if all classes were extracted!
            if freq == 1
                load(cacheFile, varName);
                trialsData = eval(varName);
                if size(trialsData.yVals,3) ~= size(params.triggering.classes,1)
                    disp([' - WARNING: some classes were not found for subj = ' subjTag ', nClz = ' num2str(size(trialsData.yVals,3))]);
                    includeThisSubj = false;
                    break;
                end
                clear trialsData;
                clear(varName);
            end            
        else
            includeThisSubj = false;
        end
    end

    % includeThisSubj ?
    if includeThisSubj
        subjList{n,1} = subjTag;
        n = n+1;   
    else
        disp([' - WARNING: excluding subj = ' subjTag]);
        if exist('trialsData','var')
            disp(['    - trialsData found, but nClz = ' num2str(size(trialsData.yVals,3))]);
        else
            disp('     - no trialsData found.');
        end
    end
end
subjList = unique(subjList);
if isempty(subjList)
    disp('WARNING: no data found!');
end
