function params = changeDir_outputResults(params, dir_outputResults)
% changes paths of output results to a new directory specified in dir_outputResults
% rewrites: params.storage.dir_results
% the dir_outputResults does not include the outName (nasme of the output
% directory), just the root directory
% for example:  
% dir_outputResults = 'D:\outputResults';
% dir_outputResults = 'G:\dox\proj_switching_EI';
% useful in case data were moved (or copied) elsewhere than in original job

% (c) Jiri, Oct21

%% change paths to output directory of results
if ~strcmp(dir_outputResults, params.storage.dir_results)
    disp(['WARNING: changing outputResults to: ' dir_outputResults]);    
    params.storage.dir_results = dir_outputResults;     % rename dir results
    assert(exist([params.storage.dir_results filesep params.storage.outName],'dir') == 7);
end
disp(['Output dir: ' params.storage.dir_results filesep params.storage.outName]); 
