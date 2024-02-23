function i_base = stftBaselineIndices(params, nSamples_sess)
% returns session indices used for STFT baseline
% nSamples_sess = number of samples of the processed session
% (c) Jiri, Aug16

% indices of game state "search"
if strfind(params.paradigm.usedParadigm, 'uboat')
    gameState = 'search';
    
    % load "raw" (unprocessed) paradigm states (-> paradigmStates)
    clear time_raw paradigmStates paradigmSettings;
    load(params.storage.sessionCacheFiles{params.thisSess}, 'time_raw', 'paradigmStates', 'paradigmSettings');
    assert(size(time_raw,1) == size(paradigmStates,1));

    % downsample? (-> gameStates)
    timeStep = params.stft_freq.timeStep;
    params.downsample.dsRate = timeStep;
    gameStates = filterData(params, paradigmStates, {'downsample'});
    assert(size(gameStates,1) == nSamples_sess);
    
    % get rejected indices
    st_search = getParadigmStateCode(gameState, 'code', paradigmSettings.stateCode_table);
    i_base = find(gameStates == st_search);

else
    i_base = 1:nSamples_sess;       % default: all samples
end

