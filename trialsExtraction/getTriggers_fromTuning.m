function triggers = getTriggers_fromTuning(params)
% uses as triggers the binned predictor from tuning analysis

% (c) Jiri, Dec12

predCh = 1;
position = params.tuning.predictor.selectedChnls{predCh};               % TO DO: for multiple predictors...
selCh = getDecodingChannels(params.connectionTable, params.tuning.predictor.signalType, position);
info.selCh = setdiff(selCh, params.tuning.predictor.excludedChnls);    % expects non-empty row vector of selected channels
info.rejCh = params.tuning.predictor.excludedChnls;
info.name = params.tuning.predictor.signalType;
info.chPos = position;   
info.nGroups = 1;                                                       % TO DO: for mulitple groups...    

%% load predictor signal
assert(params.tuning.doTuning);
load(params.storage.cacheFile, 'binPred', 'sampleRateAfterFiltering');
assert(params.amp.srate == sampleRateAfterFiltering);

%% bin values
if strcmp(params.tuning.rejection{1}, 'excludeZeroVelocity')
    binVals = cat(1, -666, binPred.binValues);
else
    binVals = binPred.binValues;
end

nSess = length(params.simulation.trainSession);
assert(size(binPred.binsArray,2) == nSess);
sessSizes = nan(nSess,1);
i_trigs = [];                                                           % indices  of triggers 
t_trigs = [];                                                           % time  of triggers, in [s]
d_trigs = [];                                                           % duration of triggers, in [samples]
i_sess  = [];                                                           % session indices
i_clz   = [];                                                           % class labels
for sess = 1:nSess
    sessSizes(sess) = size(binPred.binsArray{sess},1);
    sessTime = [1:sessSizes(sess)]'./params.amp.srate;                  % in [s], 0 = sessBeg

    % extract tuning "events"
    for b = 1:size(binVals,1)
        thisBin = binVals(b);
        [i_event, t_event] = getContValsPatches(binPred.binsArray{sess}, thisBin, 'beg');           % beg & duration of pause, in [samples]
        
        % pool over sessions
        i_trigs = cat(1, i_trigs, i_event);                                 % event sample index 
        t_trigs = cat(1, t_trigs, sessTime(i_event));                       % time of event, in [s]
        d_trigs = cat(1, d_trigs, t_event);                                 % event duration, in [samples]
        i_sess  = cat(1, i_sess, sess*ones(size(i_event,1),1));
        i_clz   = cat(1, i_clz, thisBin*ones(size(i_event,1),1));
    end
end

%% output, TO DO:
triggers{predCh}.sampleSessLabel = cat(2, i_trigs, i_sess, i_clz);      % mandatory field, triplets = [sample, sess, clz]
triggers{predCh}.timeSessLabel = cat(2, t_trigs, i_sess, i_clz);        % mandatory field, triplets = [time, sess, clz], time in [s] w.r.t. sessBeg = 0
triggers{predCh}.eventDuration = d_trigs./params.amp.srate;                 % event duration, in [s]
triggers{predCh}.info = info;
triggers{predCh}.info.binning = binPred.tunInfo;
triggers{predCh}.info.sessSizes = sessSizes;                                % mandatory field
triggers{predCh}.info.dataType = 'timeSeries';                              % mandatory field ?
triggers{predCh}.info.clzColors = colorPalette(length(binPred.binValues));
