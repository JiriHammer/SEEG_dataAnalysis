function paraTimes = getParadigmTimes(params)
% returns paradigm times defined in 'params.plot_triggering.paraTimes'
% selected trials in pred_prc(tr).wasSelected (see: getTriggers_..._timesLabels.m)
% paraTimes = mean over all selected trials
% paraTimes = 2D: [class labels x selected time points]

% (c) Jiri, Jun20

if ~isfield(params.plot_triggering, 'paraTimes')
    paraTimes = [];
    return;
end

if isempty(params.plot_triggering.paraTimes)
    paraTimes = [];
    return;
end

%% define paradigm times (mean over selected trials)        
nClz = size(params.triggering.classes,1);    % number of all classes (some may not be found)  
t_all = [];
for sess = 1:size(params.storage.sessionCacheFiles,2)
    % load predictor (trials): pred_prc (see: getTriggers_..._timesLabels.m)
    clear pred_prc classNamesLabels
    load(params.storage.sessionCacheFiles{sess}, 'pred_prc', 'classNamesLabels'); 
    assert(nClz == size(classNamesLabels,1));

    % extract reaction times from selected trials
    for tr = 1:size(pred_prc,1)
        if isfield(pred_prc(tr), 'wasSelected')
            if pred_prc(tr).wasSelected
                t_tr = [];
                if ismember('rt', params.plot_triggering.paraTimes)
                    if strcmp(params.paradigm.usedParadigm, 'emoRecognition')
                        t_tr = cat(1, t_tr, pred_prc(tr).answer.t_answer);       % time of reaction
                    else
                        t_tr = cat(1, t_tr, pred_prc(tr).t_reaction);       % time of reaction
                    end
                end
                if ismember('go', params.plot_triggering.paraTimes)
                    if strcmp(params.paradigm.usedParadigm, 'emoRecognition')
                        t_tr = cat(1, t_tr, pred_prc(tr).stimulus(1).t_stim); % time of stimulus (go! cue)
                    else
                        t_tr = cat(1, t_tr, pred_prc(tr).t_go);             % time of go! cue
                    end
                end
                if ismember('nextTrial_rt', params.plot_triggering.paraTimes)
                    t_tr = cat(1, t_tr, pred_prc(tr+1).t_reaction);     % next trial: time of reaction
                end

                t_tr = t_tr - pred_prc(tr).cutPoint;            % w.r.t. cutPoint
                t_tr = cat(1, t_tr, pred_prc(tr).clzLabel);     % last row = class label !!!
                t_all = cat(2, t_all, t_tr);
            end
        end
    end
end

%% paradigm times (mean over trials for each class)
paraTimes = nan(nClz, size(t_all,1)-1);                     % 2D: class labels x selected time points
for clz = 1:nClz
    i_clz = find(t_all(end,:) == classNamesLabels{clz,2});  % last row = class label !!!
    if ~isempty(i_clz)                
        paraTimes(clz,:) = nanmean(t_all(1:end-1,i_clz),2)';
    end
end

