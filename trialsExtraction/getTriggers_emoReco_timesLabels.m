function getTriggers_emoReco_timesLabels(params)
% extracts trigger times and class labels
% saves to session cache files

% (c) Jiri, Nov16

n_selected = 0;     % number of selected trials
n_total = 0;        % number of all trials
for sess = 1:size(params.storage.sessionCacheFiles,2)
    
    t_trigger = [];
    clz_label = [];
    
    % load trials
    clear trials pred_raw pred_prc
    load(params.storage.sessionCacheFiles{sess}, 'pred_raw', 'pred_prc');
    trials = pred_raw;
    n_total = n_total+size(trials,1);
    
    % load selected trials for given class
    clear triggers_clz
    load(params.storage.sessionCacheFiles{sess}, 'triggers_clz');
    
    % get trigger times and assign labels
    for tr = 1:size(triggers_clz,1)
        if triggers_clz(tr) == 1
            n_selected = n_selected+1;
        
            % trigger times
            t_trig = [];
            switch params.triggering.cutPoint
               case 'onShowStimuli'
                   t_trig = trials(tr).stimulus(1).t_stim;   % hack: time of sound assumed same as those of objects!
               case 'onGetReady'
                   t_trig = trials(tr).t_ready;
               case 'onPause'
                   t_trig = trials(tr).t_pause;
               case 'onPromptAnswer'
                   t_trig = trials(tr).answer.t_choice;
               case 'onButtonPress'
                   t_trig = trials(tr).answer.t_answer;                                           
               otherwise
                  error(['unknown trigger time point: ' params.triggering.cutPoint]);
            end
                    
            % class labels
            if ~isempty(t_trig)
                t_trigger = cat(2, t_trigger, t_trig); 
                clz_label = cat(2, clz_label, params.thisClz);   
                
                % update trials structure
                pred_prc(tr).wasSelected = true;
                pred_prc(tr).clzLabel = params.thisClz;
                pred_prc(tr).clzName = params.triggering.classes{params.thisClz,1};
                pred_prc(tr).cutPoint = t_trig;                
            end
            
        end
    end
    assert(size(t_trigger,2) == size(clz_label,2));
    
    % class names: append class name in sessionCaecheFile
    clzName = params.triggering.classes{params.thisClz,1};
    clear classNamesLabels
    load(params.storage.sessionCacheFiles{sess}, 'classNamesLabels'); 
    classNamesLabels{end+1,1} = clzName;             % class name
    classNamesLabels{end,2} = params.thisClz;       % class label
    save(params.storage.sessionCacheFiles{sess}, 'classNamesLabels', '-append');
    
    % triggers
    clear triggers
    load(params.storage.sessionCacheFiles{sess}, 'triggers');   
    triggers = cat(2, triggers, cat(1, t_trigger, clz_label));
    save(params.storage.sessionCacheFiles{sess}, 'triggers', 'pred_prc', '-append');
end

if n_selected == 0      % no trials were selected
    display(['WARNING: no trials selected for class: ' clzName]);
else
    display([' - selected trials for class: ' clzName ' = ' num2str(n_selected) ' (' num2str(round(100*n_selected/n_total)) ' %)']);
end
