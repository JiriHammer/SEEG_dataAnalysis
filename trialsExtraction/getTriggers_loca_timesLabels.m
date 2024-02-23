function getTriggers_loca_timesLabels(params)
% extracts trigger times and class labels, saves/updates as 'triggers'
% saves/updates trials (predictor) structure: pred_prc
% with fields:
%   .wasSelected = 0,1
%   .clzLabel
%   .cutPoint
% selected trials for given class are in triggers_clz
% saves 'triggers' & 'pred_prc' to session cache files

% (c) Jiri, Nov16, Jun20

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
    assert(size(triggers_clz,1) == size(trials,1));
    
    % get trigger times and assign labels
    for tr = 1:size(triggers_clz,1)
        if triggers_clz(tr) == 1
            n_selected = n_selected+1;
        
            % trigger time
            t_trig = getTimeOfTrigger(params, trials(tr));
%             t_trig = [];
%             switch params.triggering.cutPoint
%                case 'onGetReady'
%                    if isfield(trials(tr), 't_prep')
%                         t_trig = trials(tr).t_prep;                
%                    end
%                case 'onShowStimuli'
%                    if isfield(trials(tr), 't_stim')
%                         t_trig = trials(tr).t_stim; 
%                    end
%                case 'onGo'
%                    if isfield(trials(tr), 't_go')
%                         t_trig = trials(tr).t_go;
%                    end
%                case 'onButtonPress'
%                    if isfield(trials(tr), 't_reaction')
%                        if ~isnan(trials(tr).t_reaction)     % only if answer was provided
%                             t_trig = trials(tr).t_reaction; 
%                        end
%                    end                   
%                    
%                otherwise
%                   error(['unknown trigger time point: ' params.triggering.cutPoint]);
%             end
                    
            % cat trigger times & class labels
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
