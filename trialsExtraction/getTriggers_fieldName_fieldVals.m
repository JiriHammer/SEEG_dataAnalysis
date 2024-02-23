function getTriggers_fieldName_fieldVals(params)
% determines trials which contained the selected keyword (in params.thisCategory)
% MUST match trials 'tr' with the same keyword in trials(tr).stim_type
% saves result into respective column of 'trials_events' in each sessionCacheFile
%   - rows = trials
%   - values = zeros (event not found) and ones (event found)

% (c) Jiri, Mar16, Apr22

fieldName = params.triggering.classes{params.thisClz,2}{params.thisCtg,1};
fieldVals = params.triggering.classes{params.thisClz,2}{params.thisCtg,2};
display([' - finding events for tuple: ' fieldName]);

%% per session: load trials, find events and save
for sess = 1:size(params.storage.sessionCacheFiles,2)
    
    % load trials
    clear trials pred_raw
    load(params.storage.sessionCacheFiles{sess}, 'pred_raw');
    trials = pred_raw;
    
    % init
    nTr = size(trials,1);
    trial_event = zeros(nTr, 1);  % default: no events found
    
    % find events -> trial_event
    if isfield(trials(1), fieldName)
        for tr = 1:nTr
            if ~isempty(trials(tr).(fieldName))
                if ismember(trials(tr).(fieldName), fieldVals)
                    trial_event(tr) = 1;
                end
            end
        end    
    else
        disp([' - getTriggers: session = ' num2str(sess) ', trials struct does not contain field = ' fieldName]);
    end
    
    % update trials-categories -> trials_event
    clear trials_events
    load(params.storage.sessionCacheFiles{sess}, 'trials_events');
    trials_events(:,params.thisCtg) = trial_event;
    save(params.storage.sessionCacheFiles{sess}, 'trials_events', '-append');
end
