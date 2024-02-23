function getTriggers_emoReco_imageSound(params)
% finds which stimuli inlcuded Image or Sound
% concordancy paradigm stimuli:
%   Image + Image
%   Image + Sound (!!! sound is always "second" stimuli in trials structure !!!)

% (c) Jiri, Nov16

display([' - finding events for: ' params.thisCategory]);

%% per session: load trials, find events and save
for sess = 1:size(params.storage.sessionCacheFiles,2)
    
    % load trials
    clear trials
    load(params.storage.sessionCacheFiles{sess}, 'pred_raw');
    trials = pred_raw;
    nTr = size(trials,1);
    trial_event = zeros(nTr, 1);  % default: no events found
    
    % find events
    for tr = 1:nTr
        for stim = 1:size(trials(tr).stimulus,2)
            if isfield(trials(tr).stimulus(stim), 'img_snd')
                if strcmp(params.thisCategory, trials(tr).stimulus(stim).img_snd)
                    trial_event(tr) = 1;
                end
            end
        end
    end
    
    % update trials-categories
    clear trials_events
    load(params.storage.sessionCacheFiles{sess}, 'trials_events');
    trials_events(:,params.thisCtg) = trial_event;
    save(params.storage.sessionCacheFiles{sess}, 'trials_events', '-append');
end
