function t_trig = getTimeOfTrigger(params, trialStruct)
% returns trigger time of specific trial in trialStruct

% (c) Jiri, Jun20

t_trig = [];
switch params.triggering.cutPoint
   case 'onGetReady'
       if isfield(trialStruct, 't_prep')
            t_trig = trialStruct.t_prep;                
       end
   case 'onShowStimuli'
       if isfield(trialStruct, 't_stim')
            t_trig = trialStruct.t_stim; 
       end
   case 'onGo'
       if isfield(trialStruct, 't_go')
            t_trig = trialStruct.t_go;
       end
   case {'onButtonPress','onBP'}
       if isfield(trialStruct, 't_reaction')
           if ~isnan(trialStruct.t_reaction)     % only if answer was provided
                t_trig = trialStruct.t_reaction; 
           end
       end                   
   case 'onHold'
       if isfield(trialStruct, 't_hold')    % in MVEB: hold period = hold the info after text disappeared
            t_trig = trialStruct.t_hold;                
       end
       
   otherwise
      error(['unknown trigger time point: ' params.triggering.cutPoint]);
end