function getTriggers_applyLogicalOperator(params)
% executes a logical operator between two sets: AND, OR
% set 1 = triggers_clz from previous evaluation (or init)
% set 2 = trials_events of certain category
% saves result in session cache file in variable 'triggers_clz'

% (c) Jiri, Nov16

for sess = 1:size(params.storage.sessionCacheFiles,2)
    clear trials_events triggers_clz
    load(params.storage.sessionCacheFiles{sess}, 'trials_events','triggers_clz');
    
    % logical operation
    switch params.thisOperator       
       case 'AND'
           triggers_clz = triggers_clz & trials_events(:,params.thisCtg);
       case 'OR'
           triggers_clz = triggers_clz | trials_events(:,params.thisCtg);
       otherwise
          error(['unknown operator: ' params.triggering.classes{clz,2}{1,ctg}]);
    end                
    
    % save
    save(params.storage.sessionCacheFiles{sess}, 'triggers_clz', '-append');
end
