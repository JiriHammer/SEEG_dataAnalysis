function defineTriggers(params)
% defines 'triggers' = times and classLabels from trials (pred_raw) based on params.triggering
% saves them into session cache files as 'triggers'
% in more complex cases (emoReco, localizers), for each session:
%   - fill in 'trials_events' (tr x ctg) table of 1/0
%   - apply (sequentially) logical operators on odd columns of 'trials_events' table -> selected trials
%   - extract [t_trigger,clz_label] for selected trials -> 'triggers'

% (c) Jiri, Mar16


%% trigering method: for old scripts = params.paradigm.usedParadigm
if ~isfield(params.triggering, 'method')
    params.triggering.method = params.paradigm.usedParadigm;
end

%% define trigger method = fieldName_fieldVals (e.g.: NUMERICAL COGNITION, LOCAs, I-task)
if strcmp(params.triggering.method, 'fieldName_fieldVals')
    getTriggers_init(params, 'times_labels');           % initialize triggers
    
    % define classes
    nClz = size(params.triggering.classes,1);                   % number of all classes (some may not be found)
    for clz = 1:nClz                                            % classes (= rows in params.triggering.classes)
        params.thisClz = clz;
        nCtg = size(params.triggering.classes{clz,2},1);
        params.nCtg = nCtg;
        getTriggers_init(params, 'trials_events');      % initialize trials_events table
        
        % define events (of several categories)
        for ctg = 1:2:nCtg                                      % go thru event categories
            params.thisCtg = ctg;
            getTriggers_fieldName_fieldVals(params);
        end
        
        % define operators between the 2 consecutive event categories
        getTriggers_init(params, 'trials_class');       % initialize trials_class vector
        for ctg = 2:2:nCtg-1                                    % go thru logical operators
            params.thisCtg = ctg+1;                             % pointer to category to be merged
            assert(strcmp(params.triggering.classes{clz,2}{ctg,1}, 'logical_operator'));
            params.thisOperator = params.triggering.classes{clz,2}{ctg,2};
            getTriggers_applyLogicalOperator(params); 
        end
        
        % get trigger times
        getTriggers_loca_timesLabels(params); 
    end
end

%% EMO-RECO: define trigger method
if strcmp(params.triggering.method, 'emoRecognition')
    getTriggers_init(params, 'times_labels');           % initialize triggers
    
    % define classes
    nClz = size(params.triggering.classes,1);                   % number of all classes (some may not be found)
    for clz = 1:nClz                                            % classes (= rows in params.triggering.classes)
        params.thisClz = clz;
        nCtg = size(params.triggering.classes{clz,2},2);
        params.nCtg = nCtg;
        getTriggers_init(params, 'trials_events');      % initialize trials_events table
        
        % define events (of several categories)
        for ctg = 1:2:nCtg                                      % go thru event categories
            params.thisCtg = ctg;
            switch params.triggering.classes{clz,2}{1,ctg}
               case {'Fear','Happiness','Anger','Sadness','Disgust','Neutral','Violin','Piano'}
                   params.thisCategory = params.triggering.classes{clz,2}{1,ctg};
                   getTriggers_emoReco_stimuliTypes(params);
               case 'Concordant'
                   params.thisCategory = 'same';
                   getTriggers_emoReco_sameDiff(params);                   
               case 'Discordant'
                   params.thisCategory = 'diff';
                   getTriggers_emoReco_sameDiff(params);
               case {'Emotions','Objects'}
                   params.thisCategory = params.triggering.classes{clz,2}{1,ctg};
                   getTriggers_emoReco_objectEmotion(params);                  
               case {'Image','Sound'}
                   params.thisCategory = params.triggering.classes{clz,2}{1,ctg};
                   getTriggers_emoReco_imageSound(params);   
               case 'Answer_correct'
                   params.thisCategory = 'right';
                   getTriggers_emoReco_answer(params);    
               case 'Answer_wrong'
                   params.thisCategory = 'wrong';
                   getTriggers_emoReco_answer(params);   
               case 'Button_left'
                   params.thisCategory = 'L';
                   getTriggers_emoReco_button(params);    
               case 'Button_right'
                   params.thisCategory = 'R';
                   getTriggers_emoReco_button(params); 
               case 'Button_nothing'
                   params.thisCategory = 'n.a.';
                   getTriggers_emoReco_button(params);   
               case 'Paradigm_concordancy'
                   params.thisCategory = 'concordancy';
                   getTriggers_gameType(params);
               case 'Paradigm_freeflow'
                   params.thisCategory = 'freeflow';
                   getTriggers_gameType(params);                         
               otherwise
                  error(['unknown trigger stimuli: ' params.triggering.classes{clz,2}{1,ctg}]);
            end    
        end
        
        % define operators between the 2 consecutive event categories
        getTriggers_init(params, 'trials_class');       % initialize trials_class vector
        for ctg = 2:2:nCtg-1                                    % go thru logical operators
            params.thisCtg = ctg+1;                             % pointer to category to be merged
            params.thisOperator = params.triggering.classes{clz,2}{1,ctg};
            getTriggers_applyLogicalOperator(params); 
        end
        
        % get trigger times
        getTriggers_emoReco_timesLabels(params); 
    end
end
    
%% CAR-GAME: define trigger method
if strcmp(params.triggering.method, 'carDriving')
    switch params.triggering.classes
       case 'rews_obst'
          getTriggers_carGame_rewsObst(params);
       case 'obst_only'
          getTriggers_carGame_obstOnly(params);       
       case 'turns_LR'
          getTriggers_carGame_onTurns(params);                 
       otherwise
          error(['unknown trigger stimuli: ' params.triggering.classes]);
    end
end

%% SSEP: define trigger method
if strcmp(params.triggering.method, 'ssep')
    getTriggers_init(params, 'times_labels');           % initialize triggers
    nClz = size(params.triggering.classes,1);                   % number of all classes (some may not be found)
    for clz = 1:nClz      
        switch params.triggering.classes{clz,2}{1}  % only 'LH' or RH
           case {'LH','RH'}
              params.thisCategory = params.triggering.classes{clz,2}{1};
              getTriggers_ssep_LR(params);
           otherwise
              error(['unknown trigger stimuli: ' params.triggering.classes]);
        end
    end
end

%% LOCALIZERS: define trigger method
if strcmp(params.triggering.method, 'localizers')
    getTriggers_init(params, 'times_labels');           % initialize triggers
    
    % define classes
    nClz = size(params.triggering.classes,1);                   % number of all classes (some may not be found)
    for clz = 1:nClz                                            % classes (= rows in params.triggering.classes)
        params.thisClz = clz;
        nCtg = size(params.triggering.classes{clz,2},2);
        params.nCtg = nCtg;
        getTriggers_init(params, 'trials_events');      % initialize trials_events table
        
        % define events (of several categories)
        for ctg = 1:2:nCtg                                      % go thru event categories
            params.thisCtg = ctg;
            switch params.triggering.classes{clz,2}{1,ctg}
               case {'fruits','animals','faces','houses','letters','objects','scenes','scrambled','words','diff','easy'}
                   params.thisCategory = params.triggering.classes{clz,2}{1,ctg};
                   getTriggers_locaVisu_stimuliTypes(params);
%                case {'diff','easy'}
%                    params.thisCategory = params.triggering.classes{clz,2}{1,ctg};
%                    getTriggers_locaVisu_stimuliTypes(params);    
               case {'loca_arfa','loca_moto','loca_visu','loca_mveb','loca_mvis','loca_mcse'}
                   params.thisCategory = params.triggering.classes{clz,2}{1,ctg};
                   getTriggers_gameType(params);                          
               otherwise
                  error(['unknown trigger stimuli: ' params.triggering.classes{clz,2}{1,ctg}]);
            end    
        end
        
        % define operators between the 2 consecutive event categories
        getTriggers_init(params, 'trials_class');       % initialize trials_class vector
        for ctg = 2:2:nCtg-1                                    % go thru logical operators
            params.thisCtg = ctg+1;                             % pointer to category to be merged
            params.thisOperator = params.triggering.classes{clz,2}{1,ctg};
            getTriggers_applyLogicalOperator(params); 
        end
        
        % get trigger times
        getTriggers_loca_timesLabels(params); 
    end
end

%% SWITCH-EI: define trigger method
if strcmp(params.triggering.method, 'switchEI')
    getTriggers_init(params, 'times_labels');           % initialize triggers
    
    % define classes
    nClz = size(params.triggering.classes,1);                   % number of all classes (some may not be found)
    for clz = 1:nClz                                            % classes (= rows in params.triggering.classes)
        params.thisClz = clz;
        nCtg = size(params.triggering.classes{clz,2},2);
        params.nCtg = nCtg;
        getTriggers_init(params, 'trials_events');      % initialize trials_events table
        
        % define events (of several categories)
        for ctg = 1:2:nCtg                                      % go thru event categories
            params.thisCtg = ctg;
            switch params.triggering.classes{clz,2}{1,ctg}
               case {'I-task','E-task'}
                   params.thisCategory = params.triggering.classes{clz,2}{1,ctg};
                   getTriggers_switchEI_stimuliCondition(params);
               case {'E-I','I-E','E-E','I-I'}
                   params.thisCategory = params.triggering.classes{clz,2}{1,ctg};
                   getTriggers_switching_EI_IE(params);                          
               otherwise
                  error(['unknown trigger stimuli: ' params.triggering.classes{clz,2}{1,ctg}]);
            end    
        end
        
        % define operators between the 2 consecutive event categories
        getTriggers_init(params, 'trials_class');       % initialize trials_class vector
        for ctg = 2:2:nCtg-1                                    % go thru logical operators
            params.thisCtg = ctg+1;                             % pointer to category to be merged
            params.thisOperator = params.triggering.classes{clz,2}{1,ctg};
            getTriggers_applyLogicalOperator(params); 
        end
        
        % get trigger times
        getTriggers_loca_timesLabels(params); 
    end
end

%% I-TASK (questions YES/NO): define trigger method
if strcmp(params.triggering.method, 'I-task')
    getTriggers_init(params, 'times_labels');           % initialize triggers
    
    % define classes
    nClz = size(params.triggering.classes,1);                   % number of all classes (some may not be found)
    for clz = 1:nClz                                            % classes (= rows in params.triggering.classes)
        params.thisClz = clz;
        nCtg = size(params.triggering.classes{clz,2},2);
        params.nCtg = nCtg;
        getTriggers_init(params, 'trials_events');      % initialize trials_events table
        
        % define events (of several categories)
        for ctg = 1:2:nCtg                                      % go thru event categories
            params.thisCtg = ctg;
            switch params.triggering.classes{clz,2}{1,ctg}
               case {'nowFeelings','selfJudgement','selfSemantic','selfEpisodic','nowRelax','addNumbers','controlQuestions'}
                   params.thisCategory = params.triggering.classes{clz,2}{1,ctg};
                   getTriggers_stimuliTypes(params);
               case {'TODO1','TODO2','TODO3'}
                   params.thisCategory = params.triggering.classes{clz,2}{1,ctg};
                   getTriggers_gameType(params);                          
               otherwise
                  error(['unknown trigger stimuli: ' params.triggering.classes{clz,2}{1,ctg}]);
            end    
        end
        
        % define operators between the 2 consecutive event categories
        getTriggers_init(params, 'trials_class');       % initialize trials_class vector
        for ctg = 2:2:nCtg-1                                    % go thru logical operators
            params.thisCtg = ctg+1;                             % pointer to category to be merged
            params.thisOperator = params.triggering.classes{clz,2}{1,ctg};
            getTriggers_applyLogicalOperator(params); 
        end
        
        % get trigger times
        getTriggers_loca_timesLabels(params); 
    end
end

%% numbers: Add or Find7 (questions YES/NO): define trigger method
if ~isempty(strfind(params.triggering.method, 'numbersAddFind7')) % old: strcmp(params.paradigm.usedParadigm, 'numbersAddFind7')
    getTriggers_init(params, 'times_labels');           % initialize triggers
    
    % define classes
    nClz = size(params.triggering.classes,1);                   % number of all classes (some may not be found)
    for clz = 1:nClz                                            % classes (= rows in params.triggering.classes)
        params.thisClz = clz;
        nCtg = size(params.triggering.classes{clz,2},2);
        params.nCtg = nCtg;
        getTriggers_init(params, 'trials_events');      % initialize trials_events table
        
        % define events (of several categories)
        for ctg = 1:2:nCtg                                      % go thru event categories
            params.thisCtg = ctg;
            switch params.triggering.classes{clz,2}{1,ctg}
               case {'findSeven','addition'}
                   params.thisCategory = params.triggering.classes{clz,2}{1,ctg};
                   getTriggers_stimuliTypes(params);
               case {'I=0','I=3'}
                   params.thisCategory = params.triggering.classes{clz,2}{1,ctg};
                   getTriggers_DES_intensity(params);                          
               otherwise
                  error(['unknown trigger stimuli: ' params.triggering.classes{clz,2}{1,ctg}]);
            end    
        end
        
        % define operators between the 2 consecutive event categories
        getTriggers_init(params, 'trials_class');       % initialize trials_class vector
        for ctg = 2:2:nCtg-1                                    % go thru logical operators
            params.thisCtg = ctg+1;                             % pointer to category to be merged
            params.thisOperator = params.triggering.classes{clz,2}{1,ctg};
            getTriggers_applyLogicalOperator(params); 
        end
        
        % get trigger times
        getTriggers_loca_timesLabels(params); 
    end
end

