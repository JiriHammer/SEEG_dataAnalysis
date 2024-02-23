function params = SEI_job_classDefinitions(params, clzDefName)
% defines the classes which the triggers of which are extracted from trials
% (c) Jiri, Apr22

%% job-specific settings: E-I vs. I-E
if strcmp('switchin_EI_IE', clzDefName)
    % define paradigm classes
    params.triggering.classes = {...
        'E-I',{'E-I'};
        'I-E',{'I-E'};
        };
end

%% job-specific settings: E-E vs. I-I
if strcmp('switchin_EE_II', clzDefName)
    % define paradigm classes
    params.triggering.classes = {...
        'E-E',{'E-E'};    
        'I-I',{'I-I'};
        };
end

%% job-specific settings: E-E vs. I-I
if strcmp('switchin_IE_EE', clzDefName)
    % define paradigm classes
    params.triggering.classes = {...
        'I-E',{'I-E'};        
        'E-E',{'E-E'};    
        };
end

%% job-specific settings: E-I vs. I-E vs. E-E vs. I-I
if strcmp('switchin_EI_IE_EE_II', clzDefName)
    % define paradigm classes
    params.triggering.classes = {...
        'E-I',{'E-I'};
        'I-E',{'I-E'};        
        'E-E',{'E-E'};    
        'I-I',{'I-I'};
        };
end

%% job-specific settings: E-E vs. I-I
if strcmp('task_E_I', clzDefName)
    % define paradigm classes
    params.triggering.classes = {...
        'E-task',{'E-task'};    
        'I-task',{'I-task'};
        };
    params.triggering.cutPoint = 'onGo';
    params.triggering.time2cut = [-1, 5];
    params.triggering.baseline = [4, 4.75];
    params.plot_triggering.paraTimes = 'rt';
end

%% output directory = clzDefName + car/bip
if ismember('car', params.response.dataProcessing)
    params.storage.outName = [clzDefName '_car'];
elseif ismember('bip', params.response.dataProcessing)
    params.storage.outName = [clzDefName '_bip'];
else
    params.storage.outName = clzDefName;
end
