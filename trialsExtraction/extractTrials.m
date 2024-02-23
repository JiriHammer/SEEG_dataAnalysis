function trials = extractTrials(params)
% extract trials from responses based on triggers
% trials =  
%           data: [3073x126x500 double] ~ [samples x channels x trials]
%           time: [3073x1 double]       ~ time in [s], w.r.t. cutting point = 0   
%       sessTime: [3073x500 double]     ~ time in [s] within each session
%         labels: [1x500 double]        ~ trial labels
%     sessNumber: [1x500 double]        ~ session numbers of the trials

% (c) Jiri, Mar16

disp('Extracting trials ...');

%% load sampling rate
clear srate;
load(params.storage.cacheFile, 'srate')
assert(exist('srate','var') == 1);
   
%% define time range to cut in samples (indices w.r.t. the cutting point, here = 0)
i_cut = ceil(params.triggering.time2cut(1)*srate):ceil(params.triggering.time2cut(2)*srate);     % in [samples], w.r.t. cutting point = 0   
time = i_cut'./srate;    % in [s], w.r.t. cutting point = 0   

%% init trials
data = [];
rejected = [];
labels = [];
sessNumber = [];
sessTime = [];
    
%% extract trials...    
for sess = 1:size(params.storage.sessionCacheFiles,2)
    
    % load triggers
    clear triggers classNamesLabels;
    load(params.storage.sessionCacheFiles{sess}, 'triggers', 'classNamesLabels');
        
    % load responses
    clear resp_prc;
    load(params.storage.sessionCacheFiles{sess}, 'resp_prc');
    
    % load session time (from D struct)
    clear timeAxis;
    load(params.storage.sessionCacheFiles{sess}, 'timeAxis'); 
    assert(exist('timeAxis','var') == 1);
    
    % downsample time axis?
    if srate ~= params.init_srate                                            % to dowsample & reject bad recordings epochs
        disp(['downsampling the time axis of session: ' num2str(sess) ' ...']);
        params.downsample.dsUnit = 'fromCacheFile';
        params.thisSess = sess;
        timeAxis = filterData(params, timeAxis, {'downsample'}, 'from_params.init_srate');
    end        
    assert(size(resp_prc,1) == size(timeAxis,1));   
    
    % rejected indices
    clear resp_rej;
    load(params.storage.sessionCacheFiles{sess}, 'resp_rej');
    if ~exist('resp_rej','var')
        warning('No rejection performed.');
        resp_rej = false(size(resp_prc));
    end
    assert(size(resp_prc,1) == size(resp_rej,1));  
    assert(size(resp_prc,2) == size(resp_rej,2));  
    
    % extract trials data
    for tr = 1:size(triggers,2)
        cutPoint = triggers(1,tr);                          % in [s]
        i_cutPoint = closestval(timeAxis, cutPoint);        % in [samples]
        inds = i_cutPoint + i_cut;                          % samples to extract
        if inds(1) > 0 && inds(end) <= size(timeAxis,1)
            data = cat(3, data, resp_prc(inds,:));          % 3D: samples x channels x trials
            rejected = cat(3, rejected, resp_rej(inds,:));  % 3D: samples x channels x trials
            labels = cat(2, labels, triggers(2,tr));        % class labels
            sessNumber = cat(2, sessNumber, sess);          % session number
            sessTime = cat(2, sessTime, timeAxis(inds));    % session time of each trial
        else
            disp(['WARNING: skipping trial (out of session bounds): session = ' num2str(sess) ', trial = ' num2str(tr)]);
        end
    end
    
    disp([' - session ' num2str(sess) ' done.']);
end
assert(size(classNamesLabels,1) == size(params.triggering.classes,1));

%% output trials structure
trials.data = data;
trials.time = time;
trials.rejected = rejected;
trials.sessTime = sessTime;
trials.labels = labels;
trials.sessNumber = sessNumber;
trials.clzNames = classNamesLabels;
trials.freqBandName = 'erp';
