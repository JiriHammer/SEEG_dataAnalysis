function params = SEI_job_commonSettings()
% trigger analysis: SWITCH-EI
% class 1: switching: E -> I
% class 2: switching: I -> E
% computes spectra first (baseline = entire sessions)
% extracts trials (4D: freq x time x channels x trials)
% (c) Jiri, May20

%% default params
clear params;
params_default;
params.storage.dir_results = 'G:\dox\ms_switch_EI\data\v2_stft_sessions';

%% predictor (= trials)
params.predictor = struct(...                           % predictor signal info (used in regression, tuning, triggering, ..)
    'D_field', 'trials', ...                                % MUST match a field in D-struct
    'signalType', 'paradigm', ...                            % selects channels of a signal type of interest
    'rejection', {{''}}, ...                                % rejection of data?
    'dataProcessing', {{'raw'}}, ...    % preprocessing of predictor signal
    'saveName', 'pred' ...
    );
params.load_D = struct(...                  % field names to load from D-struct
    'fieldNames', {{'timeAxis','gameType','trials'}} ...
    );

%% response (= iEEG) processing filters
params.stft_freq = struct(...           % short time Fourier transform, avg over freq. band
    'windowSize', 0.50, ...                  % in [s]
    'timeStep', 0.10, ...                    % in [s]
    'log10_trafo', true, ...                % log10 of PSD ?
    'allFreqBands', [], ...                 % field for listing multiple frequency bands to be extracted during 1 run of STFT computation
    'freqBand', [55 150] ...                  % stft, amp, avg over freq band: [loF, hiF] in Hz  
    );
params.bp_freq = struct(...
    'filterType', 'band', ...
    'butterOrder', 3, ...                   % butterworth filter order
    'freqBand', [8 13] ...                 % band pass: [loF, hiF] in Hz  
    );

%% response (= iEEG)
params.response = struct(...                            % predictor signal info (used in regression, tuning, triggering, ..)
    'D_field', 'ampData', ...                                   % MUST match a field in D-struct
    'signalType', 'iEEG', ...                                   % selects channels of a signal type of interest
    'rejection', {{'chnls_broken','chnls_epi_IIA','chnls_epi_SOZ','normOutliers'}}, ...    % rejection of data?
    'dataProcessing', {{'car','bp_hiPass'}}, ...              % preprocessing of predictor signal
    'saveName', 'resp' ...
    );
params.rejection.normOutliers = struct(...   % rejects outliers from normal distributions (larger than 6 sigma)
    'sigmaThr', 5 ...                       % in [std]
    );

%% triggering routine
params.triggering = struct(...
    'doTriggering', true, ...
    'response', params.response, ...        % extracted (response) signal 
    'classes', 'defined below', ...         % defines the classes, options: see defineTriggers.m
    'cutPoint', 'onButtonPress', ...                 % choices: onGo, onButtonPress, onShowStimuli
    'time2cut', [-4, 4], ...                % in [s], w.r.t. cutPoint (before & after trigger time)
    'doSpectra', 'stft_sessions', ...                % do spectral analysis, choices: stft, mtft, none
    'baseline', [-4.0, -3.0], ...           % in [s], w.r.t. cutPoint 
    'freqBands', [], ...                     % defined below
    'behaviour',{{'RT','RWM'}} ...          % behaviour from paradigm: RT = reaction times, RWM = right/wrong/missed
    );
params.triggering.freqBands = {...
    'delta',[0, 3];
    'theta',[4, 7];
    'alpha',[8, 12];
    'beta',[13, 30];
    'loGamma',[30, 45];
    'hiGamma',[55, 120];
    };
params.plot_triggering = struct(...
    'time2plot', [-3.5, 3.5], ...           % time to plot, in [s] w.r.t. cutting point (= 0 s)
    'verLines', [-3:3], ...                 % vertical dashed lines
    'horLines', [0], ...                    % horizontal dashed lines
    'brainTopo', false, ...                 % brain topology in 3D transparent brain
    'trials_avg', true, ...                 % trials means
    'spectra_avg', true, ...                % spectra means
    'trials_all', false, ...                % all extracted trials
    'spectra_all', false, ...               % TO DO: spectra all trials
    'stackedTrials', false, ...             % stacked trials for each channel
    'paraTimes', {{'go','nextTrial_rt'}} ...          % paradigm times: e.g. (reaction times), choices: go, rt, nextTrial_rt  
    );

%% job-specific settings: loca_visu: objects x faces
% define paradigm type
params.paradigm.usedParadigm = 'switchEI';
params.paradigm.specificType = 'switchEI';

%% set subjects
[params.storage.pathBeg, params.storage.subjList] = get_subjectList(params);
params.storage.subjList = {
%         '20_PR3'; ...  
...        '20_PR4'; ...    % crashes!
%         '20_PR5'; ...
%         '20_PR6'; ...
%         '20_PR7'; ...
%         '20_PR8'; ...
%         '20_PR9'; ...    
%         '20_PR10'; ... 
%         '20_PR11'; ... 
%         '20_PR12'; ...
...        '20_PR13'; ...   % crashed in 4clz
        '20_PR14'; ...   % crashed
        '21_PR1'; ...
        '21_PR2'; ...
        '21_PR4'; ...
        '21_PR5'; ...
        '21_PR6'; ... 
        '21_PR7'; ... 
        '22_PR1'; ... 
    };
