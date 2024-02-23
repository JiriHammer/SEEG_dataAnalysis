% parameters file for iEEG experiments
% previous: paramsEcog.m

% (c) Jiri, Dec10 - now!   :o)

params = struct;

%% set paths
localComputerSettings;
addpath(genpath(params.storage.dir_neuroHammer));
dbstop if error;

%% paradigm information
params.paradigm = struct(...
    'usedParadigm', 'carDriving', ...    % options: 'carDriving', 'emoRecognition', 'uboat'
    'specificType', 'all' ...
    );

params.load_D = struct(...              % field names to load from D-struct
    'fieldNames', {{'fieldName_1','fieldName_2','etc.'}} ...
    );
% for example:
% params.load_D = struct(...              % field names to load from D-struct
%     'fieldNames', {{'tracker'}} ...
%     );

%% quality checks
params.qualityChecks = struct(...
    'trajectoryChecks', false, ...                  % TO DO: not ready to use ...
    'doSignalChecks', false, ...                    % session means, FT, autocorrelation, cross-correlation
    'loadIntoEeglab', false, ...                    % option to load data into EEGLAB for visual checks
    'electrodeLocalization', false, ...             % for iEEG signal only, plots MRI + CT + MNI coors of each channel
    'variable2check', 'resp' ...                    % options: resp_prc, resp_raw, pred_prc, pred_raw
    );

%% filtering: list of different filters for feature extraction
params.raw = struct(...                 % raw data per se
    );
params.zeroMean = struct(...            % norm. to zero mean
    );
params.zeroOut_negative = struct(...    % zero-out all negative values
    );
params.zeroOut_positive = struct(...    % zero-out all positive values
    );
params.norm2std = struct(...            % norm. to std. deviation = 1
    );
params.norm2iqr = struct(...            % norm. to inter-quartile range = 1
    );
params.z_score = struct(...             % norm. to zero mean & std. deviation = 1
    );
params.downsample = struct(...          % only downsamples (takes every n-th sample)
    'dsUnit', 'fromCacheFile', ...          % choices: samples, seconds, fromCacheFile (loads timeStep) 
    'dsRate', NaN ...                       % in [samples]
    );
params.decimate = struct(...            % anti-aliasing (low-pass filters) & downsamples
    'dsUnit', 'fromCacheFile', ...          % choices: samples, seconds, fromCacheFile (loads timeStep) 
    'dsRate', NaN ...                       % in [samples]
    );
params.resample = struct(...            % resamples (interpolates and decimates)
    'new_srate', 500, ...                   % new sampling rate, in [Hz]
    'filt_order', 100 ...                   % filter order, filt_order = 0 => downsampling
    );
params.fft = struct(...                 
    );
params.sgolay = struct(...              % Savitzky-Golay low pass filter
    'order', 2, ...
    'wSize', 0.25 ...                      % sgolay window size, in [s]
    );
params.asym_sgolay = struct(...         % asymmetric Savitzky-Golay low pass filter
    'order', 2, ...
    'wSize', 1.5, ...                       % sgolay window size, in [s]
    'rOffset', 0.05 ...                     % asymmetry point offset from right side of the window, in [s]
    );
params.smooth_Hann = struct(...         % smoothing (convolution) with Hann window
    'wSize', 0.50 ...                       % Hann window size, in [s]
    );
params.bp_freq = struct(...             % filtering using filtfilt.m
    'filterType', 'band', ...               % band-pass filter
    'butterOrder', 5, ...                   % butterworth filter order
    'freqBand', [60 64] ...                 % band pass: [loF, hiF] in Hz  
    );
params.bp_hiPass = struct(...           % filtering using filtfilt.m
    'filterType', 'high', ...               % high-pass filter
    'butterOrder', 3, ...                   % butterworth filter order
    'freqBand', [0.10 Inf] ...              % high pass: [hiCut, Inf] in Hz  
    );
params.bp_loPass = struct(...           % filtering using filtfilt.m
    'filterType', 'low', ...                % low-pass filter
    'butterOrder', 4, ...                   % butterworth filter order
    'freqBand', [0 4] ...                   % low pass: [loF, hiF] in Hz 
    );
params.bp_notch = struct(...            % filtering using filtfilt.m            
    'filterType', 'stop', ...               % notch filter @ 50 Hz
    'butterOrder', 3, ...                   % butterworth filter order
    'freqBand', [48 52] ...                 % low pass: [loF, hiF] in Hz     
    );
params.acausal_filter = struct(...      % filtering using filter.m
    'filterType', 'high', ...               % filter type, choices: high, low, stop, band
    'butterOrder', 3, ...                   % butterworth filter order
    'freqBand', [0.10 Inf] ...              % high pass: [hiCut, Inf] in Hz  
    );
params.stft_freq = struct(...           % short time Fourier transform, avg over freq. band
    'windowSize', 2.0, ...                  % in [s]
    'timeStep', 0.1, ...                    % in [s]
    'log10_trafo', true, ...                % log10 of PSD ?
    'allFreqBands', [], ...                 % field for listing multiple frequency bands to be extracted during 1 run of STFT computation
    'freqBand', [8 13] ...                  % stft, amp, avg over freq band: [loF, hiF] in Hz  
    );
params.stft_freqComp = struct(...       % short time Fourier transform -> complex frequency components
    'windowSize', 2.0, ...                  % in [s]
    'timeStep', 0.01, ...                   % in [s]
    'freqBins', 1:10, ...                   % obligatory field name. bin indices of STFT, corr. to freq = 1/windowSize*freqBins
    'keepDim', 'complex' ...                % options: 'complex', 'amp', 'phase'
    );
params.stft_freqComp_pred = struct(...  % short time Fourier transform -> complex frequency components (specific for predictors)
    'windowSize', 2.0, ...                  % in [s]
    'timeStep', 0.05, ...                   % in [s]
    'freqBins', 1:10, ...                   % obligatory field name. bin indices of STFT, corr. to freq = 1/windowSize*freqBins
    'keepDim', 'complex' ...                % options: 'complex', 'amp', 'phase'
    );
params.stft_freqComp_resp = struct(...  % short time Fourier transform -> complex frequency components (specific for responses)
    'windowSize', 2.0, ...                  % in [s]
    'timeStep', 0.05, ...                   % in [s]
    'freqBins', 1:10, ...                   % obligatory field name. bin indices of STFT, corr. to freq = 1/windowSize*freqBins
    'keepDim', 'complex' ...                % options: 'complex', 'amp', 'phase'
    );
params.mtft_freq = struct(...           % multi taper (short time) Fourier transform
    'windowSize', 0.5, ...                  % in [s]  
    'timeStep', 0.1, ...                    % in [s]    
    'freqBand', [8 13], ...                 % stft, amp, avg over freq band: [loF, hiF] in Hz  
    'getSpectra', 'loadFromFile', ...       % options: 'computeAndSave', 'loadFromFile'
    'fileName', '' ...                      % for load|save: /full/path/to/fileDir + automatically adds: /session_N/mtft_spectra.mat
    );
params.hilb = struct(... 
    'keepDim', 'complex' ...                % options: 'real', 'imag', 'complex', 'amp', 'phase'
    );
params.ampl = struct(... 
    );
params.n_power = struct(...             % Y = X^n    
    'n', 2 ...
    );
params.addValue = struct(...             % Y = X + c    
    'val2add', 0 ...                        % = c
    );
params.angle = struct(... 
    );
params.cosine = struct(... 
    );
params.log10 = struct(... 
    );
params.direction = struct(...   % direction (of movement), returns x (& y) components od directional unit vector
    );
params.random = struct(... 
    );
params.channelSelection = struct(...            % selects channels
    'channels', 'all' ...                    % which channels to select
    );
params.sineWave = struct(...
    'frequency', 10, ...         % in [Hz]
    'phase', 0, ...             % in [rad]
    'amplitude', 1 ...          % amplitude of sine wave
    );
params.runavg = struct(...              % running average detrending
    'wSize', 0.01 ...                        % wsize in [s]  
    );
params.car = struct(...
    'chGroups', 'perHeadbox' ...        % CAR: channel groups = 'perHeadbox', 'perElectrode'
    );
params.bip = struct(...                 % bipolar reference
    'chGroups', '12_23' ...        % CAR: channel groups = '12_23' (all neighbours), '12_34' (next neighbours)
    );
params.bip_eog = struct(...                 % bipolar reference for EOG
    );
params.nan = struct(...
    );
params.pca = struct(...
    'preprocess', {{'car', 'sg'}} ...   % preprocessing of data for filter design 
    );
params.subtractEcg = struct(...         % subtracts ECG-related signal
    'getResponse', true, ...
    'preprocess', {{'car','bp_hiPass', 'bp_freq', 'norm2std'}} ...  % preprocessing of response signal
    );
params.subtractZeroSpeed = struct(...   % to delete: subtracts mean ECoG corresp. to zero speed
    'loadFromFile', true, ...
    'filePath', '/export/jiritmp/webdavmirror/zeroSpeedBaseline' ...  % preprocessing of response signal
    );
params.derivative = struct(...
    'nPointStencils', 5 ...             % choices: 1 (then diff is used), 3, 5, 7, 9
    );
params.plot = struct(...
    );
params.save2cache = struct(...          % stores current 'filterData' into sessionCacheFile
    'addTag', [] ...                    % suffix to add (a string)
    );
params.laplace = struct(...             % computes laplacian transformation of each sample
    );
params.surrogate = struct(...           % computes surrogate timeseries (FFT, random phases, iFFT)
    );
params.shift_subtract = struct(...       % circshift the signal by shiftLag
    'timeLag', 0.5, ...                   % in [s], positive => shift in the future, negative = shift into the past
    'withSubtraction', true ...                % subtract the original (unshifted) signal
    );
params.detect_threshold = struct(...       % detect samples when activity crosses a threshold
    'threshold', 0.05, ...                   % in [std] of the signal. x(t) > thr -> t = 1
    'edge', 'L' ...                          % L (left), R (right), B (both)
    );
params.combine2complex = struct(...     % combines real 'x(t)' with 'y(t)' to complex: x(t) + iy(t)
    'chnl_1', 1, ...                       % specifies 1st channel 'x' (e.g. x-coor)
    'chnl_2', 2  ...                       % specifies 2nd channel 'y' (e.g. y-coor)
    );
params.getDir = struct(...              % combines selChnls to complex, normalizes to unit length and extracts real & imag parts <-> direction
    'signalType', 'Traj', ...
    'selChnls', {{'xVel','yVel'}}, ...
    'speedThr', 0.05 ...
    );
params.trafo2trinary = struct(...       % discrete "trinary signal" {-1, 0, 1}
    'threshold', 0.1 ...                   % in [std] of the signal. x(t) > thr -> +1, x(t) < thr -> -1, else -> 0
    );
params.trafo2artificialChnl = struct(...    % Scales xVel(left) by 1000, xVel(right) by 0.001
    'signalType', 'Traj', ...
    'chName', 'xVel'  ...    
    );
params.fakeTriggers = struct(...        % creates fake random signal for checking the expected output
    'spikes_chnls', 1:10, ...                 % channels, where to insert spikes = single value at trigger point
    'activityHigh_chnls', 11:20, ...          % channels, where to insert activation = higher activity
    'activityLow_chnls', 21:30  ...           % channels, where to insert de-activation = lower activity 
    );
params.stft2magnInfo = struct(...           % to delete: test triviality of ms1. Scales STFT magnitude by velocity
    );
params.model_copyNoise = struct(...     % model = copy + noise. Typically: copy = xVel, noise = whiteNoise (~ "copy+noise" model)
    'noise', 'whiteNoise', ...                  % options: 'whiteNoise', (to do: 'ecogNoise')
    'copy', {{'absVel','absPos'}}, ...          % if more channels => copy = their weighted sum!
    'weights', [1, 0.5], ...                    % weights for the sum of "copy" signal
    'shiftLag', [0, -0.4], ...                  % in [s], (circ)shifts the signal in time by shiftLag (if shiftLag > 0 => shift into the future)
    'withPrep', 'noPrep', ...                   % add preparatory activity?, choices: onlyPrep, mixPrep, noPrep
    'scaleFactor', 10 ...                       % in [std], inverse SNR (noise/copy)
    );
params.model_activPattern = struct(...  % model = noise + activPattern. Typically: copy = xVel,absVel, noise = whiteNoise (~ "copy+noise" model)
    'scaleFactor', 15 ...                   % in [std], inverse SNR (noise/copy)
    );
params.dipole2ecog = struct(...         % uses the potential of the forward models & velocity to create ecog datasets
    'voltageFile', '/export/jiritmp/webdavmirror/WhiteNoiseECoG/WN1/dipole_forwardModels/potentials_5mm.mat', ...
    'signalType', 'Traj', ...
    'chName', 'xVel'  ...    
    );
params.broadband = struct(...           % short time Fourier transform, avg over freq. band
    'windowSize', 1.0, ...                  % in [s]
    'timeStep',   0.5, ...                  % in [s] 
    'pscPlot', false ...
    );   
params.suaSimulation_populations = struct(...           % simulated SUA firing rate model (different populations sizes)
    'popSizes', [1e0,1e1,1e2,1e3,1e4,1e5], ...         % population sizes
    'modelType', {{'dist','pos','spd','vel','noise'}} ... % covariates of discharge model
    ); 
params.suaSimulation_repetitions = struct(...           % simulated SUA firing rate model (different repetitions, PD distributions)
    'popSize', 1000, ...         % population size
    'nRepetitions', 100, ...     % number of repetitions ("channels" with different PDs & noise)
    'modelType', {{'dist','pos','spd','vel','noise'}} ... % covariates of discharge model
    );   
params.dist_car_road = struct(...                   % computes distance between car and the road (cont. error signal)
    );

params.plot_filtering = struct(...
    'rawChannels',  false, ...                          % is this used anywhere?
    'filterDesign', false  ...                          % plot low/hi band pass filter freqz 
    );

%% data processing: filter lists for predictor & response
params.predictor = struct(...                   % predictor variable info (used in regression, tuning, triggering, ..)
    'D_field', 'ampData', ...                                % MUST match a field name in D-struct (amplifier data, tracker, ...)
    'signalType', 'iEEG', ...                                % selection of signal type
    'rejection', {{''}}, ...                                 % rejection of data?
    'dataProcessing', {{'car','bp_hiPass', 'norm2std'}}, ... % signal processing
    'saveName', 'pred' ...
    );
params.response = struct(...                    % response variable info (used in regression, tuning, triggering, ..)
    'D_field', 'tracker', ...                                % MUST match a field name in D-struct (amplifier data, tracker, ...)
    'signalType', 'xVel', ...                                % selection of signal type
    'rejection', {{''}}, ...                                 % rejection of data?
    'dataProcessing', {{'bp_loPass','norm2std'}}, ...        % signal processing
    'saveName', 'resp' ...
    );

%% decoding: multiple linear regression
params.regression = struct(...
    'doLinRegression', false, ...
    'timeLags', [-3.5:0.1:3.5], ... % in [s], distance of predictor from response samples, predictive <=> timeLag < 0 !!!
    'extractionStep', 0.02, ...     % in [s], step length in extracting the features  
    'sessionOffset', 3.6, ...       % distance in [s] for session beginning / end (ignores this period)
    'methodName', 'mlr_basic', ...  % options: 'mlr', 'ridge', 'mlr_stats'
    'nSecPerFold', 30  ...          % number of seconds per validation fold
    );

params.channels4decoding = struct(...
    'select', 'singleChnl' ...      % options: 'allChnls', 'singleChnl' 
    );

%% rejection
% data in time: 'normOutliers','epiSpikes', 'manualObservation','movtPauses','uboat_sessionState'
% channels: 'chnls_broken', 'chnls_epi_IIA', 'chnls_epi_SOZ', 'chnls_epiSpikes'
params.rejection.general = struct(...
    'timeMargins', 0.1 ...      % in [s], time offset included into rejection (appends additional data before & after)
    );
params.rejection.epiSpikes = struct(...
    'mode', 'compute' ...      % choices: 'compute' or 'load_fromD'
    );
params.rejection.movtPauses = struct(...    % TO DO ...
    'pauseTime', 0.200, ...     % in [s]
    'speedThr',  0.050 ...      % in [std]
    );
params.rejection.normOutliers = struct(...   % rejects outliers from normal distributions (larger than 6 sigma)
    'sigmaThr', 6 ...           % in [std]
    );
params.rejection.certainValue = struct(...    % TO DO ...
    'val2Reject', 0.000, ...     % value to be rejected
    'tolerance',  0.005 ...      % +/- limits
    );
params.rejection.sessionEdges = struct(...    % TO DO ...
    'edgeTime', 2.0 ...     % in [s], from session beginning & end
    );
params.uboat_sessionState = struct();           % rejects all epochs appart from 'search'
params.rejection.uboat_trialType = struct(...   % rejects uboat game type (normal, wild movements)
    'paradigmType', 1 ...           % 0 = normal game, 1 = wild movements
    );
params.rejection.chnls_broken = struct();       % broken, not working or artifact channels, identified by Jiri
params.rejection.chnls_epi_IIA = struct();      % heavy inter-ictal activity (IIA) identified by neurologists. Stored in the header files.
params.rejection.chnls_epi_SOZ = struct();      % seizure onset zone (SOZ) identified by neurologists. Stored in the header files.
params.rejection.chnls_epi_LESION = struct();   % lesion identified by neurologists. Stored in the header files.
params.rejection.chnls_epiSpikes = struct();    % heavy IIA identified by detector (e.g. Janca et al.). Stored in the header files.
params.rejection.chnls_explicit = struct(...    % just an example (add rows for more subjects)
    'subj_chnls', {{'S1', [32, 64];             % tuple: subjTag + list of channels to reject
                    'S2', [1 2 3];
                    'S3', [7 8 9]}} ...
    );

%% tuning
params.tuning = struct(...
    'doTuning', false, ...
    'method', 'simple', ...                             % response extraction, options: 'simple','predResp_thirdDim','respOnly_thirdDim', 'freq_ch2ch', 'nothing'    
    'timeLags', [-3.5:0.1:3.5], ...
    'predictor', params.predictor, ...
    'response', params.response, ...
    'ampBinMethod', 'constSamplesPerBin', ...           % options: 'fixedBinSize', 'constSamplesPerBin', 'defBinSizes' 
    'phiBinMethod', 'constSamplesPerBin', ...           % options: 'fixedBinSize', 'constSamplesPerBin', 'defBinSizes'
    'nAmpBins', 20, ...                                 % number of bins for the real (amplitude) part
    'nPhiBins', 36, ...                                 % number of bins for the imag (anglular)  part
    'nSurrogates', 100 ...                              % number of surrogate datasets for null distributions 
    );
%params.freqBands = [10 24; 30 45; 55 200];             % used only if tuning  method set to: 'freq'
params.binSizes = struct(...                            % used only if binning method set to: 'defBinSizes'
    'ampDefinition', [], ...                            % for real & magnitude of complex
    'phiDefinition' , [] ...                            % for phase of complex
    );
params.plot_tuning = struct(...
    'hist_pred', false, ...                            % histogram of predictors
    'hist_resp', false,  ...                           % histogram of responses
    'resp_results', true ...                           % tuning results of responses
    );

%% mutual information
params.mutualInfo = struct(...
    'doMutualInfo', false, ...
    'timeLags', [-3.5:0.1:3.5], ...
    'predictor', params.predictor, ...
    'response', params.response, ...
    'binMethod', 'constSamplesPerBin', ...            % options: 'fixedBinSize', 'constSamplesPerBin', 'defBinSizes' 
    'nAmpBins', 9 ...                                 % number of bins for the real (amplitude) part
    );

%% triggering
params.triggering = struct(...
    'doTriggering', false, ...
    'response', params.response, ...        % extracted (response) signal 
    'classes', 'defined below', ...         % defines the classes, options: see defineTriggers.m
    'cutPoint', 'onGo', ...                 % choices: onGo, onButtonPress, onShowStimuli
    'method','fieldName_fieldVals', ...     % defines the trigger method to extract the trials defined in 'classes'
    'time2cut', [-2, 3], ...                % in [s], w.r.t. cutPoint (before & after trigger time)
    'doSpectra', 'stft', ...                % do spectral analysis, choices: stft, mtft, none
    'baseline', [-1.0, -0.5], ...           % in [s], w.r.t. cutPoint 
    'freqBands', [], ...                     % defined below
    'behaviour',{{'RT','RWM'}} ...          % behaviour from paradigm: RT=reaction times, RWM=right/wrong/missed, BP_RedGrn=red/green button press
    );
params.plot_triggering = struct(...
    'time2plot', [-3.5, 3.5], ...           % time to plot, in [s] w.r.t. cutting point (= 0 s)
    'verLines', [-3:3], ...                 % vertical dashed lines
    'horLines', [0], ...                    % horizontal dashed lines
    'brainTopo', false, ...                 % brain topology in 3D transparent brain
    'trials_avg', true, ...                 % trials means
    'spectra_avg', true, ...                % spectra means
    'trials_all', false, ...                % single trials: time domain data
    'spectra_all', false, ...               % single trials: spectral data
    'stackedTrials', false, ...             % stacked trials for each channel
    'paraTimes', {{'go','rt'}} ...          % paradigm times: e.g. (reaction times), choices: go, rt, nextTrial_rt  
    );

%% connectivity
params.connectivity = struct(...
    'computeConnectivity', false, ...
    'connectivityMeasure', {{'DTF'}}, ...   % choice: 'GGC','DTF','ffDTF','PDC','COH'
    'MVAR_order', 10, ...
    'subBands', [1 140], ... settings by RJ: [2 12;13 25;26 48;52 75; 76 98; 102 148];
    'windowSize', 0.5, ...      % in [s]
    'timeStep', 0.125, ...      % in [s], corresponds to N samples at 512 Hz (if N = 64, tStep = 1/(512/64) = 0.125)
    'selectedChnls', 'see below', ... % see below: params.coherence.selectedChnls
    'maxDistToROI', 5, ...      % in [mm], maximum distance (tolerance) to region of interest
    'plotEC_ch2ch', true, ...        % plot coherence: each channel to all channels
    'plotEC_selGroups', true, ...    % plot coherence: selected groups of channels (each group to each group)
    'time2plot', [-3,3] ...     % in [s]
    );
params.connectivity.selectedChnls = {...    % each row = tuple: atlasName, areaName (must match names in anatomicalAreas_getList.m) 
    'Yeo7','Default'; ...
    'Yeo7','Dorsal Attention'; ...  
    };  

%% results saving
params.saveResults = struct(...
    'makeClean', true, ...
    'exportTriggering_trials', [], ...  % [] = no export, else: path to directory  
    'clearSessionCache', true, ...
    'clearCacheFile', false, ...
    'clearOutputDir', false ...
    );

%% plotting
% brain slices (useful for single-channel plots together with activations)
params.plot_brainSlices = struct(...           % user interface structure: holds most (but not all!) of the user settings
    'plot_slices', false, ...                            % 2D brain slices
    'volumes2plot', {{...       % volumes to plot (coregistered .nii files)
        'wT1', gray(128), [1, 1];  ...     % name & colormap & transparency limits
        'wrCT', copper(128), [0 0.8]; ...
        'ch-spec_Yeo7', autumn(128), [0, 0.4]; ...
        }}, ...             
    'slicePlanes', {{'axial','sagittal','coronal'}},... % planes to plot, options: 'axial','sagittal','coronal'
    'size_interpolate', 1.0, ...                        % in [mm], voxel size to which the brain is interpolated
    'printFormats', {{'fig', 'png'}}, ...
    'printResolution', 0 ...                           % choices: 0 (= screen resolution) or 600 (= dpi). Resolution of the figures.                          
); 

% 3D brain model (useful for highlighting selected networks and many channels from multiple subjects)
params.plot_brain3D = struct(...        % user interface structure: holds most (but not all!) of the user settings, see also inside brain3D_plot.m
    'plot_brain3D', false, ...           % 3D brain model
    'volumes2plot', {{...               % volumes to plot (coregistered .nii files)
        'wc1T1_colin27', [1 1 1], 0.1;  ...     % 1x3 cell: name, color, transparency
        'wDefault', [0 1 0], 0.2; ...
        'wDorsal Attention', [0 0 1], 0.2; ...
        }}, ...          
    'chVals_asCircles', true, ...       % plots channel values as color circles
    'chVals_asPatches', true, ...       % each channel value is mapped onto surface (but may be overwritten by latter points!), see also distanceThreshold for focality inside brain3D_plot.m
    'chVals_colorMap', jet(256), ...    % colormap of channel values
    'chVals_patchTransp', 0.2, ...      % transparency of channel values
    'circleSizeLims', [5 60], ...       % in points, maps chnl values to circle size
    'size_interpolate', 1.0, ...        % in [mm], voxel size to which the brain is interpolated
    'backgroundColor', 'k', ...         % figure / axes background color, good choices: 'w', 'k'
    'text_LR', true, ...                % marks L / R side
    'text_tag', 'P1_forExample', ...    % string
    'viewAngle', [0, 90], ...           % axial = 0,90; sagittal = 100,10; coronal front = 0,180; coronal back = 0,0
    'visible_axis', false ...           % if axis & their labels are visible
); 



%% old (kept for compatibility with older scripts): TO DO - remove at some point!
params.plot = struct(...    
    'scrnsz_pixels', get(0, 'Screensize'), ...          % screensize in [pixels]
    'scrnsz_inches', [0, 0, 21, 13], ...                % screensize in [inches], use: set(0,'units','inches'); get(0, 'Screensize')
    'rawChannels', false, ...
    'filterDesign', false, ...                          % plot low/hi band pass filter freqz 
    'rejection', false, ...
    'errorDistr', false, ...
    'responseDistr', false, ...
    'predictedResponse', false, ...
    'scatterFeaturesFit', false, ...
    'predVsResp', false, ...
    'brainMaps', false ...
    );

params.plot_general = struct(...
    'export_resolution', 600 ...                        % in [dpi], figures are exported in this resolution as png (or tif)
    );
params.plot_brainTopo = struct(...           % user interface structure: holds most (but not all!) of the user settings
    'time2plot', [-0.1:0.1:1.0], ...                    % in [s], w.r.t. 0 = cutting point (for example)
    'plot_slices', true, ...                            % 2D brain slices
    'plot_projections', true, ...                       % 3D brain model projections (axial, sagittal, coronal)
    'plot_animation', false, ...                        % 3D brain model GIF animation (takes longer time)
    'printResolution', 0, ...                           % choices: 0 (= screen resolution) or 600 (= dpi). Resolution of the figures.                          
    'colorMap', jet(128), ...                           % colormap for channel values
    'size_interpolate', 1.0, ...                        % in [mm], voxel size to which the brain is interpolated
    'size_coloredCube', 3.0, ...                        % in [mm], "voxel" size of the colored channel values
    'slicePlanes', {{'axial','sagittal','coronal'}},... % planes to plot, options: 'axial','sagittal','coronal'
    'model_views', [0 0; 0 89.999; 90 0]  ...           % camera view angles at which the snapshot of the 3D model are taken. For some odd reason (camera light?), view(0,90) makes grey background...
);

