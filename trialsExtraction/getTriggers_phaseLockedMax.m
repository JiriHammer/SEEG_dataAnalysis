function triggers = getTriggers_phaseLockedMax(params)
% TO DO...
% predicts signal (xVel) by MLR
% computes single-trial phase locking value (SPLV) of predicted and response signal
% triggeres on maxima of SPLV

% (c) Jiri, May12

% TO DO, is it working? (now: Mar13)

% set up MLR prediction
params.regressMethod.doLinRegression = true;
params.regressMethod.xvMode = 'sessionSplit';
params.regressMethod.nCV = length(params.simulation.trainSession);
params.exploreVals.doExploration = false;
params.template = struct(...
    'nLags', 0, ...            % number of past points to form a feature (currLag = 0)
    'lagDist', 0.10, ...       % distance in [s] between neighbouring lags
    'templStep', 0.02, ...       % distance in [s] of the template extracting the features  
    'zeroLagOffset', 0.00, ...   % distance in [s] between current point of response and 1st feature (zero lag)
    'beginOffset', 1e-5, ...      % distance in [s] for beg. of extraction point of response and 1st  feature (-zero lag)
    'endOffset', 1e-5, ...        % distance in [s] for end  of extraction point of response and last feature
    'sameLagsPerCh', true ...
    );
params.stft_freqComp_pred = struct(...   % short time Fourier transform -> complex frequency components
    'windowSize', 2.0, ...              % in [s]
    'timeStep', 0.02, ...               % in [s]
    'freqBins', 2, ...                % obligatory field name. bin indices of STFT, corr. to freq = 1/windowSize*freqBins
    'keepDim', 'complex', ...           % options: 'complex', 'amp', 'phase'
    'ccSingleChAtLag', 0 ...            % lag in [s] (positive<=>prediction) for which to compute singleCh CC
    );
params.stft_freqComp_resp = struct(...   % short time Fourier transform -> complex frequency components
    'windowSize', 2.0, ...              % in [s]
    'timeStep', 0.02, ...               % in [s]
    'freqBins', 2, ...                % obligatory field name. bin indices of STFT, corr. to freq = 1/windowSize*freqBins
    'keepDim', 'complex', ...           % options: 'complex', 'amp', 'phase'
    'ccSingleChAtLag', 0 ...            % lag in [s] (positive<=>prediction) for which to compute singleCh CC
    );    
params.predictor = struct(...
    'signalType', 'EcogGrid', ...
    'selectedChnls', {{'all'}}, ... % 'all','prefrontal','broca','premotor','primsensory','ipcPF','ipcPFm'
    'excludedChnls', [], ...
    'dataProcessing', {{'car','bp_hiPass', 'norm2std', 'stft_freqComp_pred'}}, ...  % preprocessing of predictor signal
    'saveProcessing', true ...
    );
params.response = struct(...
    'signalType', 'Traj', ...
    'selectedChnls', {{'xVel'}}, ...
    'excludedChnls', [], ...
    'dataProcessing', {{'norm2std', 'stft_freqComp_resp'}}, ...  % preprocessing of response signal
    'saveProcessing', true ...
    );
params.saveResults.addFileToCacheBase = 'yes_all';

% run MLR prediction
processingLoop(params);

% load predicted and response Trajs
fileName = [params.storage.cacheBaseDir filesep params.response.selectedChnls{1} 'Traj_allEcogGrid/fileToCacheBase_1.mat'];
assert(exist(fileName,'file') == 2);
clear models;
load(fileName, 'models');
assert(exist('models','var') == 1);

% cumulate predictions & responses, folds <=> sessions
predSig = [];
respSig = [];
for fold = 1:size(models{1}.trajec, 2)
    predSig = cat(1, predSig, models{1}.trajec{1,fold}.origPred);  % predictions, complex values
    respSig = cat(1, respSig, models{1}.trajec{1,fold}.origResp);  % responses, complex values
end

% cumulate reconstructed trajs, folds <=> sessions
predTraj = [];
respTraj = [];
for fold = 1:size(models{1}.trajec, 2)
    predTraj = cat(1, predTraj, models{1}.trajec{1,fold}.predTraj);  % predictions, complex values
    respTraj = cat(1, respTraj, models{1}.trajec{1,fold}.respTraj);  % responses, complex values
end

% single-trial phase locking value (SPLV)
load(params.storage.cacheFile, 'sampleRateAfterFiltering'); assert(exist('sampleRateAfterFiltering','var') == 1);
srate = sampleRateAfterFiltering;
wSize = round(6.0*srate);  % in [samples]
splv_signal = computeSPLV(predSig, respSig, wSize);

% smooth SPLV
Wn = 0.1/(srate/2);                          % normalized bandpass frequencies
n = 3;                                      % butterworth order
[b,a] = butter(n, Wn, 'low');               % returns polynoms of Butterw. filter
splv_smooth = filtfilt(b,a, splv_signal);

% TO DO: threshold as stat. signif. p-value over surrogate data distr.
thrs = prctile(splv_smooth, 75);  % hack for the time being...

% plot trajs and their SPLV
figure;
hold on;
plot([1:size(predTraj)]./srate, predTraj./max(predTraj), 'c');
plot([1:size(respTraj)]./srate, respTraj./max(respTraj), 'm');
plot([1:size(predTraj)]./srate, splv_signal, 'r');
plot([1:size(predTraj)]./srate, splv_smooth, 'b');
plot(get(gca, 'xlim'), [thrs thrs], '--k');

figure; 
hist(splv_smooth, 1000);

% triggers as local maxima of splv above thrs
sampleSessVals = [];
for sess = 1:size(models{1}.trajec, 2)
    sessSizes(sess) = size(models{1}.trajec{1,sess}.origPred,1);
    splvSess = computeSPLV(models{1}.trajec{1,sess}.origPred, models{1}.trajec{1,sess}.origResp, wSize);

    sampleVals = getLocalMaxima(params, splvSess);
    i_sel = find(sampleVals(:,2) >= thrs);

%         % extract local extrema
%         if strcmp(params.locExtrema.minOrMax, 'both')
%             sampleValsMax = getLocalMaxima(params, pred{sess});
%             sampleValsMin = getLocalMaxima(params, -pred{sess});
%             sampleValsMin(:,2) = -sampleValsMin(:,2);
%             sampleVals = cat(1, sampleValsMax, sampleValsMin);
%         elseif strcmp(params.locExtrema.minOrMax, 'max')
%             sampleVals = getLocalMaxima(params, pred{sess});
%         elseif strcmp(params.locExtrema.minOrMax, 'min')
%             sampleVals = getLocalMaxima(params, -pred{sess});
%             sampleVals(:,2) = -sampleVals(:,2);
%         end
    sampleSessVals = cat(1, sampleSessVals, cat(2, sampleVals(i_sel,1), sess*ones(length(i_sel),1), sampleVals(i_sel,2)));
end

%     % bin extrema values into nBins
%     yVals = sampleSessVals(:,3);
%     tunInfo = struct;
%     tunInfo.minV = min(yVals);
%     tunInfo.maxV = max(yVals);
%     tunInfo.nBins = params.locExtrema.nBins;
%     tunInfo.binMethod = 'constSamplesPerBin';
%     tunInfo.binSizes = [];                        % empty -> bin sizes get computed (if needed)
%     [~, tunInfo.binSizes] = binVector(tunInfo, yVals);   
%     [binsArray, ~] = binVector(tunInfo, yVals);
%     assert(size(binsArray,1) == size(yVals,1));
%     binValues = sort(unique(binsArray));    

% output
triggers{1}.sampleSessLabel = cat(2, sampleSessVals(:,1), sampleSessVals(:,2), ones(size(sampleSessVals,1),1));     % mandatory field
%triggers{1}.info = info;
%triggers{1}.info.binning = tunInfo;
triggers{1}.info.sessSizes = sessSizes;                                                        % mandatory field
triggers{1}.info.dataType = 'timeSeries';                                                      % mandatory field
triggers{1}.info.chPos = params.response.selectedChnls{1};
triggers{1}.info.name = params.response.signalType;
triggers{1}.info.clzColors = colorPalette(1); 
