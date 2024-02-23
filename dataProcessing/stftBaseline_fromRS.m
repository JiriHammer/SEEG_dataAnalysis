function [B_avg, B_std] = stftBaseline_fromRS(params)
% determines baseline from resting state (RS)
%   - baseline for normalization (normalization factors) -> B_avg (& B_std)
%   - baseline activity -> BA_avg & BA_std (saves 'baseData' struct to cache file)
%
% (1) find RS sessions (based on gameType)
% (2) compute STFT of the RS
% (3) rejection of bad data segments
% (4) define baseline mean & SD -> B_avg & B_std = 3D: freq x 1 x ch
% (5a) whiten (normalize) the RS spectra by B_avg -> B_val
% (5b) z-score of the RS spectra B_val by B_std -> B_val
% (6) compute mean & SD over time of the normalized RS spectra -> BA_avg & BA_std = 3D: freq x 1 x ch
% (7) plot baseline for all channels (B_avg)

% (c) Jiri, Apr22, Sep23

%% (1) find RS sessions (based on gameType) -> i_sess
i_sess = [];
for sess = 1:size(params.storage.sessionCacheFiles,2)
    clear gameType
    load(params.storage.sessionCacheFiles{sess}, 'gameType');
    if strcmp(gameType, 'restingState')
        i_sess = [i_sess, sess];
    end
end

% if no RS found, return empty vals
if isempty(i_sess)
    B_avg = [];
    B_std = [];
    return;
end

%% (2) compute STFT of the RS -> B_val & I_rej
% load srate from cache file
clear srate;
load(params.storage.cacheFile, 'srate');   % load from cache file!
fs = srate;

% STFT settings
nfft = 2^nextpow2(params.stft_freq.windowSize * fs);            % number of FFT points, windowSize in [samples]
timeStep = ceil(params.stft_freq.timeStep * fs);                % in [samples], if any downsampling occurs, 'timeStep' MUST exist
freqAxis = fs/2 * linspace(0,1,nfft/2+1);                     % in [Hz], up to Nyquist freq.
i_fr = closestval(freqAxis,params.stft_freq.freqBand(1)):closestval(freqAxis,params.stft_freq.freqBand(2)); % indices of selected frequencies
t_cutEdges = 2;     % in [s], removes the edges from the RS

% STFT of RS
B_val = [];
R_val = [];
for sess = 1:size(i_sess,2)
    clear resp_prc resp_rej time_prc srate;
    load(params.storage.sessionCacheFiles{sess}, 'resp_prc', 'resp_rej', 'time_prc', 'srate');
    assert(fs == srate);    % srate in sessionCache should be the same as in cache file
    rawData = resp_prc;     % 2D: samples x channels
    T = downsample(time_prc, timeStep); % downsampled time axis
    
    % processing STFT, single freq. band, PSD
    P_vals = nan(size(freqAxis(i_fr),2), size(T,1), size(rawData,2));  % 3D: freq x time x chnls
    for ch = 1:size(rawData,2)
        x = cat(1, flipdim(rawData(1:floor(nfft/2)-1,ch),1), rawData(:,ch), flipdim(rawData(end-floor(nfft/2)+1:end,ch),1));
        [~,~,~,P] = spectrogram(x, hann(nfft), nfft-timeStep, nfft, fs);
        assert(size(P,1) == size(freqAxis,2));
        assert(size(P,2) == size(T,1));
        P_vals(:,:,ch) = P(i_fr,:);
    end

    % downsample rejection inds
    R = downsample(resp_rej, timeStep);
    assert(size(R,1) == size(T,1));

    % reject edges of baseline
    i_t = 1:closestval(T, T(1)+t_cutEdges);                 % = first 2 s
    R(i_t,:) = true;
    i_t = closestval(T, T(end)-t_cutEdges):size(T,1);       % = last 2 s
    R(i_t,:) = true;

    % cat sessions over time
    B_val = cat(2, B_val, P_vals);      % 3D: freq x time x chnls
    R_val = cat(1, R_val, R);           % 2D: time x chnls  
end
assert(size(B_val,2) == size(R_val,1));     % time
assert(size(B_val,3) == size(R_val,2));     % chnls

% reshape rejected indices to 3D: freq x time x chnls
RRR(1,:,:) = R_val;
I_rej = repmat(RRR, [size(B_val,1),1,1]);
assert(size(B_val,1) == size(I_rej,1));     % freq
assert(size(B_val,2) == size(I_rej,2));     % time
assert(size(B_val,3) == size(I_rej,3));     % chnls

%% log-transform
B_val = 10*log10(B_val);           % 3D: freq x time x chnls          % log of PSD
    
%% (3) rejection of bad data segments -> I_rej
% rejection criteria: outliers of normal distribution per freq. bands (FB)
sigmaThr = params.rejection.normOutliers.sigmaThr;
for fb = 1:size(B_val,1)
    for ch = 1:size(B_val,3)
        d = squeeze(B_val(fb,:,ch));    % vector of samples
        [muhat,sigmahat] = normfit(d);
        i_rej = abs(d-muhat) > sigmaThr*sigmahat;     % = samples, whose distance from distribution mean is greater than N sigma
        I_rej(fb,i_rej,ch) = true;
    end
end

% reject
I_rej = logical(I_rej);
B_val(I_rej) = nan;             % where I_rej=true, B_val=nan

%% (4) define baseline mean & SD -> B_avg & B_std
B_avg = nanmean(B_val,2);       % mean over time (with rejection)
B_std = nanstd(B_val, 0, 2);    % SD over time (with rejection)
if any(isnan(B_avg))
    disp(' - baseline from RS: NaN values encountered, please check!');
    disp('   -- returning empty values');
    B_avg = [];
    B_std = [];
end

%% (5a) whiten (normalize) the RS spectra B_val by B_avg -> B_val
assert(size(B_val,1) == size(B_avg,1));
B_val = B_val - repmat(B_avg(:,1,:), [1,size(B_val,2),1]);   % minus because data were log transformed!

%% (5b) z-score of the RS spectra B_val by B_std -> B_val
assert(size(B_val,1) == size(B_std,1));
B_val = B_val ./ repmat(B_std(:,1,:), [1,size(B_val,2),1]);

%% (6) compute mean & SD over time of the normalized RS spectra -> BA_avg & BA_std
BA_avg = nanmean(B_val,2);      % mean over time (with rejection)
BA_std = nanstd(B_val, 0, 2);    % SD over time (with rejection)

%% (7) plot baseline for all channels (B_avg)
plot_stftBaseline(params, freqAxis(i_fr)', T, squeeze(B_avg), squeeze(B_std));

%% save baseData struct (BA_avg & BA_std)
baseData = struct;
baseData.xVals = freqAxis(i_fr)';
baseData.yVals = squeeze(BA_avg);
baseData.yErrs = squeeze(BA_std);
params.nCh = size(baseData.yVals,2);
baseData.info.chNames = getChannelNames(params, lower(params.response.signalType), 'resp');
baseData.info.xlabel = 'freq [Hz]';
baseData.info.verLines = [10, 30, 50, 100];
baseData.info.text = ['subj = ' params.storage.subjTag ', baseline from resting state, duration = ' num2str(T(end)/60) ' min'];
baseData.info.outDir = [params.storage.outputDir filesep 'baseline_RS'];
baseData.info.figName = 'baseline_RS';
save(params.storage.cacheFile, 'baseData', '-append');

