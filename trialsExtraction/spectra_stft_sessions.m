function spectra = spectra_stft_sessions(params, trials)
% computes spectra for each trial using STFT
% trials = 
%           data: [3073x125x500 double]
%           time: [3073x1 double]
%       sessTime: [3073x500 double]
%         labels: [1x500 double]
%     sessNumber: [1x500 double]

% (c) Jiri, May16

disp('Computing spectra from sessions ...');

%% load sampling rate
load(params.storage.cacheFile, 'srate');
fs = srate;

%% STFT settings
windowSize = params.stft_freq.windowSize;
timeStep = params.stft_freq.timeStep;
nfft = 2^nextpow2(windowSize * fs);            % number of FFT points, windowSize in [samples]
timeStep = ceil(timeStep * fs);                % in [samples], if any downsampling occurs, 'timeStep' MUST exist
freqAxis = fs/2 * linspace(0,1,nfft/2+1);      % in [Hz], up to Nyquist freq.
i_fr = closestval(freqAxis,params.stft_freq.freqBand(1)):closestval(freqAxis,params.stft_freq.freqBand(2)); % indices of selected frequencies
srate_stft = fs/timeStep;                      % in [Hz] !!! adjust srate -> step-wise processing 

%% allocation: 4D = [freq x time x chnl x trials]
% i_cut = ceil(params.triggering.time2cut(1)*srate_stft):ceil(params.triggering.time2cut(2)*srate_stft);     % in [samples], w.r.t. cutting point = 0   
i_cut = round(params.triggering.time2cut(1)*srate_stft):round(params.triggering.time2cut(2)*srate_stft);     % in [samples], w.r.t. cutting point = 0   
t_spectra = downsample(trials.time, timeStep);
assert(length(i_cut) == size(t_spectra,1));
data = nan(length(i_fr), length(i_cut), size(trials.data,2), size(trials.data,3));
n_tr = 1;
labels = [];
rejected = nan(length(i_cut), size(trials.data,2), size(trials.data,3)); % 3D: samples x channels x trials

%% processing STFT (short-time Fourier transform)
for sess = 1:size(params.storage.sessionCacheFiles,2)
    
    % load triggers
    clear triggers classNamesLabels;
    load(params.storage.sessionCacheFiles{sess}, 'triggers', 'classNamesLabels');
        
    % load responses
    clear resp_prc;
    load(params.storage.sessionCacheFiles{sess}, 'resp_prc');
    rawData = resp_prc;
    
    % load session time (from D struct)
    clear timeAxis;
    load(params.storage.sessionCacheFiles{sess}, 'timeAxis'); 
    assert(exist('timeAxis','var') == 1);
    
    % downsample time axis? (should not be needed)
    if srate ~= params.init_srate                                            % to dowsample & reject bad recordings epochs
        disp(['WARNING: downsampling the time axis of session: ' num2str(sess) ' ...']);
        params.downsample.dsRate = 'fromCacheFile';
        params.thisSess = sess;
        timeAxis = filterData(params, timeAxis, {'downsample'});
    end        
    timeAxis = downsample(timeAxis, timeStep);  % adjust to STFT

    % rejected indices
    clear resp_rej;
    load(params.storage.sessionCacheFiles{sess}, 'resp_rej');
    if ~exist('resp_rej','var')
        warning('No rejection performed.');
        resp_rej = false(size(resp_prc));
    end
    assert(size(resp_prc,1) == size(resp_rej,1));  
    assert(size(resp_prc,2) == size(resp_rej,2));  
    resp_rej = downsample(resp_rej, timeStep);  % adjust to STFT
    
    % spectral whitening: normalization factors
    B_avg = [];     % 2D: freq x ch
    B_std = [];     % 2D: freq x ch
    
    % >>> STFT <<<
    for ch = 1:size(rawData,2)
        x = cat(1, flipdim(rawData(1:floor(nfft/2)-1,ch),1), rawData(:,ch), flipdim(rawData(end-floor(nfft/2)+1:end,ch),1));
        [~,~,~,P] = spectrogram(x, hann(nfft), nfft-timeStep, nfft, srate);
        assert(size(P,1) == length(freqAxis));
        assert(size(P,2) == length(timeAxis));
        amps = 10*log10(P(i_fr,:));

        % allocation of stftData
        if ch == 1
            stftData = nan(size(amps,1),size(amps,2),size(rawData,2));  % 3D: freq x samples x chnls
        end
        
        % baseline definition
        i_base = stftBaselineIndices(params, size(amps,2));
        assert(size(i_base,2) == size(resp_rej,1));
        i_base(resp_rej(:,ch)) = [];
        
        % spectral whitening: normalization to baseline -> B_avg & B_std
        %normFactors = median(amps(:,i_base), 2);
        B_avg_ch = mean(amps(:,i_base), 2);
        B_std_ch = std(amps(:,i_base), 0, 2);
        %amps = amps./repmat(normFactors, 1, size(amps,2));
        amps = amps - repmat(B_avg_ch, 1, size(amps,2));
        B_avg = cat(2, B_avg, B_avg_ch);
        B_std = cat(2, B_std, B_std_ch);
        
        % z-score
        amps = amps ./ repmat(B_std_ch, 1, size(amps,2));
        
        % cat to stft data
        stftData(:,:,ch) = amps;
    end
    assert(size(stftData,2) == size(resp_rej,1));
    
    % plot spectral whitening normalization factors
    plot_stftBaseline(params, freqAxis(i_fr)', timeAxis, B_avg, B_std, ['baseline_RS_sess' num2str(sess)]);
        
    % extract trials data -> 4D: freq x samples x channels x trials
    for tr = 1:size(triggers,2)
        cutPoint = triggers(1,tr);                          % in [s]
        i_cutPoint = closestval(timeAxis, cutPoint);        % in [samples]
        inds = i_cutPoint + i_cut;                          % samples to extract
        if inds(1) > 0 && inds(end) <= size(timeAxis,1)
            data(:,:,:,n_tr) = stftData(:,inds,:);          % 4D: freq x samples x channels x trials
            rejected(:,:,n_tr) = resp_rej(inds,:);  % 3D: samples x channels x trials
            labels = cat(2, labels, triggers(2,tr));        % class labels
            n_tr = n_tr+1;
        else
            disp(['WARNING: skipping trial (out of session bounds): session = ' num2str(sess) ', trial = ' num2str(tr)]);
        end
    end
    
    disp([' - session ' num2str(sess) ' done.']);
end

%% log-normal distribution (in dB = 10*log10)
% data = 10*log10(data);

%% rejected indices
assert(size(rejected,1) == size(data,2));
rrr(1,:,:,:) = rejected;
rejected = repmat(rrr, [size(i_fr,2),1,1,1]);   % same for all freq
assert(size(rejected,1) == size(data,1));
assert(size(rejected,2) == size(data,2));
assert(size(rejected,3) == size(data,3));
assert(size(rejected,4) == size(data,4));

%% spectral whitening (baseline normalization)
% selTime = params.triggering.baseline;
% i_t = closestval(timeAxis,selTime(1)):closestval(timeAxis,selTime(2));
% for ch = 1:size(trials.data,2)
%     
%     % cat baseline values over time periods
%     bs_vals = [];               % baseline values
%     rj_inds = [];               % rejected indices
%     for tr = 1:size(trials.data,3)
%         bs_vals = cat(2, bs_vals, data(:,i_t,ch,tr));
%         rj_inds = cat(2, rj_inds, rejected(:,i_t,ch,tr));
%     end
%     
%     % rejection
%     inds = ~rj_inds;            % indices to use (= not rejected indices)
%     vals = bs_vals;             % all baseline values
%     vals(inds == 0) = NaN;      % those samples that are rejected are set to NaN   
%     
%     % baseline definition
%     %base = nanmedian(vals,2);    % ~ over time
%     base = nanmean(vals,2);    % ~ over time
%     assert(any(~isnan(base(1))));
%     
%     % baseline normalization
%     %data(:,:,ch,:) = data(:,:,ch,:)./repmat(base, [1,size(data,2),1,size(data,4)]);
%     data(:,:,ch,:) = data(:,:,ch,:) - repmat(base, [1,size(data,2),1,size(data,4)]);
%     
%     if mod(ch,10) == 0
%         display([' - whitening: channel ' num2str(ch) ' of ' num2str(size(trials.data,2)) ' - done.']);
%     end        
% end

%% log-normal distribution (in dB = 10*log10) (as of 03.07.2018)
%data = 10*log10(data);

%% output
spectra = struct;
spectra.data = data;
spectra.freq = freqAxis(i_fr);
spectra.time = t_spectra;
spectra.labels = labels;
spectra.clzNames = trials.clzNames;
spectra.rejected = rejected;
