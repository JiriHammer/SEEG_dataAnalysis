function spectra = spectra_stft_trials(params, trials)
% computes spectra for each trial using STFT
% trials = 
%           data: [3073x125x500 double]
%           time: [3073x1 double]
%       sessTime: [3073x500 double]
%         labels: [1x500 double]
%     sessNumber: [1x500 double]

% (c) Jiri, May16

display('Computing spectra from trials ...');

%% load sampling rate
load(params.storage.cacheFile, 'srate');
fs = srate;

%% STFT settings
windowSize = params.stft_freq.windowSize;
timeStep_sec = params.stft_freq.timeStep;
nfft = 2^nextpow2(windowSize * fs);            % number of FFT points, windowSize in [samples]
timeStep = ceil(timeStep_sec * fs);                % in [samples], if any downsampling occurs, 'timeStep' MUST exist
freqAxis = fs/2 * linspace(0,1,nfft/2+1);      % in [Hz], up to Nyquist freq.
t_full = downsample(trials.time, timeStep);
i_t = closestval(t_full, t_full(1)+windowSize/2):closestval(t_full, t_full(end)-windowSize/2);  % cut-out half of window size from each end
timeAxis = t_full(i_t);

%% selected frequency
selFreq = [0 150];      % in [Hz]
i_freq = closestval(freqAxis, selFreq(1)):closestval(freqAxis, selFreq(2));
freqAxis_sel = freqAxis(i_freq);

%% allocation: 4D = [freq x time x chnl x trials]
data = nan(length(i_freq), length(timeAxis), size(trials.data,2), size(trials.data,3));

%% processing STFT (short-time Fourier transform)
if strcmp(params.triggering.doSpectra,'stft') || strcmp(params.triggering.doSpectra,'stft_trials')
    for ch = 1:size(trials.data,2)
        for tr = 1:size(trials.data,3)
            x = cat(1, flipdim(trials.data(1:floor(nfft/2)-1,ch,tr),1), trials.data(:,ch,tr), flipdim(trials.data(end-floor(nfft/2)+1:end,ch,tr),1));
            [~,~,~,P] = spectrogram(x, hann(nfft), nfft-timeStep, nfft, fs);
            assert(size(P,1) == length(freqAxis));
            assert(size(P,2) == length(t_full));
            data(:,:,ch,tr) = P(i_freq,i_t);
        end
        if mod(ch,10) == 0
            display([' - STFT: channel ' num2str(ch) ' of ' num2str(size(trials.data,2)) ' - done.']);
        end          
    end
end

%% processing MTFT (multi-taper Fourier transform)
if strcmp(params.triggering.doSpectra,'mtft_trials')
    nw = 2;     % number of Slepian tapers (higher number makes worse freq. resolution, but smoother spectra)
    for ch = 1:size(trials.data,2)
        for tr = 1:size(trials.data,3)
            x = cat(1, flipdim(trials.data(1:floor(nfft/2)-1,ch,tr),1), trials.data(:,ch,tr), flipdim(trials.data(end-floor(nfft/2)+1:end,ch,tr),1));
            i_window = 1:ceil(windowSize*fs);
            P = nan(length(freqAxis), length(t_full));
            n = 1;
            while i_window(end) <= size(x,1)
                y = x(i_window);
                [py,f] = pmtm(y,nw,length(y),fs);   % multi-taper method
                P(:,n) = py;
                i_window = i_window + timeStep;
                n = n+1;
            end
            %[~,~,~,P_s] = spectrogram(x, hann(nfft), nfft-timeStep, nfft, fs);
            assert(size(P,1) == length(freqAxis));
            assert(size(P,2) == length(t_full));
            data(:,:,ch,tr) = P(i_freq,i_t);
        end
        disp([' - MTFT: channel ' num2str(ch) ' of ' num2str(size(trials.data,2)) ' - done.']);
    end
end

%% processing PSD (STFT PSD implemented by Radek Janca) - !BUG: time axis is wrong!
if strcmp(params.triggering.doSpectra,'psd_trials')
%     timeStep = 0.03125;   % in [s], corresponds to integer N samples at 512 Hz (if N = 8, tStep = 1/(512/8) = 0.015625)
%     windSize = 0.5;   % in [s]
    windSize = params.stft_freq.windowSize;     % in [s]
    timeStep_sec = params.stft_freq.timeStep;       % in [s]
    setting.ww=windSize; % segmentation window (sec.)
    setting.nn=(windSize-timeStep_sec)/windSize; % negative overlap x100%

    % === PSD from Radek Janca toolbox=== 
    [Px, freqAxis, timeAxis_bug] = PSD_v1(trials.data, trials.time, trials.labels, setting);
    data = Px(i_freq,:,:,:);
    clear Px;

    % time step (re)calculation
    timeStep = ceil(setting.ww * (1-setting.nn) * fs); % time step in [samples]
    t_full = downsample(trials.time, timeStep);
    i_t = closestval(t_full, t_full(1)+windSize/2):closestval(t_full, t_full(end)-windSize/2);  % cut-out half of window size from each end
    timeAxis = t_full(i_t);
    assert(size(data,2) == size(timeAxis,1));
%     timeAxis = timeAxis_bug;
    freqAxis_sel = freqAxis(i_freq)';   % assumes a row vector 1 x N 
end

%% processing CWT (continuous wavelet transform)
if strcmp(params.triggering.doSpectra,'cwt_trials')
    setting.ww=0.5; % segmentation window (sec.)
    setting.nn=0.9; % negative overlap x100%
    
    % === CWT ===
    setting.w_freg_lim = selFreq;
    [data, FW, TW] = CWT_v1(trials.data, trials.time, trials.labels, setting);  % CWTx = 4D: F x T x Ch x Tr
    data = flip(data,1);            % reverse the order along freq. dim
    freqAxis_sel = flip(FW,1)';     % assumes a row vector 1 x N
    assert(freqAxis_sel(1) < freqAxis_sel(end));

    % downsample to timeStep
    data = permute(data,[2 1 3 4]);     % 4D = T x F x Ch x Tr
    data = downsample(data, timeStep);  % 4D = T x F x Ch x Tr  (or decimate???)
    data = permute(data,[2 1 3 4]);     % 4D = F x T x Ch x Tr
    t_full = downsample(TW', timeStep);
    assert(size(data,2) == size(t_full,1));
    data = data(:,i_t,:,:);             % exclude window edges
    
    % no downsampling (runs out of memory...)
%     data = CWTx;
%     timeAxis = TW;
%     timeStep = 1;       % in [samples]
%     i_t = 1:size(data,2);
%     clear CWTx;
end

%% log-normal distribution (in dB = 10*log10)
data = 10*log10(data);

%% rejected indices
r = downsample(trials.rejected, timeStep);
r = r(i_t,:,:);
assert(size(r,1) == size(data,2));
rrr(1,:,:,:) = r;
rejected = repmat(rrr, [size(freqAxis_sel,2),1,1,1]);
assert(size(rejected,1) == size(data,1));
assert(size(rejected,2) == size(data,2));
assert(size(rejected,3) == size(data,3));
assert(size(rejected,4) == size(data,4));

%% spectral whitening (baseline normalization)
selTime = params.triggering.baseline;
i_t = closestval(timeAxis,selTime(1)):closestval(timeAxis,selTime(2));
B_avg = nan(size(data,1),1,size(data,3));   % 3D: freq x 1 x ch
B_std = nan(size(data,1),1,size(data,3));   % 3D: freq x 1 x ch
for ch = 1:size(trials.data,2)
    
    % cat baseline values over time periods
    bs_vals = [];               % baseline values
    rj_inds = [];               % rejected indices
    for tr = 1:size(trials.data,3)
        bs_vals = cat(2, bs_vals, data(:,i_t,ch,tr));
        rj_inds = cat(2, rj_inds, rejected(:,i_t,ch,tr));
    end
    
    % rejection
    inds = ~rj_inds;            % indices to use (= not rejected indices)
    vals = bs_vals;             % all baseline values
    vals(inds == 0) = NaN;      % those samples that are rejected are set to NaN   
    
    % baseline definition -> B_avg, B_std
    %base = nanmedian(vals,2);    % ~ over time
    B_avg_ch = nanmean(vals,2);    % ~ over time
    B_std_ch = nanstd(vals,0,2);    % ~ over time
    assert(any(~isnan(B_avg_ch(1))));
    assert(any(~isnan(B_std_ch(1))));
    
    % baseline normalization -> z-score
    %data(:,:,ch,:) = data(:,:,ch,:)./repmat(base, [1,size(data,2),1,size(data,4)]);
    data(:,:,ch,:) = data(:,:,ch,:) -  repmat(B_avg_ch, [1,size(data,2),1,size(data,4)]);
    data(:,:,ch,:) = data(:,:,ch,:) ./ repmat(B_std_ch, [1,size(data,2),1,size(data,4)]);
    
    if mod(ch,10) == 0
        display([' - whitening: channel ' num2str(ch) ' of ' num2str(size(trials.data,2)) ' - done.']);
    end       
    B_avg(:,1,ch) = B_avg_ch;
    B_std(:,1,ch) = B_std_ch;
end
save(params.storage.cacheFile, 'B_avg', 'B_std', '-append');

%% output
spectra = struct;
spectra.data = data;
spectra.freq = freqAxis_sel;
spectra.time = timeAxis;
spectra.labels = trials.labels;
spectra.clzNames = trials.clzNames;
spectra.rejected = rejected;
