function spectra = computeSpectra(params, trials)
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
timeStep = params.stft_freq.timeStep;
nfft = 2^nextpow2(windowSize * fs);            % number of FFT points, windowSize in [samples]
timeStep = ceil(timeStep * fs);                % in [samples], if any downsampling occurs, 'timeStep' MUST exist
freqAxis = fs/2 * linspace(0,1,nfft/2+1);      % in [Hz], up to Nyquist freq.
t_full = downsample(trials.time, timeStep);
i_t = closestval(t_full, t_full(1)+windowSize/2):closestval(t_full, t_full(end)-windowSize/2);  % cut-out half of window size from each end
timeAxis = t_full(i_t);

%% allocation: 4D = [freq x time x chnl x trials]
data = nan(length(freqAxis), length(timeAxis), size(trials.data,2), size(trials.data,3));

%% processing STFT (short-time Fourier transform)
if strcmp(params.triggering.doSpectra,'stft')
    for ch = 1:size(trials.data,2)
        for tr = 1:size(trials.data,3)
            x = cat(1, flipdim(trials.data(1:floor(nfft/2)-1,ch,tr),1), trials.data(:,ch,tr), flipdim(trials.data(end-floor(nfft/2)+1:end,ch,tr),1));
            [~,~,~,P] = spectrogram(x, hann(nfft), nfft-timeStep, nfft, fs);
            assert(size(P,1) == length(freqAxis));
            assert(size(P,2) == length(t_full));
            data(:,:,ch,tr) = P(:,i_t);
        end
        if mod(ch,10) == 0
            display([' - STFT: channel ' num2str(ch) ' of ' num2str(size(trials.data,2)) ' - done.']);
        end          
    end
end

%% processing MTFT (multi-taper Fourier transform)
if strcmp(params.triggering.doSpectra,'mtft')
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
            data(:,:,ch,tr) = P(:,i_t);
        end
        display([' - MTFT: channel ' num2str(ch) ' of ' num2str(size(trials.data,2)) ' - done.']);
    end
end

%% log-normal distribution (in dB = 10*log10)
data = 10*log10(data);

%% rejected indices
r = downsample(trials.rejected, timeStep);
r = r(i_t,:,:);
assert(size(r,1) == size(data,2));
rrr(1,:,:,:) = r;
rejected = repmat(rrr, [size(freqAxis,2),1,1,1]);
assert(size(rejected,1) == size(data,1));
assert(size(rejected,2) == size(data,2));
assert(size(rejected,3) == size(data,3));
assert(size(rejected,4) == size(data,4));

%% spectral whitening (baseline normalization)
selTime = params.triggering.baseline;
i_t = closestval(timeAxis,selTime(1)):closestval(timeAxis,selTime(2));
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
    
    % baseline definition
    %base = nanmedian(vals,2);    % ~ over time
    base = nanmean(vals,2);    % ~ over time
    assert(any(~isnan(base(1))));
    
    % baseline normalization
    %data(:,:,ch,:) = data(:,:,ch,:)./repmat(base, [1,size(data,2),1,size(data,4)]);
    data(:,:,ch,:) = data(:,:,ch,:) - repmat(base, [1,size(data,2),1,size(data,4)]);
    
    if mod(ch,10) == 0
        display([' - whitening: channel ' num2str(ch) ' of ' num2str(size(trials.data,2)) ' - done.']);
    end        
end

%% log-normal distribution (in dB = 10*log10) (as of 03.07.2018)
%data = 10*log10(data);

%% output
spectra = struct;
spectra.data = data;
spectra.freq = freqAxis;
spectra.time = timeAxis;
spectra.labels = trials.labels;
spectra.clzNames = trials.clzNames;
spectra.rejected = rejected;
