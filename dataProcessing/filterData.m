function filtData = filterData(params, rawData, filterTags, init_srate)
% given 'rawData' (samples x ch) and 'filterSettings' struct containing all
% relevant information about the filtering routine, filters the raw data
% and outputs as 'filtData'
% 'filterSettings' is a cell struct, may contain several filters followed
% by one another, ex.: {'norm2std', 'band pass', 'hilbert', 'amp'}

% (c) Jiri, Jan11

%% get starting sampling rate -> fs in [Hz]
if nargin < 4
    init_srate = 'from_sessionCacheFile';   % default
end

if isnumeric(init_srate)
    fs = init_srate;    % in [Hz], defined from above
elseif strcmp(init_srate, 'from_params.init_srate')
    fs = params.init_srate;  % from initial sampling rate
elseif strcmp(init_srate, 'from_sessionCacheFile')
    % starting / current sampling frequency of raw / processed data
    clear fs srate;
    load(params.storage.sessionCacheFiles{params.thisSess}, 'srate');   % load from sessionCache!
    fs = srate;
end
disp([' - data processing, inital sampling rate = ' num2str(fs) ' Hz']);

%% processing chain of filters
for tag = 1:size(filterTags,2)
    filterName = filterTags{tag};
    filterSettings = getfield(params, filterTags{tag});

    % filter rawData
    if strcmp(filterName, 'raw')
        disp([' - data processing, method: ' filterName ' detected.']);
        % processing
        disp(' - data processing: raw data kept.');
        filtData = rawData;
        
    elseif strcmp(filterName, 'norm2std')
        disp([' - data processing, method: ' filterName ' detected.']);
        % processing
        s = std(rawData, 0, 1);
        filtData = rawData./repmat(s, size(rawData,1),1);  % norm2std

    elseif strcmp(filterName, 'norm2iqr')
        disp([' - data processing, method: ' filterName ' detected.']);
        % processing
        s = prctile(rawData, 75, 1) - prctile(rawData, 25, 1);  % inter-quartile range (iqr)
        filtData = rawData./repmat(s, size(rawData,1),1);  % norm2iqr
        
    elseif strcmp(filterName, 'z_score')
        disp([' - data processing, method: ' filterName ' detected.']);
        % processing
        filtData = zscore(rawData, 0, 1);
        
    elseif strcmp(filterName, 'zeroMean')
        disp([' - data processing, method: ' filterName ' detected.']);
        % processing
        m = mean(rawData,1);
        filtData = (rawData - repmat(m, size(rawData,1),1));   

    elseif strcmp(filterName, 'zeroOut_positive')
        disp([' - data processing, method: ' filterName ' detected.']);
        % processing
        filtData = rawData;   
        filtData(filtData > 0) = 0;   % all positive values = 0

    elseif strcmp(filterName, 'zeroOut_negative')
        disp([' - data processing, method: ' filterName ' detected.']);
        % processing
        filtData = rawData;   
        filtData(filtData < 0) = 0;   % all positive values = 0
        
    elseif strcmp(filterName, 'downsample')
        disp([' - data processing, method: ' filterName ' detected.']);
        % settings -> dsRate (in samples)
        if strcmp(filterSettings.dsUnit, 'samples')
            dsRate = filterSettings.dsRate;                             % in [samples]
            timeStep = dsRate;
            save(params.storage.cacheFile, 'timeStep', '-append');      % if any downsampling occurs, 'timeStep' in [samples] MUST exist
        elseif strcmp(filterSettings.dsUnit, 'seconds')
            dsRate = filterSettings.dsRate*fs;                          % in [samples]
            if dsRate ~= round(dsRate)                                  % not an integer
                disp(['WARNING: downsampling factor is not an integer = ' num2str(dsRate)]);
                disp(['WARNING: rounding to = ' num2str(round(dsRate))]);
                dsRate = round(dsRate);                                 % in [samples]
            end
            timeStep = dsRate;
            save(params.storage.cacheFile, 'timeStep', '-append');      % if any downsampling occurs, 'timeStep' in [samples] MUST exist            
        elseif strcmp(filterSettings.dsUnit, 'fromCacheFile')
            clear timeStep;
            load(params.storage.cacheFile, 'timeStep');
            assert(exist('timeStep','var') == 1);
            dsRate = timeStep;                                          % in [samples]
        else
            error('uknown downsampling rate');
        end
        % processing -> timeAxis, filtData, fs (updated)
        filtData = downsample(rawData, dsRate);
        fs = fs/dsRate;                       % !!! adjust srate   
        
    elseif strcmp(filterName, 'decimate')
        disp([' - data processing, method: ' filterName ' detected.']);
        % settings -> dsRate (in samples)
        if strcmp(filterSettings.dsUnit, 'samples')
            dsRate = filterSettings.dsRate;                             % in [samples]
            timeStep = dsRate;
            save(params.storage.cacheFile, 'timeStep', '-append');       % if any downsampling occurs, 'timeStep' in [samples] MUST exist
        elseif strcmp(filterSettings.dsUnit, 'seconds')
            dsRate = filterSettings.dsRate*fs;                             % in [samples]
            if dsRate ~= round(dsRate)                                  % not an integer
                disp(['WARNING: downsampling factor is not an integer = ' num2str(dsRate)]);
                disp(['WARNING: rounding to = ' num2str(round(dsRate))]);
                dsRate = round(dsRate);                                 % in [samples]
            end            
            timeStep = dsRate;
            save(params.storage.cacheFile, 'timeStep', '-append');       % if any downsampling occurs, 'timeStep' in [samples] MUST exist            
        elseif strcmp(filterSettings.dsUnit, 'fromCacheFile')
            clear timeStep;
            load(params.storage.cacheFile, 'timeStep');
            assert(exist('timeStep','var') == 1);
            dsRate = timeStep;
        else
            error('uknown downsampling rate');
        end
        % processing -> timeAxis, filtData, fs (updated)
        timeAxis = downsample(1:size(rawData,1), dsRate);
        filtData = nan(size(timeAxis,2),size(rawData,2));
        for ch = 1:size(rawData,2)
            filtData(:,ch) = decimate(rawData(:,ch), dsRate);
        end        
        fs = fs/dsRate;                       % !!! adjust srate   
        
    elseif strcmp(filterName, 'resample')
        disp([' - data processing, method: ' filterName ' detected.']);
        if ~isfield(filterSettings, 'filt_order') filterSettings.filt_order = 100; end  % default
        % allocate?
        if size(rawData,2) > 1
            tmp = resample(rawData(:,1),filterSettings.new_srate,fs,filterSettings.filt_order);
            filtData = nan(size(tmp,1),size(rawData,2));
        end
        % resample
        for ch=1:size(rawData,2)
            filtData(:,ch) = resample(rawData(:,ch),filterSettings.new_srate,fs,filterSettings.filt_order);
        end        
        fs = filterSettings.new_srate;      % change the currect sampling rate 'fs'
        
    elseif strcmp(filterName, 'sgolay')
        disp([' - data processing, method: ' filterName ' detected.']);
        % processing
        wSizeInSamples = ceil(fs * filterSettings.wSize);
        if wSizeInSamples/2 == round(wSizeInSamples/2)      % is even
            wSizeInSamples = wSizeInSamples + 1;
        end
        filtData = sgolayfilt(rawData, filterSettings.order, wSizeInSamples);

    elseif strcmp(filterName, 'asym_sgolay')
        disp([' - data processing, method: ' filterName ' detected.']);
        % processing
        filtData = fsgwindow(rawData,filterSettings.order,filterSettings.wSizeInSamples,filterSettings.lagInSamples); % (c) S.Rotter
        
%         other implementations & debug:
%         sg = sgolayfilt(rawData, filterSettings.order,filterSettings.wSizeInSamples);  % (c) MATLAB
%         xf = fsgwindow(rawData,filterSettings.order,filterSettings.wSizeInSamples,filterSettings.lagInSamples); % (c) S.Rotter
%         figure; hold on; plot(sg(:,1), 'r'); plot(filtData(:,1), 'b'); plot(xf(:,1), 'g'); legend({'sg', 'asg', 'fsg'});

%         tmp = cat(1, flipdim(rawData(1:floor(filterSettings.wSizeInSamples/2),:),1), rawData, flipdim(rawData(end-floor(filterSettings.wSizeInSamples/2):end,:),1));
%         filtData = nan(size(rawData));
%         for t = 1:size(rawData,1)
%             i_t = [1:filterSettings.wSizeInSamples] + (t-1);
%             filtData(t,:) = sum(tmp(i_t,:).*repmat(filterSettings.smoothKernel', 1,size(rawData,2)), 1);
%         end
%         assert(isempty(find(isnan(filtData),1,'first')));

    elseif strcmp(filterName, 'smooth_Hann')        % convolve with Hann window        
        disp([' - data processing, method: ' filterName ' detected.']);
        % processing
        wSizeInSamples = ceil(fs * filterSettings.wSize);
        smoothWindow = hann(wSizeInSamples);
        filtData = conv(rawData, smoothWindow, 'same');
        
    elseif ~isempty(strfind(filterName, 'bp_'))     % = filtfilt.m
        disp([' - data processing, method: ' filterName ', filter type: ' filterSettings.filterType ' pass, detected.']);
        % processing
        freqNyquist = fs/2;
        loF = filterSettings.freqBand(1);
        hiF = filterSettings.freqBand(2);
        % band pass filter
        if strcmp(filterSettings.filterType, 'low')         % low pass filtering
            Wn = hiF/freqNyquist;                               % normalized bandpass frequencies
            n = filterSettings.butterOrder;                     % butterworth order
            [b,a] = butter(n, Wn, 'low');                       % returns polynoms of Butterw. filter
            filtData = filtfilt(b, a, rawData);
            cutFreq = hiF;
        elseif strcmp(filterSettings.filterType, 'high')    % high pass filtering
            Wn = loF/freqNyquist;                               % normalized bandpass frequencies
            n = filterSettings.butterOrder;                     % butterworth order
            [b,a] = butter(n, Wn, 'high');                      % returns polynoms of Butterw. filter
            filtData = filtfilt(b, a, rawData);            
            cutFreq = loF;
        elseif strcmp(filterSettings.filterType, 'stop')    % stop pass filtering (notch)
            Wn = [loF, hiF]/freqNyquist;                        % normalized bandpass frequencies
            n = filterSettings.butterOrder;                     % butterworth order
            [b,a] = butter(n, Wn, 'stop');                              % returns polynoms of Butterw. filter
            filtData = filtfilt(b, a, rawData);
            cutFreq = hiF;            
        elseif strcmp(filterSettings.filterType, 'band')    % band pass filtering
            Wn = [loF, hiF]/freqNyquist;                        % normalized bandpass frequencies
            n = filterSettings.butterOrder;                     % butterworth order
            [b,a] = butter(n, Wn);                              % returns polynoms of Butterw. filter
            filtData = filtfilt(b, a, rawData);
            cutFreq = hiF;
        end
        
        % visualize freq. response of the used filter
        if params.thisSess == 1 && params.plot_filtering.filterDesign
            fhandle = fvtool(b,a, 'Analysis','freq');
            set(fhandle, 'Fs',fs, 'NormalizedFrequency','off');
            set(gca, 'xlim',[0 5*cutFreq]);
            legend(fhandle, filterName);
            figname = filterSettings.prefix;
            figname(end) = [];                          % to delete '_'
            if ~exist(params.storage.outputDir,'dir')
                mkdir(params.storage.outputDir);
            end
            %saveas(fhandle, [params.storage.outputDir filesep figname '.fig']);  % does not work?
            print('-dpng','-r0', [params.storage.outputDir filesep [figname '.png']]);
            close(fhandle);
        end

    elseif ~isempty(strfind(filterName, 'acausal_filter'))     % = filter.m
        disp([' - data processing, method: ' filterName ', filter type: ' filterSettings.filterType ' pass, detected.']);
        % settings
        freqNyquist = fs/2;
        loF = filterSettings.freqBand(1);
        hiF = filterSettings.freqBand(2);
        
        % if FILTER.m is used, mirror the signal
        transientOnset = 60 * fs;     % time for transient onset filter artifact
        if transientOnset > size(rawData,1), transientOnset = 20 * fs; end
        assert(transientOnset < size(rawData,1));
        temp = cat(1, flip(rawData(1:transientOnset,:),1), rawData);  
        
        % filter
        if strcmp(filterSettings.filterType, 'low')         % low pass filtering
            Wn = hiF/freqNyquist;                               % normalized bandpass frequencies
            n = filterSettings.butterOrder;                     % butterworth order
            [b,a] = butter(n, Wn, 'low');                       % returns polynoms of Butterw. filter
            filtData = filter(b, a, temp);
            cutFreq = hiF;
        elseif strcmp(filterSettings.filterType, 'high')    % high pass filtering
            Wn = loF/freqNyquist;                               % normalized bandpass frequencies
            n = filterSettings.butterOrder;                     % butterworth order
            [b,a] = butter(n, Wn, 'high');                      % returns polynoms of Butterw. filter
            filtData = filter(b, a, temp);
            cutFreq = loF;
        elseif strcmp(filterSettings.filterType, 'stop')    % stop pass filtering (notch)
            Wn = [loF, hiF]/freqNyquist;                        % normalized bandpass frequencies
            n = filterSettings.butterOrder;                     % butterworth order
            [b,a] = butter(n, Wn, 'stop');                              % returns polynoms of Butterw. filter
            filtData = filter(b, a, temp);
            cutFreq = hiF;            
        elseif strcmp(filterSettings.filterType, 'band')    % band pass filtering
            Wn = [loF, hiF]/freqNyquist;                        % normalized bandpass frequencies
            n = filterSettings.butterOrder;                     % butterworth order
            [b,a] = butter(n, Wn);                              % returns polynoms of Butterw. filter
            filtData = filter(b, a, temp);
            cutFreq = hiF;
        end
        filtData(1:transientOnset,:) = [];
        
        % visualize freq. response of the used filter
        if params.thisSess == 1 && params.plot_filtering.filterDesign
            fhandle = fvtool(b,a, 'Analysis','freq');
            set(fhandle, 'Fs',fs, 'NormalizedFrequency','off');
            set(gca, 'xlim',[0 5*cutFreq]);
            legend(fhandle, filterName);
            figname = filterSettings.prefix;
            figname(end) = [];                          % to delete '_'
            if ~exist(params.storage.outputDir,'dir')
                mkdir(params.storage.outputDir);
            end
            %saveas(fhandle, [params.storage.outputDir filesep figname '.fig']);  % does not work?
            print('-dpng','-r0', [params.storage.outputDir filesep [figname '.png']]);
            close(fhandle);
        end
        
    elseif strcmp(filterName, 'hilb')
        disp([' - data processing, method: ' filterName ' (' filterSettings.keepDim ') detected.']);
        % processing
        complexData = hilbert(rawData);
        if strcmp(filterSettings.keepDim, 'real')
            filtData = real(complexData);
        elseif strcmp(filterSettings.keepDim, 'imag')
            filtData = imag(complexData);
        elseif strcmp(filterSettings.keepDim, 'complex')
            filtData = complexData;
        elseif strcmp(filterSettings.keepDim, 'amp')
            filtData = abs(complexData);
        elseif strcmp(filterSettings.keepDim, 'phase')  % complex, phase kept, amp = 1 (scaled by amplitude envelope)
            ampEnvelope = abs(complexData);
            filtData = complexData./ampEnvelope;  
        elseif strcmp(filterSettings.keepDim, 'realAndImag')
            filtData = cat(2, real(complexData), imag(complexData));              
        else
            error('unknown feature');
        end
        clear complexData;
        
    elseif strcmp(filterName, 'fft')
        disp([' - data processing, method: ' filterName ' detected.']);
        % processing FFT, amplitude spectra
        NFFT = 2^nextpow2(size(rawData, 1));                    % number of FFT points
        fftData = fft(rawData, NFFT)/(0.5*size(rawData, 1));    % normalizes to windowSize = all samples
        filtData = abs(fftData(1:NFFT/2+1,:,:)).^2;             % PSD

    elseif strcmp(filterName, 'stft_freq')
        disp([' - data processing, method: ' filterName ' detected.']);
        if ~isfield(filterSettings, 'log10_trafo')
            filterSettings.log10_trafo = false;
        end
        if isempty(params.stft_freq.allFreqBands)
            nBands = 1;
        else
            nBands = size(params.stft_freq.allFreqBands,1);
            selFreq = params.stft_freq.allFreqBands;
            save(params.storage.cacheFile, 'selFreq', '-append');
        end

        % settings
        nfft = 2^nextpow2(filterSettings.windowSize * fs);            % number of FFT points, windowSize in [samples]
        timeStep = ceil(filterSettings.timeStep * fs);                % in [samples], if any downsampling occurs, 'timeStep' MUST exist
        freqAxis = fs/2 * linspace(0,1,nfft/2+1);                     % in [Hz], up to Nyquist freq.
        timeAxis = downsample(1:size(rawData,1), timeStep);
        filtData = nan(length(timeAxis), size(rawData,2), nBands);      % 3D or 2D: samples x chnls x freqBands (if selected)

        i_loF = closestval(freqAxis, filterSettings.freqBand(1,1));
        i_hiF = closestval(freqAxis, filterSettings.freqBand(1,2));
        assert(i_loF <= i_hiF);
        selFreq = setdiff(i_loF:i_hiF, closestval(freqAxis, 48):closestval(freqAxis, 52));      % selected frequencies, leave out 50 Hz noise from amps !
        
        % processing STFT, single freq. band, PSD
        for ch = 1:size(rawData,2)
            x = cat(1, flipdim(rawData(1:floor(nfft/2)-1,ch),1), rawData(:,ch), flipdim(rawData(end-floor(nfft/2)+1:end,ch),1));
            %[~,~,~,P] = spectrogram(x, hann(nfft,'periodic'), nfft-timeStep, nfft, sampleRateAfterFiltering);
            [~,~,~,P] = spectrogram(x, hann(nfft), nfft-timeStep, nfft, fs);
            assert(size(P,1) == length(freqAxis));
            assert(size(P,2) == length(timeAxis));
            if isempty(params.stft_freq.allFreqBands)
                amps = P(selFreq,:);
            else
                amps = P;
            end
            
            % baseline definition
            i_base = stftBaselineIndices(params, size(amps,2));

            % normalization to baseline (=median of freq. bin over time)
            normFactors = median(amps(:,i_base), 2);
            %normFactors = mean(amps, 2);
            amps = amps./repmat(normFactors, 1, size(amps,2));

            if filterSettings.log10_trafo
                amps = 10*log10(amps);                                                                                    % log of PSD ?
            end    

            % average across frequency range
            if isempty(params.stft_freq.allFreqBands)
                filtData(:,ch) = mean(amps, 1)';
            else
                for fb = 1:size(params.stft_freq.allFreqBands,1)
                    selBand = params.stft_freq.allFreqBands(fb,:);
                    selFreq = closestval(freqAxis, selBand(1)):closestval(freqAxis, selBand(2));      % selected frequencies, leave out 50 Hz noise from amps !
                    filtData(:,ch,fb) = mean(amps(selFreq,:), 1)';
                end
            end   
            if mod(ch,10) == 0
                disp(['  - session: ' num2str(params.thisSess) ', channel: ' num2str(ch) ' done.']);
            end
        end
 
        fs = fs/timeStep;                       % !!! adjust srate -> step-wise processing        
        save(params.storage.cacheFile, 'freqAxis', 'nfft', 'timeStep', '-append');
        
    elseif ~isempty(strfind(filterName, 'stft_freqComp'))
        disp([' - data processing, method: ' filterName ' detected.']);
        
        % settings
        nfft = 2^nextpow2(filterSettings.windowSize * fs);            % number of FFT points, windowSize in [samples]
        timeStep = ceil(filterSettings.timeStep * fs);                % in [samples], if any downsampling occurs, 'timeStep' MUST exist
        freqAxis = fs/2 * linspace(0,1,nfft/2+1);                     % in [Hz], up to Nyquist freq.
        
        % processing STFT, amplitude spectra
        timeAxis = downsample(1:size(rawData,1), timeStep);
        filtData = nan(length(timeAxis), size(rawData,2), length(filterSettings.freqBins));  % output = 3D: time x chnls x freq-FC
        for ch = 1:size(rawData,2)
            x = cat(1, flipdim(rawData(1:floor(nfft/2)-1,ch),1), rawData(:,ch), flipdim(rawData(end-floor(nfft/2)+1:end,ch),1));  % mirror the ends
            %S = spectrogram(x, hann(nfft), nfft-timeStep, nfft, fs);
            [S,F,T,P] = spectrogram(x, hann(nfft), nfft-timeStep, nfft, fs);        % hack: compute & save PSD (Oct16)
            if strcmp(params.thisData, 'resp')
                save(params.storage.sessionCacheFiles{params.thisSess}, 'P', '-append');  
            end
            assert(size(S,1) == length(freqAxis));
            assert(size(S,2) == length(timeAxis));
            if strcmp(filterSettings.keepDim, 'complex')
                timeFreq = S(filterSettings.freqBins,:)';
            elseif strcmp(filterSettings.keepDim, 'amp')
                timeFreq = abs(S(filterSettings.freqBins,:)');
            elseif strcmp(filterSettings.keepDim, 'angle')
                timeFreq = angle(S(filterSettings.freqBins,:)');    
            elseif strcmp(filterSettings.keepDim, 'phase')
                timeFreq = S(filterSettings.freqBins,:)';
                if ~isempty(find(abs(timeFreq(:)) == 0))
                    x = abs(timeFreq(:));
                    i_zeroX = find(x == 0);
                    x(i_zeroX) = x(i_zeroX) + 1e-100;       % add some small number to avoid devision by zero !
                    absTimeFreq = reshape(x, size(timeFreq,1), size(timeFreq,2));
                else
                    absTimeFreq = abs(timeFreq);
                end
                timeFreq = timeFreq./absTimeFreq;
                assert(isempty(find(isnan(timeFreq),1,'first')));
            end
            filtData(:,ch,:) = timeFreq;
            if mod(ch,10) == 0
                disp(['Session: ' num2str(params.thisSess) '. Channel: ' num2str(ch) ' done.']);
            end
            clear S;
        end
        selFreq = freqAxis(filterSettings.freqBins);
        fs = fs/timeStep;                       % !!! adjust srate -> step-wise processing        
        save(params.storage.cacheFile, 'freqAxis', 'nfft', 'selFreq', 'timeStep', '-append');
        
    elseif strcmp(filterName, 'mtft_freq')
        disp([' - data processing, method: ' filterName ' detected.']);
        
        % settings
        nfft = 2^nextpow2(filterSettings.windowSize * fs);            % number of FFT points, windowSize in [samples]
        timeStep = ceil(filterSettings.timeStep * fs);                % in [samples], if any downsampling occurs, 'timeStep' MUST exist
        freqAxis = fs/2 * linspace(0,1,nfft/2+1);                     % in [Hz], up to Nyquist freq. 
        
        % spectra: load from OR save to
        outDir = [filterSettings.fileName num2str(filterSettings.windowSize*1000) filesep 'session_' num2str(params.thisSess)];
        fileName = [outDir filesep 'spectra_srate' num2str(fs) '.mat'];
                
        % processing MTFT, amplitude spectra
        if strcmp(filterSettings.getSpectra, 'loadFromFile')
            load(fileName, 'timeStep');     % time step could be different, if using precomputed results
        end
        timeAxis = downsample(1:size(rawData,1), timeStep);
        filtData = nan(length(timeAxis), size(rawData,2));      % output = 2D: time x chnls
        for ch = 1:size(rawData,2)
            if strcmp(filterSettings.getSpectra, 'computeAndSave')
                % compute
                x = cat(1, flipdim(rawData(1:floor(nfft/2)-1,ch),1), rawData(:,ch), flipdim(rawData(end-floor(nfft/2)+1:end,ch),1));
                NW = 2.5;   % magic parameters from Evariste
                k = 3;                
                [S,frq,ts]=my_mtm_spectrum(x, nfft, timeStep, fs, NW, k);       % !!! mtft !!!
                
                % save S (under 'chX' into '/filePath/session_N/mtft_spectra.mat')
                if ~exist(outDir, 'dir')
                    mkdir(outDir);
                end
                if ~(exist(fileName,'file') == 2)
                    save(fileName, '-v7.3', 'freqAxis', 'timeStep', 'nfft');
                end
                resultVarName = ['ch' num2str(ch)];
                eval([resultVarName '=S;' ]);    
                save(fileName, resultVarName, '-v7.3', '-append');                
            elseif strcmp(filterSettings.getSpectra, 'loadFromFile')
                % load (as 'chX' from '/filePath/session_N/spectra_srateXXX.mat')
                assert(exist(fileName,'file') == 2);
                resultVarName = ['ch' num2str(ch)];
                load(fileName, resultVarName);
                assert(exist(resultVarName,'var') == 1);
                if ch == 1
                    load(fileName, 'freqAxis', 'timeStep', 'nfft');
                end
                
                % evaluate to 'S'
                S = eval(resultVarName);
                clear(resultVarName);
            end
            
            % amplitude density
            assert(size(S,1) == length(freqAxis));
            assert(size(S,2) == length(timeAxis));
            amps = abs(S)./(nfft*0.5);
                
            % normalization to baseline (=median of freq. bin over time)
            normFactors = median(amps, 2);
            amps = amps./repmat(normFactors, 1, size(amps,2));
            
            % average across frequency range
            i_loF = closestval(freqAxis, filterSettings.freqBand(1,1));
            i_hiF = closestval(freqAxis, filterSettings.freqBand(1,2));
            assert(i_loF <= i_hiF);
            filtData(:,ch) = mean(amps(i_loF:i_hiF,:), 1)';
  
            disp(['Session: ' num2str(params.thisSess) ', channel: ' num2str(ch) ' done.']);
        end
        fs = fs/timeStep;                       % !!! adjust srate -> step-wise processing        
        save(params.storage.cacheFile, 'freqAxis', 'nfft', 'timeStep', '-append');                  
        
    elseif strcmp(filterName, 'ampl')
        disp([' - data processing, method: ' filterName ' detected.']);
        % processing
        filtData = abs(rawData);
        
    elseif strcmp(filterName, 'n_power')
        disp([' - data processing, method: ' filterName ' detected.']);
        % processing
        n = filterSettings.n;
        filtData = rawData.^n;        
        
    elseif strcmp(filterName, 'addValue')
        disp([' - data processing, method: ' filterName ' detected.']);
        % processing
        c = filterSettings.val2add;
        filtData = rawData + c;                

    elseif strcmp(filterName, 'angle')
        disp([' - data processing, method: ' filterName ' detected.']);
        % processing
        filtData = angle(rawData);
        
    elseif strcmp(filterName, 'cosine')
        disp([' - data processing, method: ' filterName ' detected.']);
        % processing
        filtData = cos(rawData);        
        
    elseif strcmp(filterName, 'log10')
        disp([' - data processing, method: ' filterName ' detected.']);
        % processing
        filtData = log10(rawData);       
        
    elseif strcmp(filterName, 'random')
        disp([' - data processing, method: ' filterName ' detected.']);
        % processing
        filtData = rand(size(rawData));   
        data_rnd = filtData;
        save(params.storage.sessionCacheFiles{params.thisSess}, 'data_rnd', '-append');  
        clear data_rnd;
        
    elseif strcmp(filterName, 'channelSelection')
        disp([' - data processing, method: ' filterName ' detected.']);
        % processing
        selMethod = filterSettings.channels;
        if strcmp(selMethod, 'all')
            filtData = rawData;
        elseif ismember(selMethod, {'xpos_2D','xvel_2D','xacc_2D'})
            assert(size(rawData,2) == 2);   % cols = x,y
            filtData = rawData(:,1);
        elseif ismember(selMethod, {'ypos_2D','yvel_2D','yacc_2D'})
            assert(size(rawData,2) == 2);   % cols = x,y
            filtData = rawData(:,2);            
        elseif strcmp(selMethod, 'motor')
            clear H;
            load(params.storage.cacheFile, 'H');
            varName = ['selCh_H_' params.thisData];
            clear(varName);
            load(params.storage.cacheFile, varName);
            selCh_H = eval(varName);
            assert(size(rawData,2) == size(selCh_H,2));
            esmLabels = {'hand motor', 'arm motor'};
            chnls_raw = [];
            chnls_H = [];
            for ch = 1:size(rawData,2)
                ch_H = selCh_H(ch);     % index in H.channels
                for lbl = 1:size(esmLabels,2)
                    if ~isempty(strfind(H.channels(ch_H).esm,esmLabels{lbl}))
                        chnls_raw = cat(2, chnls_raw, ch);      % = ch index in rawData
                        chnls_H = cat(2, chnls_H, ch_H);        % = ch index in H.channels
                    end
                end
            end
            chnls_raw = unique(chnls_raw);
            filtData = rawData(:,chnls_raw);
            
            % update selected channels (selCh_H_pred or selCh_H_resp) in H.channels
            if params.thisSess == size(params.storage.sessionCacheFiles,2)      % only after processing the last session!
                chnls_H = unique(chnls_H);           
                eval([varName '=chnls_H;' ]); 
                save(params.storage.cacheFile, varName, '-append');      
            end
        end
        
    elseif strcmp(filterName, 'runavg')
        disp([' - data processing, method: ' filterName ' detected.']);
        n_kernel = ceil(fs * filterSettings.wSize);
        c_kernel = ones(n_kernel,1);            % convolution kernel
        filtData = nan(size(rawData)); 
        for ch = 1:size(rawData,2)
            filtData(:,ch) = conv(rawData(:,ch),c_kernel,'same')./n_kernel;     % avg after convolution
        end
                
    elseif ismember(filterName, {'car','car_1', 'car_2', 'bip', 'bip_eog', 'nan', 'pca'})
        disp([' - data processing, method: ' filterName ' detected.']);
        filterSettings.name = filterName;
        filterMatrix = createSpatialFilter(params, size(rawData,2), filterSettings);
        assert(size(rawData,2) == size(filterMatrix,1));

        % spatially filter the signal
        filtData = rawData * filterMatrix;  
        assert(size(filtData,1) == size(rawData,1));
%         if params.thisSess == size(params.storage.sessionCacheFiles,2)
%             rename_chNames_H(params, filterName, params.thisData);
%         end
        
    elseif strcmp(filterName, 'subtractEcg')
        disp([' - data processing, method: ' filterName ' detected.']);
        
        % subtract mean ECoG response triggered by ECG
        resultVarName = ['ecgTriggeredResponseOn' params.usedSignalType];
        load(params.storage.sessionCacheFiles{params.thisSess}, resultVarName);  % try to load the ecg triggered response
        if exist(resultVarName, 'var') == 1
            ecgResp = eval(resultVarName);
        else
            subtractEcgResponse(params);    % call the method to compute & store the ecg triggered response (for all sessions)
            load(params.storage.sessionCacheFiles{params.thisSess}, resultVarName);
            assert(exist(resultVarName,'var') == 1);
            ecgResp = eval(resultVarName);
        end
        assert(size(ecgResp,2) == size(rawData,2));
        
        filtData = rawData - ecgResp;
        assert(size(filtData,1) == size(rawData,1));

    elseif strcmp(filterName, 'subtractZeroSpeed')
        disp([' - data processing, method: ' filterName ' detected.']);
        fileName = ['/export/jiritmp/webdavmirror/zeroSpeedBaseline/' params.job.experimentId '/' params.job.analysisPrefix '/baseline_zeroSpeed.mat'];
        assert(exist(fileName,'file') == 2);
        load(fileName, 'zeroSpeedBase');
        assert(size(zeroSpeedBase,1) == size(rawData,2));
        filtData = rawData - repmat(zeroSpeedBase', [size(rawData,1),1]);
        
    elseif strcmp(filterName, 'derivative')
        disp([' - data processing, method: ' filterName ' detected.']);
        n_stencils = filterSettings.nPointStencils;
        % processing
        if n_stencils == 1
            filtData = diff(rawData)./(1/fs);
            filtData = cat(1, filtData(1,:), filtData);         % same number of elements
        else
            %filtData = estFirstDerivative(rawData, 1/fs);      % used until (21.07.2017)
            filtData = nan(size(rawData));
            for ch = 1:size(rawData,2)
                filtData(:,ch) = cent_diff_n(rawData(:,ch), 1/fs, n_stencils);
            end
        end
        assert(size(filtData,1) == size(rawData,1));
        
    elseif strcmp(filterName, 'surrogate')
        disp([' - data processing, method: ' filterName ' detected.']);
        % processing
        filtData = surrogate_data(rawData, 0, 0);
        
     elseif strcmp(filterName, 'laplace')
        disp([' - data processing, method: ' filterName ' detected.']);
        % processing
        filtData = nan(size(rawData));
        nRows = length(getRowOrColNames(params.connectionTable, params.usedSignalType, 'Row'));
        nCols = length(getRowOrColNames(params.connectionTable, params.usedSignalType, 'Col'));
        nCh = size(rawData,2);
        for t = 1:size(rawData,1)
            allCh = rawData(t,:);
            gridCh = reshape(allCh, nRows, nCols);
            lapCh = del2(gridCh);
            filtData(t,:) = reshape(lapCh, 1, nCh);
        end
        assert(isempty(find(isnan(filtData),1,'first')));

    elseif strcmp(filterName, 'plot')
        disp([' - data processing, method: ' filterName ' detected.']);
        clear trialsData;
        trialsData = struct;
        trialsData.yVals = rawData;
        trialsData.xVals = [1:size(rawData,1)]'./fs;
        sigma = std(rawData(:), 0, 1);
        trialsData.info.horLines = [-3*sigma, 0, 3*sigma];
        trialsData.info.yLims = getYLims(rawData);        
        trialsData.info.outliers = getOutliers(rawData);
        trialsData.info.text = ['subject: ' params.storage.subjTag ', session: ' num2str(params.thisSess)];
        trialsData.info.outDir = [params.storage.outputDir filesep 'plot_processing'];
        trialsData.info.figName = ['filterNumber' num2str(tag) '_session' num2str(params.thisSess)];
        plotTrials(params, trialsData);
        filtData = rawData;

    elseif strcmp(filterName, 'save2cache')
        disp([' - data processing, method: ' filterName ' detected.']);
        if isempty(filterSettings.addTag)
            tag2add = '_1';
        else
            tag2add = filterSettings.addTag;
        end
        outputVarName = ['filtData' tag2add];
        eval([outputVarName '=rawData;' ]);         
        save(params.storage.sessionCacheFiles{params.thisSess}, outputVarName, '-append');      
        filtData = rawData;
        
    elseif strcmp(filterName, 'combine2complex')
        disp([' - data processing, method: ' filterName ' detected.']);
        % processing
        ch_1 = filterSettings.chnl_1;
        ch_2 = filterSettings.chnl_2;
        filtData = complex(rawData(:,ch_1), rawData(:,ch_2));
        %assert(~isreal(filtData(1)));
        
    elseif strcmp(filterName, 'direction')     % returns angle of direction, in [radians]
        disp([' - data processing, method: ' filterName ' detected.']);
        % processing
        if size(rawData,2) == 2          % ~ 2D task
            x = rawData(:,1);               % typically, rawData = velocity
            y = rawData(:,2);                                                          % 2D speed in [std]            
        elseif size(rawData,2) == 1      % ~ 1D task
            x = rawData(:,1);
            y = zeros(size(rawData,1),1);
        else
            error('unknown dataset. Is it tracker? For 1D or 2D task?');
        end
        xy = complex(x,y);      % complex velocity by adding yVel (or zeros)
        
        % point (0,0) has no direction (avoid division by zero)
        i_zero = abs(xy) == 0;  % avoid division by zero
        
        % unit length complex "directional" vector
        d = zeros(size(rawData,1),1);
        d(~i_zero) = xy(~i_zero)./abs(xy(~i_zero));
        
        % real & imag for x- & y- components of direction
        filtData = zeros(size(rawData,1),size(rawData,2));
        filtData(:,1) = real(d);           % -> x-direction!
        if size(rawData,2) == 2          % ~ 2D task
            filtData(:,2) = imag(d);       % -> y-direction!
        end
                
    elseif strcmp(filterName, 'trafo2trinary')
        disp([' - data processing, method: ' filterName ' detected.']);
        assert(size(rawData,2) == 1);      % ~ 1D task (e.g. xvel for car game)
        % processing
        filtData = zeros(size(rawData));
        i_hi = find(rawData >  filterSettings.threshold);   % threshold in [std]
        i_lo = find(rawData < -filterSettings.threshold);   % threshold in [std]
        filtData(i_hi) =  1;
        filtData(i_lo) = -1;        

    elseif strcmp(filterName, 'shift_subtract')
        disp([' - data processing, method: ' filterName ' detected.']);       
        % processing
        shiftLag = round(filterSettings.timeLag*fs);    % in [s], shiftLag > 0 -> shift of the data into the future
        filtData = circshift(rawData, shiftLag,1);
        if filterSettings.withSubtraction
            filtData = filtData - rawData;              % subtract the raw data -> difference (negative shiftLag gives you the position, where you will be in tau seconds from now)
        end
        
    elseif strcmp(filterName, 'trafo2artificialChnl')
        disp([' - data processing, method: ' filterName ' detected.']);
        % processing
        load(params.storage.sessionCacheFiles{params.thisSess}, filterSettings.signalType);
        Y = eval(filterSettings.signalType);
        i_ch = getDecodingChannels(params.connectionTable, filterSettings.signalType, filterSettings.chName);
        y = Y(:,i_ch);
        y = y./std(y,0,1);
        assert(size(rawData,1) == size(y,1));
        filtData = rawData;
        
        chnls_1 = 10;
        ecogChnls_1 = zeros(size(rawData,1), length(chnls_1));
        i_hi = find(y > 0);
        ecogChnls_1(i_hi,:) = repmat(y(i_hi), 1, length(chnls_1));
        filtData(:,chnls_1) = ecogChnls_1;
        
        chnls_2 = 35;
        ecogChnls_2 = zeros(size(rawData,1), length(chnls_2));
        i_lo = find(y < 0);
        ecogChnls_2(i_lo,:) = repmat(abs(y(i_lo)), 1, length(chnls_2));
        filtData(:,chnls_2) = ecogChnls_2;
        
    elseif strcmp(filterName, 'fakeTriggers')
        filtData = randn(size(rawData));
        assert(strcmp(params.paradigm.usedParadigm,'carDriving'));
        load(params.storage.sessionCacheFiles{params.thisSess}, 'pred_raw');
        trials = pred_raw;
        timeAxis = [1:size(rawData,1)]'./fs;
        
        % REWARDS: higher activation
        for tr = 1:size(trials.rewards,1)
            if trials.rewards{tr}.wasCollected
                t_trigger = trials.rewards{tr}.time; 
                i_spike = closestval(timeAxis, t_trigger);
                i_modul = closestval(timeAxis, t_trigger):closestval(timeAxis, t_trigger+0.2);
%                 if ~isempty(filterSettings.spikes_chnls)
%                     i_chnls = filterSettings.spikes_chnls;
%                     filtData(i_spike,i_chnls) = 1;              % spike
%                 end
                if ~isempty(filterSettings.activityHigh_chnls)
                    i_chnls = filterSettings.activityHigh_chnls;
                    filtData(i_modul,i_chnls) = 2*randn(size(i_modul,2), size(i_chnls,2));  % higher activity  
                end
%                 if ~isempty(filterSettings.activityLow_chnls)
%                     i_chnls = filterSettings.activityLow_chnls;
%                     filtData(i_modul,i_chnls) = 0.5*randn(size(i_modul,2), size(i_chnls,2));  % lower activity  
%                 end               
            end
        end

        % OBSTACLES: spikes & de-activation
        for tr = 1:size(trials.penalty,1)
            if trials.penalty{tr}.wasCollected
                t_trigger = trials.penalty{tr}.time; 
                i_spike = closestval(timeAxis, t_trigger);
                i_modul = closestval(timeAxis, t_trigger):closestval(timeAxis, t_trigger+0.2);
                if ~isempty(filterSettings.spikes_chnls)
                    i_chnls = filterSettings.spikes_chnls;
                    filtData(i_spike,i_chnls) = 1;              % spike
                end
%                 if ~isempty(filterSettings.activityHigh_chnls)
%                     i_chnls = filterSettings.activityHigh_chnls;
%                     filtData(i_modul,i_chnls) = 2*randn(size(i_modul,2), size(i_chnls,2));  % higher activity  
%                 end                
                if ~isempty(filterSettings.activityLow_chnls)
                    i_chnls = filterSettings.activityLow_chnls;
                    filtData(i_modul,i_chnls) = 0.5*randn(size(i_modul,2), size(i_chnls,2));  % lower activity  
                end               
            end
        end         
        
    elseif strcmp(filterName, 'stft2magnInfo')
        disp([' - data processing, method: ' filterName ' detected.']);
        % processing
        i_ch = 2;       % ~ xVel
        v_min = nan(length(params.simulation.trainSession),1);
        for sess = 1:length(params.simulation.trainSession)
            theSess = params.simulation.trainSession(sess);
            load(params.storage.sessionCacheFiles{theSess}, 'Traj');
            v_min(sess) = min(Traj(:,i_ch));
        end
        v_min = min(v_min);
        
        % STFT settings
        ws = 2.0;                                                           % in [s]
        ts = 0.00001;                                                       % in [s]
        nfft = 2^nextpow2(ws * fs);                   % number of FFT points, windowSize in [samples]
        timeStep = ceil(ts * fs);                     % in [samples], if any downsampling occurs, 'timeStep' MUST exist
        freqAxis = fs/2 * linspace(0,1,nfft/2+1);     % in [Hz], up to Nyquist freq.
        freqBins = 1:4*ws;                                                  % which freq. components to keep!
        selFreqBins = 1; %1:4*ws;                                           % which freq. components modify (v_magn info)
        
        % processing STFT
        timeAxis = downsample(1:size(rawData,1), timeStep);
        ch = 1;
        x = cat(1, flipdim(rawData(1:floor(nfft/2)-1,ch),1), rawData(:,ch), flipdim(rawData(end-floor(nfft/2)+1:end,ch),1));  % mirror the ends
        S = spectrogram(x, hann(nfft,'periodic'), nfft-timeStep, nfft, fs);
        assert(size(S,1) == length(freqAxis));
        assert(size(S,2) == length(timeAxis));
        timeFreq = S(freqBins,:)';                          % normalize the ampEnv = 1
        if ~isempty(find(abs(timeFreq(:)) == 0))
            x = abs(timeFreq(:));
            i_zeroX = find(x == 0);
            x(i_zeroX) = x(i_zeroX) + 1e-100;               % add some small number to avoid devision by zero !
            absTimeFreq = reshape(x, size(timeFreq,1), size(timeFreq,2));
        else
            absTimeFreq = abs(timeFreq);
        end
        timeFreq = timeFreq./absTimeFreq; 

        load(params.storage.sessionCacheFiles{params.thisSess}, 'Traj');
        p_amp = 1.05*abs(v_min) + Traj(:,i_ch);             % create new ampEnv = 1.05*|v_min| + v(t)
        assert(size(p_amp,1) == size(timeFreq,1));
        p_ampFC = ones(size(p_amp,1), size(timeFreq,2));
        p_ampFC(:,selFreqBins) = repmat(p_amp, 1, length(selFreqBins));
        y = timeFreq.*p_ampFC;         
        
        % iSTFT
        cumRes = zeros(size(timeFreq,1),1);
        for i_fc = 1:size(timeFreq,2)
            if i_fc == 1                                    % !!! zero DC !!!
                %cumRes = cumRes + real(y(:,i_fc));

            else                                          % !!! reconstruction from STFT !!!
                x = nan(size(y,1),1);                               
                for t = 1:size(y,1)
                    % construct FT-vector 'X' for iFT
                    YS = zeros(length(freqAxis),1);                 % mimic of spectrogram output @ freq.component fc
                    YS(i_fc) = y(t,i_fc);
                    X = cat(1, YS, conj(flipdim(YS(2:end-1),1)));   % mimic of fft output @ freq.component fc: conjugated values flipped around center freq
                    assert(size(X,1) == nfft);

                    % iFT
                    x_tf = ifft(X, nfft);                           % mimic of the time window (input to STFT @ time step t)
                    x(t) = x_tf(nfft/2);                            % take only the middle point
                end
                assert(isempty(find(isnan(x),1,'first')));
                cumRes = cumRes + x;
            end       
        end
        filtData = cumRes./std(cumRes, 0, 1);
        
        clear S;
        fs = fs/timeStep;                       % !!! adjust srate -> step-wise processing
        
    elseif strcmp(filterName, 'model_copyNoise')
        disp([' - data processing, method: ' filterName ' detected.']);
        
        % --- copy 'c' ---
%         c = zeros(size(rawData,1),1);
        for ch = 1:size(filterSettings.copy,2)
            % process copy -> v
            params.shift_subtract = struct(...       % shift the signal by shiftLag
                'timeLag', filterSettings.shiftLag(ch), ...                   % in [s], positive => shift in the future, negative = shift into the past
                'withSubtraction', false ...                % subtract the original (unshifted) signal
                );        
%             copy_filterTags = cat(2, {'bp_loPass','resample','decimate'}, get_filterTags(filterSettings.copy{ch}), {'smooth_Hann','shift_subtract', 'norm2std'});
            copy_filterTags = cat(2, {'bp_loPass','resample','decimate'}, get_filterTags(filterSettings.copy{ch}), {'shift_subtract', 'norm2std'});
%             copy_filterTags = cat(2, {'bp_loPass'}, get_filterTags(filterSettings.copy{ch}), {'smooth_Hann','shift_subtract', 'norm2std'});
            clear tracker
            load(params.storage.sessionCacheFiles{params.thisSess},'tracker');
            v = filterData(params, tracker, copy_filterTags, params.init_srate);
            assert(size(v,1) == size(rawData,1));
            
            % init copy -> c
            if ch == 1
                c = zeros(size(v,1),1);         
            end
            
            if strcmp(filterSettings.withPrep, 'noPrep')
                c = c + filterSettings.weights(ch)*v;               % only the copy activity
            elseif strcmp(filterSettings.withPrep, 'onlyPrep')
                v_prep = preparatoryActivity(params, v, ch, fs);
                c = c + filterSettings.weights(ch)*v_prep;          % only the prep activity
            elseif strcmp(filterSettings.withPrep, 'mixPrep')
                v_prep = preparatoryActivity(params, v, ch);
                v = filterData(params, v, {'smooth_Hann'}, params.init_srate);
                c = c + filterSettings.weights(ch)*(v + v_prep);    % add prep + copy
            elseif strcmp(filterSettings.withPrep, 'absRoad')       % abs val of road = distance to target
                load(params.storage.sessionCacheFiles{params.thisSess},'road');
                copy_filterTags = {'bp_loPass','ampl','shift_subtract', 'norm2std'};
                v_road = filterData(params, road(:,1), copy_filterTags, params.init_srate);
                c = c + filterSettings.weights(ch)*v_road;    
            elseif strcmp(filterSettings.withPrep, 'absRoad_withSpeed')       % abs val of road = distance to target
                load(params.storage.sessionCacheFiles{params.thisSess},'road');
                copy_filterTags = {'bp_loPass','ampl','shift_subtract', 'norm2std'};
                v_road = filterData(params, road(:,1), copy_filterTags, params.init_srate);
                c = c + filterSettings.weights(ch)*v + filterSettings.weights(ch)*v_road;   % add also 'v' (e.g. = speed)                                                      
            end       
        end
%         c = (c-mean(c,1))./std(c,0,1);                      % copy: zeroMean & norm2std
        c = cat(2, c, zeros(size(c,1),1));                  % copy, add zeros as a dummy second channel
        
        % --- noise 'n' ---
        if strcmp(filterSettings.noise, 'ecogNoise')    % TO DO
            error('not done now.');           
        elseif strcmp(filterSettings.noise, 'whiteNoise')
            y = rand(size(c));                                                % white noise <0, 1>
            y = (y-repmat(mean(y,1),[size(y,1),1]))./repmat(std(y,0,1),[size(y,1),1]);  % zeroMean & norm2std
            n = filterSettings.scaleFactor*y;                                           % scale the noise
        elseif strcmp(filterSettings.noise, 'gaussNoise')
            n = filterSettings.scaleFactor*randn(size(c));    % scale the noise                                                        % white noise <0, 1>
        else
            error(['Unknown noise source: ' filterSettings.noise]);
        end
        
        % --- model = copy+noise ---
        filtData = c+n;                                     % ch1: copy+noise, ch2: noise only

    elseif strcmp(filterName, 'model_chnlsCopyNoise')          % LFC(ch,t) = A[ch,veloc(t)]*speed(t) + Noise(sigma)
        disp([' - data processing, method: ' filterName ' detected.']);
        
        % copy 'c'
        load(params.storage.sessionCacheFiles{params.thisSess}, 'Traj');
        ch_copy = getDecodingChannels(params.connectionTable, 'Traj', filterSettings.copy);
        v = Traj(:,ch_copy);
        v = (v-mean(v,1))./std(v,0,1);                                  % zeroMean & norm2std
        c = zeros(size(rawData,1),size(rawData,2));
        c(:,1:filterSettings.nChnls) = repmat(v, [1,filterSettings.nChnls]);
            
        % noise 'n'
        if strcmp(filterSettings.noise, 'whiteNoise')                       % whiteNoise
            n = rand(size(rawData,1),size(rawData,2));                          % white noise <0, 1>                                 % scale the noise
        elseif strcmp(filterSettings.noise, 'ecogNoise_regressOut')         % ECoG noise (regress-out copy variable)
            n = nan(size(rawData,1),size(rawData,2));
            for ch = 1:size(rawData,2)
                beta = regress(rawData(:,ch), cat(2, ones(size(v,1),1), v));        % regression model: predict ecog based on trajCopy
                pred_ecog = cat(2, ones(size(rawData,1),1), v) * beta;              % predicted ecog
                n(:,ch) = rawData(:,ch) - pred_ecog;                                % subtract
            end
        elseif strcmp(filterSettings.noise, 'ecogNoise_surrogate')          % ECoG noise (surrogate distr.)
            n = surrogate_data(rawData, 0, 0);
        else
            error(['Unknown noise: ' filterSettings.noise]);
        end
        
        % scale noise amplitude
        n = (n-repmat(mean(n,1), [size(n,1),1]))./repmat(std(n,0,1), [size(n,1),1]);        % zero mean, unit std     
        n = filterSettings.scaleFactor*n;  
        
        % model = copy+noise (same size as 'ECoG')
        filtData = c+n;
        
    elseif strcmp(filterName, 'model_activPattern')          % LFC(ch,t) = A[ch,veloc(t)]*speed(t) + Noise(sigma)
        disp([' - data processing, method: ' filterName ' detected.']);
        
        % activation pattern 'A(ch)': 
        if strcmp(params.job.experimentId, 'WN1') || strcmp(params.job.experimentId, 'CS4')         % hack !
            ch_L_pos = [3,4,18,19,20,27,28];
            ch_L_neg = [6,7,13,14,22];
            A_L = zeros(1,size(rawData,2));
            A_L(ch_L_pos) = 1;
            A_L(ch_L_neg) = -1;
            ch_R_pos = [18,19,27,28];
            ch_R_neg = [3,4,6,7,13,14,20,22];
            A_R = zeros(1,size(rawData,2));
            A_R(ch_R_pos) = 1;
            A_R(ch_R_neg) = -1;
            
        else
            error(['not defined for subj: ' params.job.experimentId]);
        end
        
        % copy 'c(ch,t)'
        load(params.storage.sessionCacheFiles{params.thisSess}, 'Traj');
        ch_vel = getDecodingChannels(params.connectionTable, 'Traj', 'xVel');
        v = Traj(:,ch_vel);
        v = (v-mean(v,1))./std(v,0,1);                                  % zeroMean & norm2std
        i_L = find(v < 0);
        i_R = find(v > 0);
        c = zeros(size(rawData,1),size(rawData,2));
        c(i_L,:) = abs(v(i_L)) * A_L;
        c(i_R,:) = abs(v(i_R)) * A_R;
            
        % whiteNoise 'n'
        y = rand(size(rawData,1),size(rawData,2));                          % white noise <0, 1>
        y = (y-repmat(mean(y,1), [size(rawData,1),1]))./repmat(std(y,0,1), [size(rawData,1),1]);       % zeroMean & norm2std
        n = filterSettings.scaleFactor*y;                                   % scale the noise

        % model = copy+noise (same size as 'ECoG')
        filtData = c+n;

    elseif strcmp(filterName, 'dipole2ecog')
        disp([' - data processing, method: ' filterName ' detected.']);
        
        % load 2 dipoles
        assert(exist(filterSettings.voltageFile,'file') == 2);
        load(filterSettings.voltageFile, 'voltage');
        tmp = voltage(:,:,1)';      % rows -> cols
        dip1 = tmp(:)';             % cats cols (A1,..,A8,B1,..,B8,C1, ... ,H8) <=> rawData channels
        tmp = voltage(:,:,2)';
        dip2 = tmp(:)';
        
        % polar velocity
        copy_filterTags = {'sgolay','derivative','bp_loPass','norm2std'};
        clear tracker
        load(params.storage.sessionCacheFiles{params.thisSess},'tracker');
        assert(size(tracker,1) == size(rawData,1));
        v = filterData(params, tracker, copy_filterTags);        
        assert(size(v,1) == size(rawData,1));
        
        % modeled ecog activity (sum of 2 dipoles, 1st active for PD, 2nd active 2 anti-PD)
        phi_PD = 0;  %  pi/4;
        filtData = zeros(size(rawData,1), size(dip1,2));
        filtData = filtData + repmat(dip1, [size(rawData,1),1]) .* repmat(abs(v), [1,size(dip1,2)]) .* repmat(0.5*(cos(angle(v)-phi_PD)+1), [1, size(dip1,2)]);
        filtData = filtData + repmat(dip2, [size(rawData,1),1]) .* repmat(abs(v), [1,size(dip2,2)]) .* repmat(0.5*(cos(angle(v)-phi_PD-pi)+1), [1,size(dip2,2)]);
        
    elseif strcmp(filterName, 'broadband')
        % STFT -> PSD(f,t), PCA[PSD(f,t)] -> PSC, projection onto 1st PSC 
        nfft = 2^nextpow2(filterSettings.windowSize * fs);            % number of FFT points, windowSize in [samples]
        timeStep = ceil(filterSettings.timeStep * fs);                % in [samples], if any downsampling occurs, 'timeStep' MUST exist
        freqAxis = fs/2 * linspace(0,1,nfft/2+1);                     % in [Hz], up to Nyquist freq.
        timeAxis = downsample(1:size(rawData,1), timeStep);
        filtData = nan(length(timeAxis), size(rawData,2));
        for ch = 1:size(rawData,2)
            x = cat(1, flipdim(rawData(1:floor(nfft/2)-1,ch),1), rawData(:,ch), flipdim(rawData(end-floor(nfft/2)+1:end,ch),1));
            [~,~,~,P] = spectrogram(x, hann(nfft,'periodic'), nfft-timeStep, nfft, fs);
            assert(size(P,1) == length(freqAxis));
            assert(size(P,2) == length(timeAxis));
            normFactors = mean(P, 2);                                       % baseline
            P = P./repmat(normFactors, 1, size(P,2));                       % norm2base
            P = log10(P);                                                   % log of normalized PSD
            C = cov(P');                                                    % PSD covariance, [freq x freq]
            [eigenvec,eigenval,varExplained] = pcacov(C);                   % PCA -> matrix of principal spectral components
            psc_1 = eigenvec(:,1);                                          % 1st PSC = 1st column
            bb = P' * psc_1;                                                % projection onto PSC_1
            assert(size(bb,1) == length(timeAxis));
            %bb_sgolay = sgolayfilt(bb, 4, 0.02*sampleRateAfterFiltering+1);   % smooth with SGolay
            %bb_gauss = filter(gausswin(0.05*sampleRateAfterFiltering),1,bb);
            bb = zscore(bb);                                                % z-score
            %bb = exp(bb);                                                   % exponentiate
            filtData(:,ch) = bb;
            
            if filterSettings.pscPlot
                f = figure('visible', 'on', 'Position', [1, 1, 1900, 1200], 'Color','w');
                set(f,'PaperUnits','inches','PaperPosition',[0 0 19 12]);                
                clrs = {'c','m','r','b','g','y','k'};
                chNames = getChannelNames(params.connectionTable, params.signalType);
                subplot(2,2,1); hold on;
                legendNames = [];
                for k = 1:4
                    plot(freqAxis, eigenvec(:,k), clrs{k});
                    legendNames = [legendNames, {['PSC ' num2str(k)]}];
                end           
                axis tight;
                plot([50, 50], get(gca,'ylim'), '--k');
                plot(get(gca,'xlim'), [0, 0], '--k');
                xlabel('freq [Hz]');
                legend(legendNames);
                title('PSCs');
                
                subplot(2,2,2);
                hist(bb, 100);
                title('distribution of PSC 1');
                
                subplot(2,2,[3,4]); hold on;
                plot(timeAxis./(fs/timeStep)/60, bb, 'r');
                selFreq = closestval(freqAxis, 60):closestval(freqAxis, 200);
                hiGamma = mean(P(selFreq,:),1);
                plot(timeAxis./(fs/timeStep)/60, zscore(hiGamma), 'b');       % (hiGamma-mean(hiGamma,2))./std(hiGamma,0,2)
                legend([{'bb'},{'hi \gamma'}]);
                xlabel('time [s]');
                ylabel('log(PSD)')
                title('PSC_1');
                
                % text
                tx =axes('visible','off', 'position',[0 0 1 1]);
                mytitle = ['Broadband signal extraction. signal type: ' params.signalType ', channel: ' num2str(ch) '. channelName = ' chNames{ch}];
                mytitle = strrep(mytitle, '_','\_');
                text(0.016, 0.97, mytitle, 'fonts', 14, 'fontw', 'bold');
                
                % save 
                outDir = [params.storage.outputDir '/' 'broadband_psc'];
                if ~exist(outDir, 'dir')
                    mkdir(outDir);
                end                  
                figname = ['ch' num2str(ch) '_sess' num2str(params.thisSess)];
                %saveas(f, [outDir filesep [figname '.fig']]);
                %print(f, '-dtiff', '-noui', '-r0', [outDir filesep figname]);
                print('-dpng','-r0', [outDir filesep [figname '.png']]);
                close(f);                   
            end
            if mod(ch,10) == 0
                disp(['Session: ' num2str(params.thisSess) '. Channel: ' num2str(ch) ' done.']);
            end
        end
        fs = fs/timeStep;                       % !!! adjust srate -> step-wise processing        
        save(params.storage.cacheFile, 'timeStep', '-append');      
        
    elseif strcmp(filterName, 'dist_car_road')
        disp([' - data processing, method: ' filterName ' detected.']);
        % processing
        load(params.storage.sessionCacheFiles{params.thisSess}, 'road', 'gameSpeed');
        assert(exist('road','var') == 1);
        if params.srate ~= fs
            params.downsample.dsRate = 'fromCacheFile';
            road = filterData(params, road, {'downsample'});
        end
        assert(size(road,1) == size(rawData,1));
        y_pixels = [1:size(rawData,1)]'./fs * gameSpeed;    % in [pixels], = time [s] * gameSpeed [pixels/s]
        filtData = mindist_carFromRoad([rawData,y_pixels], [road(:,1),y_pixels], fs);      
        assert(size(filtData,1) == size(rawData,1));
        
    elseif strcmp(filterName, 'suaSimulation_populations')  
        disp([' - data processing, method: ' filterName ' detected.']);
        
        % kinematic parameters
        abspos = abs(rawData./repmat(std(rawData,0,1),[size(rawData,1),1]));    % distance, in [SD]
        phipos = filterData(params, rawData, {'direction','angle'});            % direction of position, in [rad]
        vel = filterData(params, rawData, {'derivative'});                      % velocity
        absvel = abs(vel./repmat(std(vel,0,1),[size(vel,1),1]));                % speed, in [SD]
        phivel = filterData(params, vel, {'direction','angle'});                % direction of position, in [rad]
        
        % neuronal population parameters
        N = max(filterSettings.popSizes);                                       % maximum number of neurons    
        if params.thisSess == 1    
            pos_PD = linTransform(rand(1,N), [0,1],[-pi,pi]);                   % random PDs of position tuning, in [rad]
            vel_PD = linTransform(rand(1,N), [0,1],[-pi,pi]);                   % random PDs of velocity tuning, in [rad]
            save(params.storage.cacheFile, 'pos_PD','vel_PD', '-append');
        else
            clear pos_PD pos_PD;
            load(params.storage.cacheFile, 'pos_PD','vel_PD');
            assert(size(pos_PD,2) == N);
        end
            
        % SUA simulated data
        filtData = zeros(size(rawData,1),size(filterSettings.popSizes,2));
        for ch = 1:size(filterSettings.popSizes,2)
            N = filterSettings.popSizes(ch);                % number of neurons
            if ismember('dist',filterSettings.modelType)    % distance tuning
                G_dist = 0.01;                  % scale factor, in [spikes / s / cm]
                r = G_dist * abspos;            % firing rate of 1 neuron
                P = N * r;                      % population firing
                assert(size(P,1) == size(rawData,1));
                filtData(:,ch) = filtData(:,ch) + P;
            end
            
            if ismember('pos',filterSettings.modelType)    % position tuning
                G_pos = 1;                      % scale factor, in [spikes / s / cm]              
                P = zeros(size(rawData,1),1);
                for n = 1:N
                    r = G_pos * abspos .* cos(phipos - pos_PD(n));  % firing rate of 1 neuron
                    P = P + r;                  % population firing
                end
                assert(size(P,1) == size(rawData,1));
                filtData(:,ch) = filtData(:,ch) + P;
            end
            
            if ismember('spd',filterSettings.modelType)    % distance tuning
                G_spd = 0.01;                    % scale factor, in [spikes / s / (cm/s)]
                r = G_spd * absvel;             % firing rate of 1 neuron
                P = N * r;                      % population firing
                assert(size(P,1) == size(rawData,1));
                filtData(:,ch) = filtData(:,ch) + P;
            end    
            
            if ismember('vel',filterSettings.modelType)    % velocity tuning
                G_vel = 1;                      % scale factor, in [spikes / s / (cm/s)]              
                P = zeros(size(rawData,1),1);
                for n = 1:N
                    r = G_vel * absvel .* cos(phivel - vel_PD(n));  % firing rate of 1 neuron
                    P = P + r;                  % population firing
                end
                assert(size(P,1) == size(rawData,1));
                filtData(:,ch) = filtData(:,ch) + P;
            end
            
            if ismember('noise',filterSettings.modelType)    % noise
                sigmaNoise = 1;
                P = zeros(size(rawData,1),1);
                for n = 1:N
                    P = P + sigmaNoise*randn(size(rawData,1),1);    % population firing
                end
                filtData(:,ch) = filtData(:,ch) + P;
            end
            disp([' - SUA simulation for N = ' num2str(N) ' neurons done.']);
        end

    elseif strcmp(filterName, 'suaSimulation_repetitions')      % repetitions = channels
        disp([' - data processing, method: ' filterName ' detected.']);
        
        % kinematic parameters
        if size(rawData,2) == 1     % car-game (1-D) paradigm  TO DO: !!! check direction filter !!! (was not working with the 2-D RTP)
            abspos = abs(rawData./repmat(std(rawData,0,1),[size(rawData,1),1]));    % distance, in [SD]
            phipos = filterData(params, rawData, {'direction','angle'});            % direction of position, in [rad]
            vel = filterData(params, rawData, {'derivative'});                      % velocity
            absvel = abs(vel./repmat(std(vel,0,1),[size(vel,1),1]));                % speed, in [SD]
            phivel = filterData(params, vel, {'direction','angle'});                % direction of position, in [rad]
        elseif size(rawData,2) == 2     % RTP (2-D) paradigm
            abspos = filterData(params, rawData, {'combine2complex','ampl','norm2std'});               % distance, in [SD]
            phipos = filterData(params, rawData, {'combine2complex','angle'});                    % direction of position, in [rad]
            absvel = filterData(params, rawData, {'derivative','combine2complex','ampl','norm2std'});  % speed
            phivel = filterData(params, rawData, {'derivative','combine2complex','angle'});       % direction of velocity, in [rad]
        end
        
        % neuronal population parameters
        N = filterSettings.popSize;                                             % number of neurons   
        nCh = filterSettings.nRepetitions;                                      % number of repetitions of the experiment with different neuronal populations, realizes as "channels"
        if params.thisSess == 1    
            pos_PD = linTransform(rand(nCh,N), [0,1],[-pi,pi]);                   % random PDs of position tuning, in [rad]
            vel_PD = linTransform(rand(nCh,N), [0,1],[-pi,pi]);                   % random PDs of velocity tuning, in [rad]
            save(params.storage.cacheFile, 'pos_PD','vel_PD', '-append');
        else
            clear pos_PD pos_PD;
            load(params.storage.cacheFile, 'pos_PD','vel_PD');
            assert(size(pos_PD,2) == N);
        end
            
        % gain factors / modulations of firing rate 'r'
        if strcmp('prepAct', filterSettings.gainFactors)    % implements "anticipatory" distance & "executive" speed CN model
            SNR_factor = 0;     % no noise to prep. activity
            params = msDist_get_cnModelSettings(params, 'D_prep_2D', SNR_factor);
            G_dst = filterData(params, rawData(:,1), {'model_copyNoise','norm2std'});
            G_dst = G_dst(:,1);
            params = msDist_get_cnModelSettings(params, 'S_prep_2D', SNR_factor);
            G_spd = filterData(params, rawData(:,1), {'model_copyNoise','norm2std'}); 
            G_spd = G_spd(:,1);
        elseif strcmp('5S_co+S_oc', filterSettings.gainFactors)
            SNR_factor = 0;     % no noise to prep. activity
            params = msDist_get_cnModelSettings(params, '2D_5S_co+S_oc', SNR_factor); 
            G_spd = filterData(params, rawData(:,1), {'model_copyNoise','norm2std'});
            G_spd = G_spd(:,1);
            G_dst = G_spd;      % same gain for position and distance tuning
        elseif strcmp('S_co+S_oc', filterSettings.gainFactors)
            SNR_factor = 0;     % no noise to prep. activity
            params = msDist_get_cnModelSettings(params, '2D_S_co+S_oc', SNR_factor); 
            G_spd = filterData(params, rawData(:,1), {'model_copyNoise','norm2std'});
            G_spd = G_spd(:,1);
            G_dst = G_spd;      % same gain for position and distance tuning            
        else%if strcmp('const', filterSettings.gainFactors)
            G_dst = abspos;
            G_spd = absvel;            
%             G_dst = ones(size(rawData,1),1);
%             G_spd = ones(size(rawData,1),1);
        end
        
        % SUA simulated data
        filtData = zeros(size(rawData,1),nCh);          % population signal: init
        for ch = 1:nCh
            if ismember('D',filterSettings.modelType)   % distance tuning (linear)
                G_d = 0.01;                                      % scale factor, in [spikes / s / cm]
%                 r = G_d * G_dst .* abspos;                    % firing rate of 1 neuron
                r = G_d * G_dst;                                % firing rate of 1 neuron
                P_d = N * r;                                    % population firing
                assert(size(P_d,1) == size(rawData,1));
                filtData(:,ch) = filtData(:,ch) + P_d;          % add to pop. signal
            end
            
            if ismember('P',filterSettings.modelType)   % position tuning (cosine)
                G_p = 1;                                        % scale factor, in [spikes / s / cm]              
%                 P = zeros(size(rawData,1),1);
%                 for n = 1:N
%                     r = G_p * G_dst .* abspos .* cos(phipos - pos_PD(ch,n));  % firing rate of 1 neuron
%                     P = P + r;                              % population firing
%                 end
%                 assert(size(P,1) == size(rawData,1));
%                 filtData(:,ch) = filtData(:,ch) + P;        % add to pop. signal
                
                % matrix version - faster?
                M_G_dst  = repmat(G_dst, [1, size(pos_PD,2)]);   % 2D: t x N_neurons
                M_phipos = repmat(phipos, [1, size(pos_PD,2)]);   % 2D: t x N_neurons
                M_pos_PD = repmat(pos_PD(ch,:), [size(abspos,1), 1]);   % 2D: t x N_neurons
%                 M_abspos = repmat(abspos, [1, size(pos_PD,2)]);   % 2D: t x N_neurons
%                 M_r = G_p * M_G_dst .* M_abspos .* cos(M_phipos - M_pos_PD);  % firing rate of 1 neuron
                M_r = G_p * M_G_dst .* cos(M_phipos - M_pos_PD);  % firing rate of 1 neuron
                P_p = sum(M_r,2);                            % sum over neurons
                assert(size(P_p,1) == size(rawData,1));
                filtData(:,ch) = filtData(:,ch) + P_p;        % add to pop. signal
                clear M_G_dst M_abspos M_phipos M_pos_PD M_r P2;
            end
            
            if ismember('S',filterSettings.modelType)    % speed tuning (linear)
                G_s = 0.01;                                      % scale factor, in [spikes / s / (cm/s)]
%                 r = G_s * G_spd .* absvel;                      % firing rate of 1 neuron
                r = G_s * G_spd;                                % firing rate of 1 neuron
                P_s = N * r;                                    % population firing
                assert(size(P_s,1) == size(rawData,1));
                filtData(:,ch) = filtData(:,ch) + P_s;          % add to pop. signal
            end    
            
            if ismember('V',filterSettings.modelType)    % velocity tuning (cosine)
                G_v = 1;                                  % scale factor, in [spikes / s / (cm/s)]              
%                 P = zeros(size(rawData,1),1);
%                 for n = 1:N
%                     r = G_v * G_spd .* absvel .* cos(phivel - vel_PD(ch,n));  % firing rate of 1 neuron
%                     P = P + r;                              % population firing
%                 end
%                 assert(size(P,1) == size(rawData,1));
%                 filtData(:,ch) = filtData(:,ch) + P;        % add to pop. signal
                
                % matrix version - faster?
                M_G_spd  = repmat(G_spd, [1, size(vel_PD,2)]);   % 2D: t x N_neurons
                M_phivel = repmat(phivel, [1, size(vel_PD,2)]);   % 2D: t x N_neurons
                M_vel_PD = repmat(vel_PD(ch,:), [size(absvel,1), 1]);   % 2D: t x N_neurons
%                 M_absvel = repmat(absvel, [1, size(vel_PD,2)]);   % 2D: t x N_neurons
%                 M_r = G_v * M_G_spd .* M_absvel .* cos(M_phivel - M_vel_PD);  % firing rate of 1 neuron
                M_r = G_v * M_G_spd .* cos(M_phivel - M_vel_PD);  % firing rate of 1 neuron
                P_v = sum(M_r,2);                            % sum over neurons
                assert(size(P_v,1) == size(rawData,1));
                filtData(:,ch) = filtData(:,ch) + P_v;        % add to pop. signal         
                clear M_G_spd M_absvel M_phivel M_vel_PD M_r P2;
            end
            
            if ismember('N',filterSettings.modelType)   % noise
                sigmaNoise = 1;
%                 P = zeros(size(rawData,1),1);
%                 for n = 1:N
%                     P = P + sigmaNoise*randn(size(rawData,1),1);    % population firing
%                 end
                P_n = sum(sigmaNoise*randn(size(rawData,1),N),2);    % population firing
                assert(size(P_n,1) == size(rawData,1));
                filtData(:,ch) = filtData(:,ch) + P_n;        % add to pop. signal
            end
            disp(['   - SUA simulation, ch = ' num2str(ch) '/' num2str(nCh) ' done.']);
   
            % debug figure            
%             fig_make;
%             hold on;
%             t = [1:size(P_n,1)]./fs;
%             plot(t, P_n, 'k');
%             plot(t, P_d, 'b');
%             plot(t, P_s, 'r');
%             plot(t, P_p, 'c');
%             plot(t, P_v, 'm');

        end
        disp([' - SUA simulation for N = ' num2str(N) ' neurons done.']);
        
    elseif strcmp(filterName, 'sineWave')
        t = [1:size(rawData,1)]'./fs;
        f = filterSettings.frequency;
        phi = filterSettings.phase;
        A = filterSettings.amplitude;
        filtData = A*sin(2*pi*f*t-phi);
        
    else
        error(['Unknown filter: ' filterName]);
    end
    
    clear rawData;
    rawData = filtData;
    clear filtData;
end

%% final result
filtData = rawData;

%% update 'srate' ? = sampling rate after filtering (eg. in case of downsampling, step processing - stft, mtft)
clear srate;
load(params.storage.cacheFile, 'srate');
if fs > srate
    disp(['WARNING: Sampling rate after filtering = ' num2str(fs) ' IS LARGER than saved in cacheFile = ' num2str(srate)]);
    disp(' - a feature or a bug?');
end
if strcmp(init_srate, 'from_params.init_srate')  % useful in cases when this function is called outside of the loadAndFilter.m (typical session processing)
    disp('WARNING: Sampling rate after filtering will not be updated.');
else
    disp(['Updating sampling rate after filtering, saving to cacheFile: srate = ' num2str(fs)]);
    srate = fs;
    if params.thisSess == 1
        save(params.storage.cacheFile, 'srate', '-append');     % update sampling rate in cache file
    end
    save(params.storage.sessionCacheFiles{params.thisSess}, 'srate', '-append');     % update sampling rate in session cache file

    % save processed data time axis
    time_prc = [1:size(filtData,1)]'./fs;
    save(params.storage.sessionCacheFiles{params.thisSess}, 'time_prc', '-append');
end

%% debug figure (uncomment & execute in new editor)
% f = fig_make;
% t = [1:size(rawData,1)]./fs;
% ch = 10;
% plot(t, rawData(:,ch), 'r'); 
% hold on;
% plot(t, filtData(:,ch), 'b');
% xlabel('time (s)');
