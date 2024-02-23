function signalSpectra(plotInfo)
% loads training sessions and computes specified signal type the PSD
% - FFT
% - STFT

% (c) Jiri, Feb14
% known "little bug": bad recording epochs get rejected first, then PSD (-> "jumps in data")

%% settings
if ~isfield(plotInfo, 'signalType')
    plotInfo.signalType = 'EcogGrid';
end

if ~isfield(plotInfo, 'preprocess')
    plotInfo.preprocess = {'car','bp_hiPass', 'norm2std'};      % typical for ECoG
end

if ~isfield(plotInfo, 'srateFromCache')
    plotInfo.srateFromCache = false;                            % to get srate after processing
end

subjInfos = plotInfo.subjInfos;
pathBeg   = plotInfo.pathBeg;
outName  = plotInfo.outName;
signalType = plotInfo.signalType;
spectralMethod = plotInfo.spectralMethod;
wSize = plotInfo.wSize;                     % in [s]


%% load Traj, FFT or STFT & saves to new file: 'fftTraj.mat' or 'stftTraj_wS##.mat'
for subj = 1:size(subjInfos,1)
    outFile = [pathBeg '/' outName '/' subjInfos{subj,1} '/' subjInfos{subj,2} '/' signalType '_stftPSD_wS' num2str(wSize) '.mat'];
    if exist(outFile,'file') == 2
        display(['PSD file found: ' outFile]);
                
    else
        
        % load params from cache file
        cacheFile = [pathBeg '/' outName '/' subjInfos{subj,1} '/' subjInfos{subj,2} '/cache.mat'];
        display(['loading: ' cacheFile ' ...']);
        assert(exist(cacheFile,'file') == 2);
        clear params;
        load(cacheFile, 'params');
        assert(exist('params','var') == 1);
        if isstruct(params.connectionTable)
            params.connectionTable = eval(['C' params.connectionTable.name ';']);
        end
        if plotInfo.srateFromCache
            load(cacheFile, 'sampleRateAfterFiltering');
            srate = sampleRateAfterFiltering;
        else
            srate = params.amp.srate;
        end

        % load Traj from training sessions -> raw{sess}
        raw = cell(length(params.simulation.trainSession),1);
        for sess = 1:length(params.simulation.trainSession)
            thisSess = params.simulation.trainSession(sess);
            sessionCache = [pathBeg '/' outName '/' subjInfos{subj,1} '/' subjInfos{subj,2} '/session' num2str(thisSess) '/cache.mat'];
            assert(exist(sessionCache,'file') == 2);
            load(sessionCache, signalType);
            assert(exist(signalType,'var') == 1);
            raw{sess} = eval(signalType);
        end

        % ------------------------- CAR, hiPass -------------------------------
        params.signalType = signalType;
        params.bp_hiPass = struct(...
            'filterType', 'high', ...
            'butterOrder', 4, ...               % butterworth filter order
            'freqBand', [0.1 Inf] ...           % high pass: [hiCut, Inf] in Hz  
            );    
        [filterSettings ~] = getFeatureMethod(params, plotInfo.preprocess);             % also to reject bad recordings epochs !!!
        if ~isfield(params, 'rejection')                                                            % hack !!!   this field must exist (for outdated jobs)
            params.rejection = [];
        end    
        for sess = 1:length(params.simulation.trainSession)
            params.thisSess = params.simulation.trainSession(sess);
            raw{sess} = filterData(params, raw{sess}, filterSettings, false);
        end   

        % ------------------------------FFT------------------------------------
        if strcmp(spectralMethod, 'fft')
            freqAxis = [0.01:0.01:100];
            data = nan(size(freqAxis,2), size(raw,1));
            for sess = 1:size(raw,1)

                % FFT
                NFFT = 2^nextpow2(size(raw{sess}, 1));                          % number of FFT points
                fftData = fft(raw{sess}, NFFT)/size(raw{sess},1);        % normalizes to windowSize = all samples
                fftData = 2*abs(fftData(1:NFFT/2+1)).^2;                        % PSD, single-sided
                sessFreqAxis = srate/2 * linspace(0,1,NFFT/2+1);

                % smooth FFT result
                Wn = 10/(srate/2);                               % normalized bandpass frequencies
                [b,a] = butter(3, Wn, 'low');
                smoothFFT = filtfilt(b,a, fftData);
    %                 smoothFFT = fftData;

                % resample to freqAxis
                fs = timeseries(smoothFFT, sessFreqAxis);
                fsResampled = resample(fs, freqAxis);
                %figure; hold on; plot(sessFreqAxis, smoothFFT, 'b'); plot(freqAxis, fsResampled.data, 'r');
                data(:,sess) = fsResampled.data;
            end

            % output smoothed FFT results
            x.xVals = freqAxis;
            x.yVals = mean(data, 2);
            x.yErrs = std(data, 0, 2)./(size(data,2)^0.5);
            %figure; plotband(x.xVals(1:100), x.yVals(1:100), x.yErrs(1:100));

            % save to new file: fftTraj.mat
            outFile = [pathBeg '/' outName '/' subjInfos{subj,1} '/' subjInfos{subj,2} '/' signalType '_fftPSD.mat'];
            save(outFile, 'x');  
        end

        % -----------------------------STFT------------------------------------
        if strcmp(spectralMethod, 'stft')
            %wSize = 20;                                                        % in [s], defined thru plotInfo
            tStep = wSize/4;                                                    % in [s]
            nfft = 2^nextpow2(wSize * srate);                                   % number of FFT points, windowSize in [samples]
            timeStep = ceil(tStep * srate);                                     % in [samples], if any downsampling occurs, 'timeStep' MUST exist
            freqAxis = srate/2 * linspace(0,1,nfft/2+1);                        % in [Hz], up to Nyquist freq.

            % allocation
            tmp = [];
            ch = 1;
            for sess = 1:size(raw,1)  
                [~,~,~,P] = spectrogram(raw{sess}(:,ch), hann(nfft,'periodic'), nfft-timeStep, nfft, params.amp.srate);
                assert(size(P,1) == length(freqAxis));
                tmp = cat(2, tmp, P);
            end  
            segmData = nan(size(freqAxis,2), size(raw{sess},2), size(tmp,2));  % ~ [freq x ch x segments]
            segmData(:,ch,:) = tmp;
            clear tmp;

            for ch = 2:size(raw{sess},2)
                n = 1;
                for sess = 1:size(raw,1)   
                    [~,~,~,P] = spectrogram(raw{sess}(:,ch), hann(nfft,'periodic'), nfft-timeStep, nfft, params.amp.srate);
                    assert(size(P,1) == length(freqAxis));
                    segmData(:,ch,n:n+size(P,2)-1) = P;
                    n = n+size(P,2);
                end  
            end
            i_nan = isnan(squeeze(segmData(1,1,:)));
            %segmData(:,:,i_nan) = [];
            assert(sum(i_nan) == 0);

            % output smoothed FFT results
            clear x;
            x.xVals = freqAxis;
            x.yVals = mean(segmData, 3);
            x.yStds = std(segmData, 0, 3);
            x.yErrs = std(segmData, 0, 3)./(size(segmData,3)^0.5);
            %figure; plotband(x.xVals, x.yVals, x.yErrs);

            % save to new file: fftTraj.mat
            outFile = [pathBeg '/' outName '/' subjInfos{subj,1} '/' subjInfos{subj,2} '/' signalType '_stftPSD_wS' num2str(wSize) '.mat'];
            save(outFile, 'x');  

        end
        display(['Spectrum of ' signalType ' saved into: ' outFile]);
    end
    
end
display('pre-loading of PSD done!');

