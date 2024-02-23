function plotSpectralTrials(params, resp, triggers)
% given cell array resp{sess} of 2D data: [samples, ch] computes the STFT
% normalizes spectra to baseline time to compensate for the 1/f decay
% for each class & channel extracts + plots mean of trials 
% the trials are cutted from triggers{1}.sampleSessLabel
% further, plots mean psd over freq bands

% (c) Jiri, Nov12, Mar13
    
%% channel names and pos
signalType = params.triggering.response.signalType;
chNamesAndPos = getChNamesAndPos(params.connectionTable, signalType);
chNames = getChannelNames(params.connectionTable, signalType);
selCh = setdiff(1:size(chNames,2), params.triggering.response.excludedChnls);
chNames = chNames(selCh);
nCh = length(chNames);

clims = [0.8, 1.2];

%% settings & STFT
load(params.storage.cacheFile, 'sampleRateAfterFiltering');
assert(sampleRateAfterFiltering == params.amp.srate);
wSize = 0.25;                                                                   % in [s]
timeStep = ceil(0.01*sampleRateAfterFiltering);                                 % in [samples], if any downsampling occurs, 'timeStep' MUST exist
selFreq = [0, 240];

nfft = 2^nextpow2(wSize * sampleRateAfterFiltering);                            % number of FFT points, windowSize in [samples]
freqAxis = sampleRateAfterFiltering/2 * linspace(0,1,nfft/2+1);                 % in [Hz], up to Nyquist freq.
i_freq = closestval(freqAxis, selFreq(1)):closestval(freqAxis, selFreq(2));     % inds of selected frequencies
        
% allocation of spectral trials
rawTime = [1:size(resp{1},1)]./params.amp.srate;
timeAxis = downsample(rawTime, timeStep);                           % in [s]
triggerPoint = rawTime(triggers{1}.sampleSessLabel(50,1));          % dummy ansatz
i_time = closestval(timeAxis, triggerPoint+params.triggering.time2cut(1)):closestval(timeAxis, triggerPoint+params.triggering.time2cut(2))-1;

% number of labels
labels = unique(triggers{1}.sampleSessLabel(:,3));
binCenters = triggers{1}.info.binning.binCenters;
assert(length(binCenters) == length(labels) - strcmp(params.tuning.rejection{1},'excludeZeroVelocity'));

% selected labels
selLabels = labels;                                            % all
selBinCenters = binCenters;

% allocation of psdData, plots psd(freq,ch)
psdData = cell(1,length(labels));
colors = colorPalette(length(labels));

%% processing STFT, amplitude spectra, plot classMeans: ch-by-ch
if ismember('imagescClassMeans', params.triggering.spectralTrails)
    for clz = 1:length(selLabels)
        clzLabels = find(triggers{1}.sampleSessLabel(:,3) == selLabels(clz));
        % figure
        f = figure('visible', 'on', 'Position', [1, 1, 1920, 1200]);
        set(f,'PaperUnits','inches','PaperPosition',[0 0 19 12]);

        % text
        tx =axes('visible','off', 'position',[0 0 1 1]);
        mytitle = ['Triggering on ' params.triggering.trigMethod ' of: ' triggers{1}.info.chPos '. Response signalType: ' signalType '. Bin center: ' num2str(selBinCenters(clz),3) '. Number of trials: ' num2str(length(clzLabels))];
        mytitle = strrep(mytitle, '_','\_');
        text(0.016, 0.97, mytitle, 'fonts', 14, 'fontw', 'bold'); 

        psdData{clz}.xVals = freqAxis(i_freq);
        psdData{clz}.yVals = nan(length(i_freq), nCh);
        psdData{clz}.yErrs = nan(length(i_freq), nCh);
        psdTimeRange = [-0.1, 0.1];                                                     % in [s], w.r.t. extraction time point
        spectra_avg = nan(length(i_freq),length(i_time),size(resp{1},2));               % 3D: [freq x time x chnls]
        spectra_std = nan(length(i_freq),length(i_time),size(resp{1},2));               % 3D: [freq x time x chnls]
        for ch = 1:size(resp{1},2)
            spectraTrials = nan(length(i_freq), length(i_time), length(clzLabels));      % 4D: [freq, time, trial]
            skipTrial = [];
            nTr = 1;

            for sess = 1:size(resp,2)
                rawData = resp{sess};
                rawTime = [1:size(rawData,1)]./params.amp.srate;
                timeAxis = downsample(rawTime, timeStep);   % in [s]

                x = cat(1, flipdim(rawData(1:floor(nfft/2)-1,ch),1), rawData(:,ch), flipdim(rawData(end-floor(nfft/2)+1:end,ch),1));
                S = spectrogram(x, hann(nfft,'periodic'), nfft-timeStep, nfft, sampleRateAfterFiltering);
                clear x;
                assert(size(S,1) == length(freqAxis));
                assert(size(S,2) == length(timeAxis));
                amps = abs(S(i_freq,:))./(nfft*0.5);                                        % amplitude density
                clear S;

                % normalization to baseline (=median of freq. bin over time)
                i_base = closestval(timeAxis, 5.0):closestval(timeAxis, timeAxis(end)-5.0);
                %normFactors = median(amps(:,i_base), 2);
                normFactors = mean(amps(:,i_base), 2);
                amps = amps./repmat(normFactors, 1, size(amps,2));

                % cutting spectral trials
                i_sess = intersect(find(triggers{1}.sampleSessLabel(:,2) == sess), clzLabels);
                for tr = 1:length(i_sess)
                    triggerPoint = rawTime(triggers{1}.sampleSessLabel(i_sess(tr),1));
                    if triggerPoint+params.triggering.time2cut(1) < 2.0 || triggerPoint+params.triggering.time2cut(2) > timeAxis(end) - 2.0
                        skipTrial = [skipTrial, nTr];
                        %display(['skipping trial: ' num2str(tr) ' from session: ' num2str(sess)]);
                    else
                        i_trial = closestval(timeAxis, triggerPoint+params.triggering.time2cut(1)):closestval(timeAxis, triggerPoint+params.triggering.time2cut(2))-1;
                        spectraTrials(:,:,nTr) = amps(:,i_trial);
                    end
                    nTr = nTr+1;
                end
            end 
            if ~isempty(skipTrial)
                display(['skipped trials: ' num2str(length(skipTrial))]);
                spectraTrials(:,:,skipTrial) = [];
            end

            % channel name and position
            chName = chNames{ch};
            if strcmp(signalType, 'EcogGrid') && ~isempty(findstr(params.job.experimentId, 'WN'))
                aaClr = getChClrOfAa(params.connectionTable, signalType, selCh(ch));
            else
                aaClr = [0 0 0];  % ~ black
            end
            [left bottom width height] = getBoxPosForName(chNamesAndPos, chName);

            % new subplot
            axes('OuterPosition',[left bottom width height], 'visible','on', 'Color','w');
            set(gca, 'XColor',aaClr, 'YColor',aaClr, 'Box','on', 'LineWidth',1.5);        
            hold on;

            % !!! plot !!!
            trialTime = timeAxis(i_trial) - timeAxis(i_trial(1)) + params.triggering.time2cut(1);
            meanSpectrum = mean(spectraTrials,3);
            imagesc(trialTime, freqAxis(i_freq), meanSpectrum, clims);
            %imagesc(trialTime, freqAxis(i_freq), meanSpectrum);
            if ch == size(resp{1},2)
                colorbar;
            end
            axis tight;
            plot([0 0], get(gca,'ylim'), '--k');        % ver. bars @ 0 s
            title([chName '(' num2str(selCh(ch)) ')']);

            % psd over freq
            i_t = closestval(trialTime, psdTimeRange(1)):closestval(trialTime, psdTimeRange(2));
            psdData{clz}.yVals(:,ch) = mean(meanSpectrum(:,i_t), 2);
            psdData{clz}.yErrs(:,ch) = std(squeeze(mean(spectraTrials(:,i_t,:), 2)), 0, 2)./(size(spectraTrials,3).^0.5);     % mean over time, std over trials + SEM
            spectra_avg(:,:,ch) = mean(spectraTrials,3);
            spectra_std(:,:,ch) = std(spectraTrials,0,3);
            display(['Channel: ' num2str(ch) ' done.']);
        end

        % save
        outDir = [params.storage.outputDir '/triggersSpectra_' triggers{1}.info.chPos];
        if ~exist(outDir, 'dir')
            mkdir(outDir);
        end      
        figname = ['classLabel_' num2str(clz)];
        saveas(f, [outDir filesep figname '.fig']);
        %print(f, '-dpng','-r100', [outDir filesep figname '.png']);
        print(f, '-dtiff', '-noui', '-r0', [outDir filesep figname]);
        close(f);  

        % spectral trials
        psdData{clz}.spectra_avg = spectra_avg;
        psdData{clz}.spectra_std = spectra_std;
        psdData{clz}.time = trialTime;
        psdData{clz}.freq = freqAxis(i_freq);
        psdData{clz}.info.trigCh = triggers{1}.info.chPos;
        psdData{clz}.info.signalType = signalType;
        psdData{clz}.info.classBinCenter = selBinCenters(clz);
        psdData{clz}.info.nTrials = length(clzLabels);
        psdData{clz}.info.clr = colors(clz,:);
    end

    display(['figures stored in: ' outDir]);
end

%% TO DO?: plot psd over freq for each class
if ismember('plotbandClassMeans', params.triggering.spectralTrails)
    psdData{1,1}.info.signalType = params.triggering.response.signalType;
    psdData{1,1}.info.selCh = selCh;        
    psdData{1,1}.info.verLines = [8, 25, 50, 200];
    psdData{1,1}.info.horLines = 0;
    %trialsData{1,1}.info.yLims = [-1.2, 1.2];
    psdData{1,1}.info.xlabel = 'freq [Hz]';
    psdData{1,1}.info.text = [params.job.experimentId '. ' params.triggering.response.signalType ' triggered on local extrema of: ' triggers{1}.info.chPos '. PSD over frequency bands. Bandwidth = 10 Hz'];
    psdData{1,1}.info.figname = ['psdFreq_' triggers{1}.info.chPos];
    psdData{1,1}.info.outDir = outDir;
    psdData{1,1}.info.yLims = [0.5, 1.5];
    %topoplotMeanWithError(params, psdData); 
    save(params.storage.cacheFile, 'psdData', 'params', '-append');       % save mean trials
end


%% cut&save all trials (based on pauses!) for each chnl
if ismember('saveSpectra_pauses', params.triggering.spectralTrails)
    time2cut = [-0.5, 2.0];
    triggerPoint = rawTime(triggers{1}.sampleSessLabel(50,1));          % dummy ansatz
    i_time = closestval(timeAxis, triggerPoint+time2cut(1)):closestval(timeAxis, triggerPoint+time2cut(2))-1;
    baselineSpectra = nan(length(i_freq), size(resp,2), size(resp{1},2));                           % 3D: [freq x sess x chnls]
    for ch = 1:size(resp{1},2)
        spectraTrials = nan(length(i_freq), length(i_time), size(triggers{1}.pauseDuration,1));         % 3D: [freq, time, trial]
        skipTrial = [];
        nTr = 1;

        for sess = 1:size(resp,2)
            rawData = resp{sess};
            rawTime = [1:size(rawData,1)]./params.amp.srate;
            timeAxis = downsample(rawTime, timeStep);   % in [s]

            x = cat(1, flipdim(rawData(1:floor(nfft/2)-1,ch),1), rawData(:,ch), flipdim(rawData(end-floor(nfft/2)+1:end,ch),1));
            S = spectrogram(x, hann(nfft,'periodic'), nfft-timeStep, nfft, sampleRateAfterFiltering);
            clear x;
            assert(size(S,1) == length(freqAxis));
            assert(size(S,2) == length(timeAxis));
            amps = abs(S(i_freq,:))./(nfft*0.5);                                        % amplitude density
            clear S;

            % normalization to baseline (=mean of freq. bin over time)
            i_base = closestval(timeAxis, 5.0):closestval(timeAxis, timeAxis(end)-5.0);
            %normFactors = median(amps(:,i_base), 2);
            normFactors = mean(amps(:,i_base), 2);
            amps = amps./repmat(normFactors, 1, size(amps,2));
            baselineSpectra(:,sess,ch) = normFactors;

            % cutting spectral trials
            i_sess = find(triggers{1}.sampleSessLabel(:,2) == sess);
            for tr = 1:length(i_sess)
                triggerPoint = rawTime(triggers{1}.sampleSessLabel(i_sess(tr),1));
                if triggerPoint+time2cut(1) < 2.0 || triggerPoint+time2cut(2) > timeAxis(end) - 2.0
                    skipTrial = [skipTrial, nTr];
                    %display(['skipping trial: ' num2str(tr) ' from session: ' num2str(sess)]);
                else
                    i_trial = closestval(timeAxis, triggerPoint+time2cut(1)):closestval(timeAxis, triggerPoint+time2cut(2))-1;
                    spectraTrials(:,:,nTr) = amps(:,i_trial);
                end
                nTr = nTr+1;
            end
        end 
        assert(size(spectraTrials,3) == nTr-1);
        if ~isempty(skipTrial)
            spectraTrials(:,:,skipTrial) = [];
        end

        % save trials
        outputVarName = ['spectraTrials_ch' num2str(ch)];
        eval([outputVarName '=spectraTrials;' ]);    
        save(params.storage.cacheFile, outputVarName, '-append');
        display(['Channel: ' num2str(ch) ' done.']);
        clear(outputVarName);
    end
    display(['skipped trials: ' num2str(length(skipTrial))]);

    % save time, freq, pauses
    frAxis = freqAxis(i_freq);                                                  % in [s], w.r.t. extraction time point
    trialTime = timeAxis(i_trial) - timeAxis(i_trial(1)) + time2cut(1);
    pauseDurations = triggers{1}.pauseDuration;
    pauseDurations(skipTrial) = [];
    assert(size(spectraTrials,3) == size(pauseDurations,1));
    save(params.storage.cacheFile, 'pauseDurations', 'frAxis', 'trialTime', 'baselineSpectra', '-append');
end

%% cut&save all trials (based on pauses!) for each chnl
if ismember('saveFreqBands_pauses', params.triggering.spectralTrails)
    time2cut = [-2.0, 2.0];
    triggerPoint = rawTime(triggers{1}.sampleSessLabel(50,1));          % dummy ansatz
    i_time = closestval(timeAxis, triggerPoint+time2cut(1)):closestval(timeAxis, triggerPoint+time2cut(2))-1;
    
    % number of trials
    nTr = 0;
    n = 0;
    skipTrial = [];
    for sess = 1:size(resp,2)
        rawData = resp{sess};
        rawTime = [1:size(rawData,1)]./params.amp.srate;
        timeAxis = downsample(rawTime, timeStep);   % in [s]

        % cutting spectral trials
        i_sess = find(triggers{1}.sampleSessLabel(:,2) == sess);
        for tr = 1:length(i_sess)
            triggerPoint = rawTime(triggers{1}.sampleSessLabel(i_sess(tr),1));
            if triggerPoint+time2cut(1) < 2.0 || triggerPoint+time2cut(2) > timeAxis(end) - 2.0
                skipTrial = [skipTrial, n];
            else
                nTr = nTr+1;
            end
            n = n+1;
        end
    end
    assert(nTr + length(skipTrial) == size(triggers{1}.pauseDuration,1));
            
    % init
    nFreq = size(params.freqBands,1);
    nCh = size(resp{1},2);
    triggering = cell(nFreq,1);
    for fr = 1:nFreq
        triggering{fr}.avg = nan(length(i_time), nCh, nTr);
        triggering{fr}.freqRange = params.freqBands(fr,:);
    end
    baselineSpectra = nan(length(i_freq), size(resp,2), nCh);                           % 3D: [freq x sess x chnls]
    
    % STFT: ch-by-ch
    for ch = 1:nCh
        spectraTrials = nan(length(i_freq), length(i_time), nTr);                       % 3D: [freq, time, trial]
        n = 1;
        skipTrial = [];
        ntr = 0;

        for sess = 1:size(resp,2)
            rawData = resp{sess};
            rawTime = [1:size(rawData,1)]./params.amp.srate;
            timeAxis = downsample(rawTime, timeStep);   % in [s]

            x = cat(1, flipdim(rawData(1:floor(nfft/2)-1,ch),1), rawData(:,ch), flipdim(rawData(end-floor(nfft/2)+1:end,ch),1));
            S = spectrogram(x, hann(nfft,'periodic'), nfft-timeStep, nfft, sampleRateAfterFiltering);
            clear x;
            assert(size(S,1) == length(freqAxis));
            assert(size(S,2) == length(timeAxis));
            amps = abs(S(i_freq,:))./(nfft*0.5);                                        % amplitude density
            clear S;

            % normalization to baseline (=mean of freq. bin over time)
            i_base = closestval(timeAxis, 5.0):closestval(timeAxis, timeAxis(end)-5.0);
            %normFactors = median(amps(:,i_base), 2);
            normFactors = mean(amps(:,i_base), 2);
            amps = amps./repmat(normFactors, 1, size(amps,2));
            baselineSpectra(:,sess,ch) = normFactors;

            % cutting spectral trials
            i_sess = find(triggers{1}.sampleSessLabel(:,2) == sess);
            for tr = 1:length(i_sess)
                triggerPoint = rawTime(triggers{1}.sampleSessLabel(i_sess(tr),1));
                if triggerPoint+time2cut(1) < 2.0 || triggerPoint+time2cut(2) > timeAxis(end) - 2.0
                    skipTrial = [skipTrial, n];
                else
                    ntr = ntr+1;
                    i_trial = closestval(timeAxis, triggerPoint+time2cut(1)):closestval(timeAxis, triggerPoint+time2cut(2))-1;
                    spectraTrials(:,:,ntr) = amps(:,i_trial);
                end
                n = n+1;
            end
        end
        assert(ntr+length(skipTrial) == size(triggers{1}.pauseDuration,1));

        % freq bands modulations
        for fr = 1:nFreq
            i_fr = closestval(freqAxis,params.freqBands(fr,1)):closestval(freqAxis,params.freqBands(fr,2));
            triggering{fr}.avg(:,ch,:) = mean(spectraTrials(i_fr,:,:),1);
            triggering{fr}.trialTime = timeAxis(i_trial) - timeAxis(i_trial(1)) + time2cut(1);
        end
        display(['Channel: ' num2str(ch) ' done.']);
    end
    display(['skipped trials: ' num2str(length(skipTrial))]);
        
    % save time, freq, pauses
    pauseDurations = triggers{1}.pauseDuration;
    pauseDurations(skipTrial) = [];
    assert(size(spectraTrials,3) == size(pauseDurations,1));
    triggering{1}.pauseDurations = pauseDurations;
    triggering{1}.base = baselineSpectra;
    save(params.storage.cacheFile, 'triggering', '-append');
end
