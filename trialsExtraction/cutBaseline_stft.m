function cutBaseline_stft(params, resp, triggers)
% gets preprocessed data: resp{sess} = [samples x ch]
% & struct triggers. 
%     sampleSessLabel: [3238x3 double]      - triplets = [sample, sess, clz]
%       timeSessLabel: [3238x3 double]      - time in [s]
%       eventDuration: [3238x1 double]      - in [s]
%                info: [1x1 struct]
% computes STFT and extracts its mean of the triggered events

% (c) Jiri, Jul13

srate = params.amp.srate;                               % assumes no prior downsampling!
labels = unique(triggers.sampleSessLabel(:,3));
b = cell(length(labels),1);                             % 3D = [freq, time, chnls] for each class

for sess = 1:size(resp,2)
    sessTime = [1:size(resp{sess},1)]./srate;
    
    % STFT settings
    nfft = 2^nextpow2(params.stft_freq.windowSize * srate);                     % number of FFT points, windowSize in [samples]
    timeStep = ceil(params.stft_freq.timeStep * srate);                         % in [samples], if any downsampling occurs, 'timeStep' MUST exist
    freqAxis = srate/2 * linspace(0,1,nfft/2+1);                                % in [Hz], up to Nyquist freq.
    timeAxis = sessTime(downsample([1:size(resp{sess},1)], timeStep));          % in [s], 0 = sessBeg
    stftData = nan(length(freqAxis), length(timeAxis), size(resp{sess},2));     % 3D = [freq, time, chnls]

    % processing STFT, PSD
    for ch = 1:size(resp{sess},2)
        x = cat(1, flipdim(resp{sess}(1:floor(nfft/2)-1,ch),1), resp{sess}(:,ch), flipdim(resp{sess}(end-floor(nfft/2)+1:end,ch),1));
        [~,~,~,P] = spectrogram(x, hann(nfft,'periodic'), nfft-timeStep, nfft, srate);
        assert(size(P,1) == length(freqAxis));
        assert(size(P,2) == length(timeAxis));        
%         S = spectrogram(x, hann(nfft,'periodic'), nfft-timeStep, nfft, srate);
%         assert(size(S,1) == length(freqAxis));
%         assert(size(S,2) == length(timeAxis));
        stftData(:,:,ch) = P;                                  % PSD
        if mod(ch,10) == 0
            display(['Session: ' num2str(sess) '. Channel: ' num2str(ch) ' done.']);
        end
    end    
            
    for clz = 1:length(labels)
    
        thisClz = labels(clz);
        i_clz  = find(triggers.sampleSessLabel(:,3) == thisClz);
        i_sess = find(triggers.sampleSessLabel(:,2) == sess);
        i_beg = intersect(i_clz, i_sess);                           % event indices in trigger.timeSessLabel corresponding to sess & clz
        i_t = [];
        for t = i_beg'
            i_ev = closestval(timeAxis, triggers.timeSessLabel(t,1)):closestval(timeAxis, triggers.timeSessLabel(t,1)+triggers.eventDuration(t));
            i_t = [i_t, i_ev];
        end
        b{clz} = cat(2, b{clz}, stftData(:,i_t,:));
    end
end

%% freq spectrum of classes (means & std)
base = cell(length(labels),1);                             % 2D = [freq, chnls] for each class
for clz = 1:length(labels)
    base{clz}.avg = squeeze(mean(b{clz}, 2));
    base{clz}.std = squeeze(std(b{clz}, 0, 2));
    base{clz}.freq = freqAxis;
end

%% save
save(params.storage.cacheFile, 'base', '-append');

