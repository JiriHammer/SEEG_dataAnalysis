function spectra = spectra_stft_baseRS(params, trials)
% computes spectra for each trial using STFT
% baseline = resting state (if resting state not found, use spectra_stft_sessions)
% trials = 
%           data: [3073x125x500 double]
%           time: [3073x1 double]
%       sessTime: [3073x500 double]
%         labels: [1x500 double]
%     sessNumber: [1x500 double]

% (c) Jiri, Apr22

disp('Computing spectra from resting state ...');

%% load sampling rate (after filtering = preprocessing)
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

%% STFT baseline (from resting state)
[B_avg, B_std] = stftBaseline_fromRS(params);
if isempty(B_avg)
    disp('STFT baseline from RS: not found.');
    disp(' - baseline method: stft_sessions');
    params.triggering.doSpectra = 'stft_sessions';
    spectra = spectra_stft_sessions(params, trials);
    return;
end
save(params.storage.cacheFile, 'B_avg', 'B_std', '-append');

%% allocation: 4D = [freq x time x chnl x trials]
% i_cut = ceil(params.triggering.time2cut(1)*srate_stft):ceil(params.triggering.time2cut(2)*srate_stft);     % in [samples], w.r.t. cutting point = 0   
i_cut = round(params.triggering.time2cut(1)*srate_stft):round(params.triggering.time2cut(2)*srate_stft);     % in [samples], w.r.t. cutting point = 0   
t_spectra = downsample(trials.time, timeStep);
assert(length(i_cut) == size(t_spectra,1));
data = nan(length(i_fr), length(i_cut), size(trials.data,2), size(trials.data,3));
n_tr = 1;
labels = [];
rejected = [];

%% processing STFT (short-time Fourier transform)
B_sess_avg = [];
B_sess_std = [];
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
        params.downsample.dsUnit = 'fromCacheFile';
        params.thisSess = sess;
        timeAxis = filterData(params, timeAxis, {'downsample'}, 'from_params.init_srate');
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
    
    % >>> STFT <<<
    stftData = [];      % 3D: freq x samples x chnls
    assert(size(rawData,2) == size(B_avg,3));
    for ch = 1:size(rawData,2)
        x = cat(1, flipdim(rawData(1:floor(nfft/2)-1,ch),1), rawData(:,ch), flipdim(rawData(end-floor(nfft/2)+1:end,ch),1));
        [~,~,~,P] = spectrogram(x, hann(nfft), nfft-timeStep, nfft, srate);
        assert(size(P,1) == length(freqAxis));
        assert(size(P,2) == length(timeAxis));
        amps = 10*log10(P(i_fr,:));
        
        % baseline of the session
%         B_sess_avg_ch = mean(amps(:,i_base), 2);
%         B_sess_std = std(amps(:,i_base), 0, 2);
        
        % normalization to baseline -> minus B_avg 
        assert(size(amps,1) == size(B_avg,1));
        amps = amps - repmat(B_avg(:,1,ch), [1, size(amps,2)]);
        
        % z-score: divide by B_std (??? is this needed ??? now = Sep23)
        assert(size(amps,1) == size(B_std,1));
        amps = amps ./ repmat(B_std(:,1,ch), [1, size(amps,2)]);
        
        % cat over chnls to stft data
        stftData = cat(3, stftData, amps);  % 3D: freq x time x chnls
    end
    assert(size(stftData,2) == size(resp_rej,1));
    
    % extract trials data -> 4D: freq x samples x channels x trials
    for tr = 1:size(triggers,2)
        cutPoint = triggers(1,tr);                          % in [s]
        i_cutPoint = closestval(timeAxis, cutPoint);        % in [samples]
        inds = i_cutPoint + i_cut;                          % samples to extract
        if inds(1) > 0 && inds(end) <= size(timeAxis,1)
            data(:,:,:,n_tr) = stftData(:,inds,:);          % 4D: freq x samples x channels x trials
            rejected = cat(3, rejected, resp_rej(inds,:));  % 3D: samples x channels x trials
            labels = cat(2, labels, triggers(2,tr));        % class labels
            n_tr = n_tr+1;
        else
            disp(['WARNING: skipping trial (out of session bounds): session = ' num2str(sess) ', trial = ' num2str(tr)]);
        end
    end
    
    disp([' - session ' num2str(sess) ' done.']);
end

%% rejected indices
assert(size(rejected,1) == size(data,2));
rrr(1,:,:,:) = rejected;
rejected = repmat(rrr, [size(i_fr,2),1,1,1]);   % same for all freq
assert(size(rejected,1) == size(data,1));
assert(size(rejected,2) == size(data,2));
assert(size(rejected,3) == size(data,3));
assert(size(rejected,4) == size(data,4));

%% output
spectra = struct;
spectra.data = data;
spectra.freq = freqAxis(i_fr);
spectra.time = t_spectra;
spectra.labels = labels;
spectra.clzNames = trials.clzNames;
spectra.rejected = rejected;

%% TO DO: plot baseline vs. mean over session
% trialsData = struct;
% trialsData.xVals = freqAxis(i_fr)';
% trialsData.yVals = squeeze(B_avg);
% trialsData.yErrs = squeeze(B_std);
% params.nCh = size(trialsData.yVals,2);
% trialsData.info.chNames = getChannelNames(params, lower(params.response.signalType), 'resp');
% trialsData.info.xlabel = 'freq [Hz]';
% trialsData.info.verLines = [10, 30, 50, 100];
% trialsData.info.text = ['subj = ' params.storage.subjTag ', baseline from resting state, duration = ' num2str(T(end)/60) ' min'];
% trialsData.info.outDir = [params.storage.outputDir filesep 'baseline_RS'];
% trialsData.info.figName = 'baseline_RS';
% plotTrials(params, trialsData);
