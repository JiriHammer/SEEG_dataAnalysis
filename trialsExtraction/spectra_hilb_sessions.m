function spectra = spectra_hilb_sessions(params, trials)
% computes spectra for each trial using band-pass filtering + Hilbert
% transform + amplitude envevelope
% extracta trials
% trials = 
%           data: [3073x125x500 double]
%           time: [3073x1 double]
%       sessTime: [3073x500 double]
%         labels: [1x500 double]
%     sessNumber: [1x500 double]

% (c) Jiri, Sep21

disp('Computing spectra from Hilbert from sessions ...');

%% load sampling rate
load(params.storage.cacheFile, 'srate');
fs = srate;

%% filter settings -> params
timeStep = 4;                           % in [samples], to reduce memory requirements
params.hilb.keepDim = 'amp';
params.decimate.dsRate = timeStep;      % in [samples]
srate_hilb = srate/timeStep;            % in [Hz]

%% frequencies
FB = [params.hilb_freq.allFreqBands(1:end-1)', params.hilb_freq.allFreqBands(2:end)'];
freqAxis = mean(FB,2)';

%% TO DO: allocation: 4D = [freq x time x chnl x trials]
i_cut = round(params.triggering.time2cut(1)*srate_hilb):round(params.triggering.time2cut(2)*srate_hilb);     % in [samples], w.r.t. cutting point = 0   
t_spectra = downsample(trials.time, timeStep);
assert(length(i_cut) == size(t_spectra,1));
data = nan(length(freqAxis), length(i_cut), size(trials.data,2), size(trials.data,3));
n_tr = 1;
labels = [];
rejected = [];

%% go thru all freq. bands
for sess = 1:size(params.storage.sessionCacheFiles,2)
    
    %% load data
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
    timeAxis = downsample(timeAxis, timeStep);  % downsample by timeStep

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
    
    %% >>> BP + HILB <<<
    freqData = [];
    for freq = 1:size(FB,1)
        disp(['Spectra as BP+HT+ampEnv over sessions, sess = ' num2str(sess) '/' num2str(size(params.storage.sessionCacheFiles,2)) ', freq: ' num2str(freq) '/' num2str(size(FB,1)) ' ...']);

        % band-pass filter + hilbAmpEnv + decimate
        params.downsample = struct(...          % only downsamples (takes every n-th sample)
            'dsUnit', 'samples', ...      % choices: samples, seconds, fromCacheFile (loads timeStep) 
            'dsRate', timeStep ...                   % in [samples]
            );   
        params.decimate = struct(...          % only downsamples (takes every n-th sample)
            'dsUnit', 'samples', ...      % choices: samples, seconds, fromCacheFile (loads timeStep) 
            'dsRate', timeStep ...                   % in [samples]
            );           
        params.bp_freq.freqBand = FB(freq,:);
        params.bp_freq.butterOrder = 3;
        params.thisSess = sess;
        filtData = filterData(params, rawData, {'bp_freq','hilb','log10','z_score'}); % z_score = spectral whitening

        % decimate by timeStep
        decData = nan(size(timeAxis,1),size(rawData,2));
        for ch = 1:size(rawData,2)
            decData(:,ch) = decimate(filtData(:,ch), timeStep);
        end        
                
        % spectral whitening: normalize to mean ampl. envelope over non-rejected samples
%         normFactors = nan(1,size(filtData,2));          % 1 x nCh
%         for ch = 1:size(filtData,2)
%             normFactors(ch) = mean(filtData(~resp_rej(:,ch)),1);
%         end
%         normData(1,:,:) = filtData./repmat(normFactors, [size(filtData,1),1]);

        % spectral whitening: normalize to mean ampl. envelope (all samples)
%         normFactors = mean(filtData,1);
%         normData(1,:,:) = filtData./repmat(normFactors, [size(filtData,1),1]);
        
        % cat to spectra
        tmpData(1,:,:) = decData;
        freqData = cat(1, freqData, tmpData);  % 3D: freq x samples x channels
        clear filtData decData
%         if freq == 10
%             why;
%         end
    end
    
    %% extract trials data -> 4D: freq x samples x channels x trials
    for tr = 1:size(triggers,2)
        cutPoint = triggers(1,tr);                          % in [s]
        i_cutPoint = closestval(timeAxis, cutPoint);        % in [samples]
        inds = i_cutPoint + i_cut;                          % samples to extract
        if inds(1) > 0 && inds(end) <= size(timeAxis,1)
            data(:,:,:,n_tr) = freqData(:,inds,:);          % 4D: freq x samples x channels x trials
            rejected = cat(3, rejected, resp_rej(inds,:));  % 3D: samples x channels x trials
            labels = cat(2, labels, triggers(2,tr));        % class labels
            n_tr = n_tr+1;
        else
            disp(['WARNING: skipping trial (out of session bounds): session = ' num2str(sess) ', trial = ' num2str(tr)]);
        end
    end
    
    disp([' - session ' num2str(sess) ' done.']);
end

%% update srate to cacheFile
srate = srate_hilb;
save(params.storage.cacheFile, 'srate', '-append');

%% log-normal distribution (in dB = 10*log10)
% data = 10*log10(data);

%% rejected indices
assert(size(rejected,1) == size(data,2));
rrr(1,:,:,:) = rejected;
rejected = repmat(rrr, [size(freqAxis,2),1,1,1]);   % same for all freq
assert(size(rejected,1) == size(data,1));
assert(size(rejected,2) == size(data,2));
assert(size(rejected,3) == size(data,3));
assert(size(rejected,4) == size(data,4));

%% output
spectra = struct;
spectra.data = data;
spectra.freq = freqAxis;
spectra.time = t_spectra;
spectra.labels = labels;
spectra.clzNames = trials.clzNames;
spectra.rejected = rejected;
