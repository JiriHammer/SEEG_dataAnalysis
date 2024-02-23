function trials_out = trialsFiltering(params, trials_in)
% applies post-processing on trials data
% postprocessing defined in: params.triggering.trialsFiltering
% trials_in = struct with fields:
%             data: [79×57×122 double]
%             time: [79×1 double]
%         rejected: [79×57×122 double]
%         sessTime: [4097×121 double]
%           labels: [1×122 double]
%       sessNumber: [1×121 double]
%         clzNames: {2×2 cell}
%     freqBandName: 'beta'
%          figname: '4_beta_20_PR3'
%          selFreq: [13 30]

% (c) Jiri, Sep22

%% default output
trials_out = trials_in;
if ~isfield(params.triggering, 'trialsFiltering')
    return;
end
if isempty(params.triggering.trialsFiltering{1})
    return;
end

%% get current sampling rate -> fs
fs = 1/mean(diff(trials_in.time));       % computed from time vector

%% apply post-processing on trials data
rawData = trials_in.data;
for fltr = 1:size(params.triggering.trialsFiltering,2)
    filterName = params.triggering.trialsFiltering{fltr};
    
    % --- low pass filtering (TBD, now=Dec23)
    if strcmp(filterName, 'lowPass')  
        freqNyquist = fs/2;
        loF = params.bp_freq.freqBand(1);
        hiF = params.bp_freq.freqBand(2);
        Wn = hiF/freqNyquist;                               % normalized bandpass frequencies
        n = 3;                                              % butterworth order
        [b,a] = butter(n, Wn, 'low');                       % returns polynoms of Butterw. filter
        filtData = filtfilt(b, a, rawData);
    end            

    
    % --- band-pass (low-pass) filtering
    if strcmp(filterName, 'bandPass')  
        freqNyquist = fs/2;
        loF = params.bp_freq.freqBand(1);                   % in [Hz], lower  band cutoff
        hiF = params.bp_freq.freqBand(2);                   % in [Hz], higher band cutoff    
        n = 3;                                              % butterworth order
        
        % --- low pass filtering
        if loF == 0  
            Wn = hiF/freqNyquist;                               % normalized bandpass frequencies            
            [b,a] = butter(n, Wn, 'low');                       % returns polynoms of Butterw. filter
            filtData = filtfilt(b, a, rawData);
        else        
            Wn = [loF, hiF]/freqNyquist;                        % normalized bandpass frequencies
            [b,a] = butter(n, Wn);                              % returns polynoms of Butterw. filter
        end
        
        % filtering
        filtData = filtfilt(b, a, rawData);
        
        % uncomment for checking the filter response
%         fhandle = fvtool(b,a, 'Analysis','freq');
%         set(fhandle, 'Fs',fs, 'NormalizedFrequency','off');        
    end
    
    % --- hilbert ampEnv
    if strcmp(filterName, 'hilbAmp') 
        complexData = hilbert(rawData);
        filtData = abs(complexData);
    end
    
    % --- hilbert complexPhase
    if strcmp(filterName, 'hilbComplexPhase') 
        complexData = hilbert(rawData);
        filtData = complexData./abs(complexData);
    end
    
    % --- hilbert: phase angle, [-pi,pi]
    if strcmp(filterName, 'hilbPhaseAngle') 
        complexData = hilbert(rawData);
        filtData = angle(complexData./abs(complexData));
    end
    
    % --- z-score (over all trials)
    if strcmp(filterName, 'z-score') 
        filtData = zscore(rawData,0,1);     % computes zscore over time (dim = 1) for each trial independently!
    end    
    
    % --- flip in time
    if strcmp(filterName, 'flipInTime') 
        filtData = flip(rawData,1);     % flips samples in time (dim = 1) for each trial
    end    
    
    % --- derivative (N point stencils)
    if strcmp(filterName, 'derivative')  
        n_stencils = 3;
        % processing
        if n_stencils == 1
            filtData = diff(rawData)./(1/fs);
            filtData = cat(1, filtData(1,:), filtData);         % same number of elements
        else
            filtData = nan(size(rawData));
            for ch = 1:size(rawData,2)
                for tr = 1:size(rawData,3)
                    filtData(:,ch,tr) = cent_diff_n(rawData(:,ch,tr), 1/fs, n_stencils);
                end
            end
        end
    end
    
    if strcmp(filterName, 'singlePrecision')
        filtData = single(rawData);
    end
    
    assert(size(filtData,1) == size(rawData,1));
    assert(size(filtData,2) == size(rawData,2));
    assert(size(filtData,3) == size(rawData,3));
    rawData = filtData;
end
trials_out.data = filtData;     % update output
