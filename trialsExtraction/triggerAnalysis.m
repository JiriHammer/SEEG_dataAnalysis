function triggerAnalysis(params)
% determines trigger events and uses those samples to extract response in
% given time range
% triggers{clz}(sample,sess) are indicated by METHOD in params.triggerOnMETHOD

% (c) Jiri, Apr12 (in whitenoise), Mar16

if ~params.triggering.doTriggering
    disp('Triggering on events not selected!');
    return;
else
    disp('Triggering on events ...');
end
if ~isfield(params.plot_triggering,'verLines'), params.plot_triggering.verLines = [-1 0 1]; end
if ~isfield(params.plot_triggering,'horLines'), params.plot_triggering.horLines = 0; end

%% get triggers
defineTriggers(params);

%% behavioural analysis
doBehaviourAnalysis(params);

%% extract trials from continuous data (sessions) -> ERP (event-related potential)
trials = extractTrials(params);
if isempty(trials.data)
    disp(['WARNING: ' params.storage.subjTag ': no triggers found !']);
    return;
end

%% ERP: plot class means + SE (plotband)
% arranges data from 'trials' to 'trialsData' structure
trials.figname = ['0_erp_' params.storage.subjTag];
trialsData = prepare2plotTrials(params, trials);
trials.info = trialsData.info;

% plot?
% if params.plot_triggering.trials_avg    
%     plotTrials(params, trialsData);
% end

%% export trials ?
if isfield(params.saveResults, 'exportTriggering_trials')
    if ~isempty(params.saveResults.exportTriggering_trials)
%         exportTrials(params, trials);   % for use for external sharing (BP, ...)
        save(params.storage.cacheFile, 'trials', '-append');
    end
end

%% ERP: plot brain topo in slices
trialsData.info.name_freqBand = 'erp';
plotBrainTopo_trialsData(params, trialsData);

%% ERP: plot single (all) trials (TO DO: add trial-specific paraTimes)
if params.plot_triggering.trials_all
    trials.text = 'ERP';
    trials.outDir = [params.storage.outputDir filesep 'singleTrials'];
    plotSingleTrials(params, trials);
end

%% coherency & imaginary coherency (assumes: trials = ERP (raw) signal)
coherenceAnalysis(params, trials);

%% efective connectivity (DTF, ...) (assumes: trials = ERP (raw) signal)
connectivityAnalysis(params, trials);

%% compute spectra: assumes trials = ERP (raw) signal
if strcmp(params.triggering.doSpectra,'stft_trials') || strcmp(params.triggering.doSpectra,'stft')
    spectra = spectra_stft_trials(params, trials);    % computes spectra by STFT on trials
elseif strcmp(params.triggering.doSpectra,'mtft_trials')
    spectra = spectra_stft_trials(params, trials);    % computes spectra by MTFT on trials (multi-tapering method)
elseif strcmp(params.triggering.doSpectra,'stft_sessions')
    spectra = spectra_stft_sessions(params, trials);  % computes spectra by STFT on sessions, then extracts trials
elseif strcmp(params.triggering.doSpectra,'stft_baseRS')
    spectra = spectra_stft_baseRS(params, trials);    % computes spectra by STFT, baseline from resting state (RS)    
elseif strcmp(params.triggering.doSpectra,'hilb_sessions')
    spectra = spectra_hilb_sessions(params, trials);  % computes spectra by band-pass + Hilbert transform + ampl.envelope on sessions, then extracts trials    
elseif strcmp(params.triggering.doSpectra,'psd_trials')
    spectra = spectra_stft_trials(params, trials);    % computes spectra by STFT on trials (PSD impl. from Radek Janca)
elseif strcmp(params.triggering.doSpectra,'cwt_trials')
    spectra = spectra_stft_trials(params, trials);    % computes spectra by continuous wavelet transform (CWT) on trials          
else
    disp('Spectra computation not selected.');
    return;
end

%% spectra: plot means
if params.plot_triggering.spectra_avg
    clzLabels = unique(spectra.labels); 
    spectra.info = trialsData.info;
    for clz = 1:length(clzLabels)

        % spectra -> struct spectralData
        spectralData = prepare2plotSpectra(params, spectra, clz);

        % plot spectralData
        plotSpectra(params, spectralData);    

        % save
        outputVarName = ['spectralData_clz' num2str(clz)];
        eval([outputVarName '=spectralData;' ]); 
        disp(['Saving spectral data, class: ' num2str(clz) ' (' spectra.clzNames{clz,1} ') ...']);
        save(params.storage.cacheFile, outputVarName, '-append');          
    end        
end

%% spectra: plot single (all) trials  (TO DO: add trial-specific paraTimes)
if params.plot_triggering.spectra_all
    
    % aux. info
    spectralData.info.outDir = [params.storage.outputDir filesep 'singleTrials'];
    spectralData.info.cLims = [-5, 5];
    i_x = closestval(spectra.time,params.plot_triggering.time2plot(1)):closestval(spectra.time,params.plot_triggering.time2plot(2));

    % plot each trial
    for tr = 1:size(spectra.data,4)
        [tf, i_c] = ismember(spectra.labels(tr),cell2mat(spectra.clzNames(:,2)));
        assert(tf);
        spectralData.info.text = ['TRIAL: ' num2str(tr) ...
            ', subject: ' params.storage.subjTag ...
            ', triggered: ' params.triggering.cutPoint ...
            ', class: ' spectra.clzNames{i_c,1}];
        spectralData.info.figName = ['trial_' num2str(tr) '_spectra'];  % use real class names
        spectralData.info.rejected = spectra.rejected(:,:,:,tr);
        spectralData.xVals = spectra.time(i_x);
        spectralData.cVals = spectra.data(:,i_x,:,tr);                  % 3D: [freq x time x chnls]
        spectralData.info.rejected = spectra.rejected(:,i_x,:,tr);                  % 3D: [freq x time x chnls]
        % uncomment to add time-domain trial data
%         i_t = closestval(trials.time,params.plot_triggering.time2plot(1)):closestval(trials.time,params.plot_triggering.time2plot(2));
%         spectralData.trial.data = trials.data(i_t,:,tr);              % add trial below spectrum
%         spectralData.trial.time = trials.time(i_t);
        plotSpectra(params, spectralData);
    end
end

%% frequency bands: plot spectral power
for fb = 1:size(params.triggering.freqBands,1)

    % select frequencies
    freqBand = params.triggering.freqBands{fb,2};
    freqBandName = params.triggering.freqBands{fb,1};
    i_fr = closestval(spectra.freq,freqBand(1)):closestval(spectra.freq,freqBand(2));

    % re-define the trials struct
    trials.time = spectra.time;
%     d = squeeze(mean(spectra.data(i_fr,:,:,:),1));          % 3D = time x chnls x clz
%     i_t = closestval(trials.time,params.triggering.baseline(1)):closestval(trials.time,params.triggering.baseline(2));
%     r = squeeze(spectra.rejected(1,:,:,:));
%     b = mean(d(i_t,:,:),1);                                 % baseline (mean over time)
%     d = d - repmat(b, [size(d,1),1,1]);                     % baseline subtraction
%     trials.data = d;
    trials.data = squeeze(mean(spectra.data(i_fr,:,:,:),1));  % 3D = time x chnls x trials
    assert(size(trials.data,1) == size(trials.time,1));
    trials.rejected = squeeze(spectra.rejected(1,:,:,:));   % assumes same rejection for all freq. bins of spectrograms
    trials.labels = spectra.labels;

    % arranges data from 'trials' to 'trialsData' structure
    trials.selFreq = freqBand;
    trials.freqBandName = freqBandName;
    trials.figname = [num2str(fb) '_' freqBandName '_' params.storage.subjTag];
%     if strcmp(freqBandName,'hiGamma')
%         why;
%     end
    trials = trialsFiltering(params, trials);
    trialsData = prepare2plotTrials(params, trials);

    % plot: AVG +/- SEM
    if params.plot_triggering.trials_avg
        plotTrials(params, trialsData);  
    end

    % plot brain topo in slices
    trialsData.info.name_freqBand = freqBandName;
    plotBrainTopo_trialsData(params, trialsData);        
    
%     if params.plot_triggering.brainTopo
%     
%         % brain info
%         if ~isfield(plotInfo, 'brain')
%             plotInfo.file2load = 'wT1.nii';      % choices: wrCT, wT1, wT2, wc1T1, wc2T2, ...
%             plotInfo.brain = getBrainData(params, plotInfo);
%         end
% 
%         % plot
%         time2plot = params.plot_brainTopo.time2plot;        % in [s], w.r.t. cutting point 
%         plotInfo.clims = [prctile(trialsData.yVals(:),0.5), prctile(trialsData.yVals(:),99.5)];
%         for clz = 1:size(trialsData.yVals,3)
%             for t = 1:size(time2plot,2)
%                 i_t = closestval(trialsData.xVals, time2plot(t));
%                 vals = trialsData.yVals(i_t,:,clz);
%                 plotInfo.text = ['subject: ' params.storage.subjTag ...
%                                  ', triggered: ' params.triggering.cutPoint ...
%                                  ', class: ' trials.clzNames{clz,1} '(' num2str(trialsData.info.nClz(clz)) ')' ...
%                                  ', time = ' num2str(time2plot(t)) ' s'];   
%                 plotInfo.outDir =[params.storage.outputDir filesep 'topoBrain_' freqBandName filesep 'class_' num2str(clz) '_' trials.clzNames{clz,1}];
%                 plotInfo.figName=['fig_' num2str(i_t)];
%                 plot_brainSlices(params, vals(:), plotInfo);
%             end
%         end
%     end     

    % plot: single (all) trials
    if params.plot_triggering.trials_all
        trials.text = params.triggering.freqBands{fb,1};
        trials.outDir = [params.storage.outputDir filesep 'singleTrialPlots' filesep params.triggering.freqBands{fb,1}];
        plotSingleTrials(params, trials);
    end        
end

%% plot stacked (single) trials for all freq. bands (ch-by-ch)
if isfield(params.plot_triggering, 'stackedTrials')
    if params.plot_triggering.stackedTrials
        plotSingleChannel_stackedTrials(params, spectra, spectralData);
    end
end






