%% single trials: spectra + raw
% - load trials
% - select channels
% - compute STFT -> spectra

% (c) Jiri, Feb24

%% settings
subjTag = '20_PR11';
selected_channels = {'Rd1-Rd2', 'Ld6-Ld7'};
pathBeg = 'F:\dox\ms_switch_EI\data\v31_stft_session_bip\switchin_EI_IE_bip';
P_level = 0.05;

%% get job-related settings -> params
% --- specify analysis to be performed -> runner
runner = struct;
runner.newJobs = false; % set to 'false' if the results are already computed
runner.newFDR = false;  % true = updates the FDR correction over all patients (and freq. bands, ...) 

% --- job settings: specific (define spectral method) -> params
spectralMethod = 'v31_stft_session_bip';           % ~ v31 (BIP), similar to v12 (STFT base = sess, no zscore, no rej, FB = loFreq0)

params = msSEI_getJobParams(runner, spectralMethod, P_level);

%% load trials
params.storage.subjTag = subjTag;

% --- cache file
cacheFile = [pathBeg filesep subjTag filesep 'cacheFile.mat'];
assert(exist(cacheFile,'file') == 2);
% params.storage.cacheFile = cacheFile;               % used in connectivityAnalysis.m

% load raw trials
clear trials
load(cacheFile, 'trials');  

% load raw trials
clear H selCh_H_resp
load(cacheFile, 'H', 'selCh_H_resp');  
selCh_H = selCh_H_resp;
assert(size(selCh_H,2) == size(trials.data,2));

%% compute STFT -> spectra
params.storage.cacheFile = cacheFile;
params.triggering.doSpectra = 'stft_trials';
% params.triggering.doSpectra = 'mtft_trials';
spectra = spectra_stft_trials(params, trials);    % computes spectra by STFT on trials
spectra.info = trials.info;

%% select channels
i_ch = [];
for c = 1:size(selected_channels,2)
    for ch = 1:size(selCh_H,2)
        if strcmp(H.channels(selCh_H(ch)).name, selected_channels{c})
            i_ch = cat(2, i_ch, ch);
        end
    end
end

%% from trigger analysis
params.plot_triggering.time2plot = [-2.3, 2.3];

spectralData = struct;
spectralData.info.marg_h = [0.1 0.1];
spectralData.info.marg_w = [0.1 0.1];
spectralData.info.gap = [0.01 0.01];
    
% aux. info
spectralData.info.outDir = [pathBeg filesep 'msSEI_fig_sCh_singleTrials' filesep params.triggering.doSpectra];
spectralData.info.cLims = [-5, 5];
i_x = closestval(spectra.time,params.plot_triggering.time2plot(1)):closestval(spectra.time,params.plot_triggering.time2plot(2));
i_f = closestval(spectra.freq,0):closestval(spectra.freq,120);
spectralData.info.chNames = trials.info.chNames(i_ch);

% plot each trial
for tr = 1:size(spectra.data,4)
    [tf, i_c] = ismember(spectra.labels(tr),cell2mat(spectra.clzNames(:,2)));
    assert(tf);
    spectralData.info.text = ['TRIAL: ' num2str(tr) ...
        ', subject: ' subjTag ...
        ', triggered: ' params.triggering.cutPoint ...
        ', class: ' spectra.clzNames{i_c,1}];
    spectralData.info.figName = ['trial_' num2str(tr) '_spectra'];  % use real class names
    spectralData.info.rejected = spectra.rejected(:,:,:,tr);
    spectralData.xVals = spectra.time(i_x);
    spectralData.yVals = spectra.freq(i_f)';
    spectralData.cVals = spectra.data(i_f,i_x,i_ch,tr);                  % 3D: [freq x time x chnls]
    spectralData.info.rejected = spectra.rejected(i_f,i_x,i_ch,tr);                  % 3D: [freq x time x chnls]
    
    % uncomment to add time-domain trial data
    i_t = closestval(trials.time,params.plot_triggering.time2plot(1)):closestval(trials.time,params.plot_triggering.time2plot(2));
    spectralData.trial.data = trials.data(i_t,i_ch,tr);              % add trial below spectrum
    spectralData.trial.time = trials.time(i_t);
    
    % >>> plot <<<
    plotSpectra(params, spectralData);
end

%% plot avg over clz
clzLabels = unique(spectra.labels); 
for clz = 1:length(clzLabels)

    % spectra -> struct spectralData
    spectralData = prepare2plotSpectra(params, spectra, clz);
    spectralData.info.outDir = [pathBeg filesep 'msSEI_fig_sCh_singleTrials' filesep params.triggering.doSpectra];
    spectralData.info.chNames = spectralData.info.chNames(i_ch);
    
    % plot spectralData
    spectralData.cVals = spectralData.cVals(:,:,i_ch);  
    plotSpectra(params, spectralData);    
end
