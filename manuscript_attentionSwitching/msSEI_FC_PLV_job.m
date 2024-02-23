%% ms SEI: PLV (phase locking value) job - connectivity between 2 selected AAs

% (c) Jiri, Dec23

%% settings: v1 - relaxed conditions (p-level, length of significance, ...)
dirSuffix = '_msFig_AA_v1_relaxCond';
t_sel = [-1.0, 2.0];    % the true sel time will be cut by 1 s on edges (= template length): [-2, 2] -> [-1, 1]
P_level = 0.05;         % significance level
t_siglen = 0.10;        % in [s], significance time threshold (ch is signif. if its significance > threshold) 
brainAtlas = 'Yeo7';

%% get job-related settings (params): jobExecutor_SEI.m, v11 -> params
% --- specify analysis to be performed -> runner
runner = struct;
runner.newJobs = false; % set to 'false' if the results are already computed
runner.newFDR = true;  % true = updates the FDR correction over all patients (and freq. bands, ...) 

% --- job settings: specific (define spectral method) -> params
spectralMethod = 'v31_stft_session_bip';           % ~ v31 (BIP), similar to v12 (STFT base = sess, no zscore, no rej, FB = loFreq0)
% spectralMethod = 'v30_stft_baseRS_bip';           % ~ v30 (BIP), similar to v24 (new= no zscore, mo rej, FB = loFreq0), STFT with RS baseline
% spectralMethod = 'v29_stft_baseRS_bip';           % ~ v29 (BIP), similar to v27 (new=zscore, rej = 10 SD, FB = loFreq0), STFT with RS baseline
% spectralMethod = 'v28_stft_baseRS_bip';         % ~ v24 (BIP), STFT with RS baseline
% spectralMethod = 'v27_stft_baseRS_bip';         % ~ v24 (BIP), STFT with RS baseline
% spectralMethod = 'v24_stft_baseRS_car';         % ~ v24 (CAR), STFT with RS baseline
% spectralMethod = 'v24_stft_baseRS_bip';         % ~ v24 (BIP), STFT with RS baseline
% spectralMethod = 'stft_sessions';       % ~ v12
% spectralMethod = 'stft_baseRS';         % ~ v11  (CAR)
% spectralMethod = 'stft_baseRS_bip';     % ~ v11b (BIP)

params = msSEI_getJobParams(runner, spectralMethod, P_level);

%% ================== settings: connectivity -> params ====================
params.connectivity = struct(...
    'computeConnectivity', true, ...
    'connectivityMeasure', {{'DTF'}}, ...   % choice: 'GGC','DTF','ffDTF','PDC','COH'
    'MVAR_order', 20, ...
    'subBands', [1 120], ... settings by RJ: [2 12;13 25;26 48;52 75; 76 98; 102 148];
    'windowSize', 0.5, ...      % in [s]
    'timeStep', 0.125, ...      % in [s], corresponds to N samples at 512 Hz (if N = 64, tStep = 1/(512/64) = 0.125)
    'selectedChnls', 'see below', ... % see below: params.coherence.selectedChnls
    'maxDistToROI', 10, ...      % in [mm], maximum distance (tolerance) to region of interest
    'plotEC_ch2ch', true, ...       % plot coherence: each channel to all channels
    'plotEC_selGroups', true, ...   % plot coherence: selected groups of channels (each group to each group)
    'plotEC_zscore', true, ...       % apply z-score in figs over time dimension
    'time2plot', [-3,3] ...     % in [s]
    );    


% --- v1: PLV
% connectivityMeasure = 'PLV';
% params.connectivity.connectivityMeasure = {connectivityMeasure};
% FC_outDirName = [connectivityMeasure '_v1'];   
% figDirSuffix = 'allFreq_DMN_DAN';
% params.connectivity.selectedChnls = {...    % each row = tuple: atlasName, areaName (must match names in anatomicalAreas_getList.m) 
%     'Yeo7',  'Default',           'was_sgnf', {'loFreq0','loGamma','hiGamma'}; ...
%     'Yeo7',  'Dorsal Attention',  'was_sgnf', {'loFreq0','loGamma','hiGamma'}; ...  
%     }; 
% FB_hilb = [0:4:116]';
% params.connectivity.freqBands = cat(2, FB_hilb, FB_hilb+4);

% --- v2: PLV - in low freq, 1-Hz bands
% connectivityMeasure = 'PLV';
% params.connectivity.connectivityMeasure = {connectivityMeasure};
% FC_outDirName = [connectivityMeasure '_v2'];   
% figDirSuffix = 'allFreq_DMN_DAN';
% params.connectivity.selectedChnls = {...    % each row = tuple: atlasName, areaName (must match names in anatomicalAreas_getList.m) 
%     'Yeo7',  'Default',           'was_sgnf', {'loFreq0','loGamma','hiGamma'}; ...
%     'Yeo7',  'Dorsal Attention',  'was_sgnf', {'loFreq0','loGamma','hiGamma'}; ...  
%     }; 
% FB_hilb = [0:1:13]';
% params.connectivity.freqBands = cat(2, FB_hilb, FB_hilb+1);

% --- v3: PLV - in low freq, 2-Hz bands, single precision
connectivityMeasure = 'PLV';
params.connectivity.connectivityMeasure = {connectivityMeasure};
% FC_outDirName = [connectivityMeasure '_v3'];   
FC_outDirName = [connectivityMeasure '_v5_tbd_test'];   
figDirSuffix = 'allFreq_DMN_DAN';
params.connectivity.selectedChnls = {...    % each row = tuple: atlasName, areaName (must match names in anatomicalAreas_getList.m) 
    'Yeo7',  'Default',           'was_sgnf', {'delta','theta','alpha','beta','loGamma','hiGamma'}; ...
    'Yeo7',  'Dorsal Attention',  'was_sgnf', {'delta','theta','alpha','beta','loGamma','hiGamma'}; ...  
    }; 
FB_hilb = [0:2:120]';
params.connectivity.freqBands = cat(2, FB_hilb, FB_hilb+2);

% --- v4: PLV - in low freq + 4-Hz bands, single precision
% connectivityMeasure = 'PLV';
% params.connectivity.connectivityMeasure = {connectivityMeasure};
% FC_outDirName = [connectivityMeasure '_v4'];   
% figDirSuffix = 'allFreq_DMN_DAN';
% params.connectivity.selectedChnls = {...    % each row = tuple: atlasName, areaName (must match names in anatomicalAreas_getList.m) 
%     'Yeo7',  'Default',           'was_sgnf', {'delta','theta','alpha','beta','loGamma','hiGamma'}; ...
%     'Yeo7',  'Dorsal Attention',  'was_sgnf', {'delta','theta','alpha','beta','loGamma','hiGamma'}; ...  
%     };
% FB_hilb_low = [
%     0 3;
%     3 7;
%     7 12;
%     12 16;
%     ];
% FB_hilb = [16:4:120]';
% params.connectivity.freqBands = cat(1, FB_hilb_low, cat(2, FB_hilb, FB_hilb+4));

%% output directory
% FC_outDir = [params.storage.dir_results filesep params.storage.outName filesep 'DTF_vRJ_outputFiles'];
FC_outDir = [params.storage.dir_results filesep params.storage.outName filesep FC_outDirName];
if ~exist(FC_outDir, 'dir')
    mkdir(FC_outDir);
end    
    
%% compute PLV for each subj
for subj = 1:size(params.storage.subjList,1)
    subjTag = params.storage.subjList{subj,1};
    disp(['FC - DTF: processing, subj = ' subjTag]);
    
    %% --- cache file
    cacheFile = [params.storage.dir_results filesep params.storage.outName filesep subjTag filesep 'cacheFile.mat'];
    assert(exist(cacheFile,'file') == 2);
    params.storage.cacheFile = cacheFile;               % used in connectivityAnalysis.m
    
    % load raw trials
    clear trials
    load(cacheFile, 'trials');  

    %% get selected channels -> chGroups
    params.P_level = P_level;
    params.t_siglen = t_siglen;
    params.nCh = size(trials.data,2);
    chGroups = ch2roi_getSelChnls(params, subjTag);   
    
    %% selected channels -> i_chnls
    i_chnls = [];
    nCh_grp = 0;
    chFoundInEachGroup = true;
    for aa = 1:size(chGroups.nn_name,1)
        if isempty(chGroups.ch_inds_trials{aa})
            chFoundInEachGroup = false;
            break;
        end        
        i_chnls = cat(2, i_chnls, chGroups.ch_inds_trials{aa});
        chGroups.ch_inds_M{aa} = nCh_grp + [1:size(chGroups.ch_inds_trials{aa},2)];
        nCh_grp = size(chGroups.ch_inds_trials{aa},2);
    end
    if ~chFoundInEachGroup
        continue;
    end
    
    %% plot single trials (for rejection)
%     params.storage.subjTag = subjTag;
%     params.H = chGroups.H;
%     params.selCh_H = cat(2, chGroups.ch_inds_H{1}, chGroups.ch_inds_H{2});  % hardcoded for 2 clz
%     params.nCh = size(i_chnls,2);
%     trials_plot = trials;
%     trials_plot.data = trials.data(:,i_chnls,:);     % channel selection
%     trials_plot.text = 'ERP';
%     trials_plot.outDir = [params.storage.dir_results filesep params.storage.outName filesep 'singleTrials' filesep subjTag];
%     plotSingleTrials(params, trials_plot);
    
    %% no trials rejection
    i_tr_ok = [1:size(trials.data,3)];
    
    %% trials filtering
    params.triggering.trialsFiltering = {'z-score'};     % , 'flipInTime', , 'derivative'
%     params.bp_freq.freqBand = params.connectivity.selectedChnls{1,5};
    trials = trialsFiltering(params, trials);
    
    %% compute FC: DTF -> M matrix (saved) = 5D: M(ch_trg,ch_src,freq,time,clz)
%     FC_fileName = [FC_outDir filesep 'M_DTF_mvar' num2str(params.connectivity.MVAR_order) '_FB0_subj' num2str(subj) '_' subjTag '.mat'];
    FC_fileName = [FC_outDir filesep 'M_PLV_' subjTag '.mat'];
    trials_plv = trials;
    trials_plv.data = trials.data(:,i_chnls,i_tr_ok);     % channel selection
%     trials_plv.data = single(trials_plv.data);            % single precision (faster?)
    trials_plv.labels = trials.labels(i_tr_ok);
    trials_plv.selCh_all = i_chnls;
    trials_plv.info.selCh_groups = chGroups;
%     connectivityAnalysis(params, trials_dtf, FC_fileName);
    msSEI_FC_PLV_impl(params, trials_plv, FC_fileName);
    
end % of subj
disp('PLV done.');

