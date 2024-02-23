%% ms SEI: job - functional connectivity (FC) between 2 selected AAs (DMN & DAN)
% figs = msSEI_FC_DTF_fig_v2_selCh.m

% (c) Jiri, Oct23

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
runner.newFDR = false;  % true = updates the FDR correction over all patients (and freq. bands, ...) 

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
    'connectivityMeasure', {{'S','h','PDC','COH','DTF','pCOH','dDTF','ffDTF', 'pCOH2', 'PDCF', 'coh','GGC','Af','GPDC'}}, ...   % choice: 'GGC','DTF','ffDTF','PDC','COH'
    'MVAR_order', 21, ...
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

% --- v9: MVAR_ord = AIC, sel_FB = 'loFreqs','loGamma','hiGamma'
% FC_outDirName = 'DTF_v9_outputFiles_LFB_HGB';
% FC_outDirName = 'DTF_v9_wiZscore';   
% params.connectivity.MVAR_order = 'AIC';
% params.connectivity.selectedChnls = {...    % each row = tuple: atlasName, areaName (must match names in anatomicalAreas_getList.m) 
%     'Yeo7',  'Default',           'was_sgnf', {'loFreq0','loGamma','hiGamma'}; ...
%     'Yeo7',  'Dorsal Attention',  'was_sgnf', {'loFreq0','loGamma','hiGamma'}; ...  
%     }; 

% --- v10: MVAR_ord = 40, sel_FB = 'loFreqs','loGamma','hiGamma'
% FC_outDirName = 'FC_v10_outputFiles_mvar119';   
% params.connectivity.MVAR_order = 119;
% params.connectivity.selectedChnls = {...    % each row = tuple: atlasName, areaName (must match names in anatomicalAreas_getList.m) 
%     'Yeo7',  'Default',           'was_sgnf', {'loFreq0','loGamma','hiGamma'}; ...
%     'Yeo7',  'Dorsal Attention',  'was_sgnf', {'loFreq0','loGamma','hiGamma'}; ...  
%     }; 

% --- v11: MVAR_ord = 17, DTF as RJ (band-pass filtered) in sel_FB = 'loFreq0','loGamma','hiGamma'
% FC_outDirName = 'FC_v11_outputFiles_RJ';   
% params.connectivity.MVAR_order = 17;
% params.connectivity.subBands = [0 30; 31 48; 52, 98; 102, 120];
% params.connectivity.selectedChnls = {...    % each row = tuple: atlasName, areaName (must match names in anatomicalAreas_getList.m) 
%     'Yeo7',  'Default',           'was_sgnf', {'loFreq0','loGamma','hiGamma'}; ...
%     'Yeo7',  'Dorsal Attention',  'was_sgnf', {'loFreq0','loGamma','hiGamma'}; ...  
%     }; 

% --- v12: MVAR_ord = 45, DTF as RJ (band-pass filtered) in sel_FB = 'loFreq0','loGamma','hiGamma'
% FC_outDirName = 'FC_v12_ord45_RJ';   
% params.connectivity.MVAR_order = 45;
% params.connectivity.subBands = [0 30; 31 48; 52, 98; 102, 120];
% params.connectivity.selectedChnls = {...    % each row = tuple: atlasName, areaName (must match names in anatomicalAreas_getList.m) 
%     'Yeo7',  'Default',           'was_sgnf', {'loFreq0','loGamma','hiGamma'}; ...
%     'Yeo7',  'Dorsal Attention',  'was_sgnf', {'loFreq0','loGamma','hiGamma'}; ...  
%     }; 

% --- v12: MVAR_ord = 45, DTF as RJ (band-pass filtered) in sel_FB = 'loFreq0','loGamma','hiGamma'
FC_outDirName = 'FC_v12_ord45_FB7_RJ';   
params.connectivity.MVAR_order = 45;
params.connectivity.subBands = [0 30; 31 48; 52, 98; 102, 120];
params.connectivity.selectedChnls = {...    % each row = tuple: atlasName, areaName (must match names in anatomicalAreas_getList.m) 
    'Yeo7',  'Default',           'was_sgnf', {'delta','theta','alpha','beta','loGamma','hiGamma'}; ...
    'Yeo7',  'Dorsal Attention',  'was_sgnf', {'delta','theta','alpha','beta','loGamma','hiGamma'}; ...  
    }; 

%% output directory
% FC_outDir = [params.storage.dir_results filesep params.storage.outName filesep 'DTF_vRJ_outputFiles'];
FC_outDir = [params.storage.dir_results filesep params.storage.outName filesep FC_outDirName];
if ~exist(FC_outDir, 'dir')
    mkdir(FC_outDir);
end    
    
%% compute DTF for each subj
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
    FC_fileName = [FC_outDir filesep 'M_mvar' num2str(params.connectivity.MVAR_order) '_FB0_' subjTag '.mat'];
    trials_dtf = trials;
    trials_dtf.data = trials.data(:,i_chnls,i_tr_ok);     % channel selection
    trials_dtf.labels = trials.labels(i_tr_ok);
    trials_dtf.selCh_all = i_chnls;
    trials_dtf.info.selCh_groups = chGroups;
    connectivityAnalysis(params, trials_dtf, FC_fileName);
        
end % of subj
disp('FC done.');

return;




%% old, MSU? (now=12.10.2023) - selects channels from 2 NNs

%% settigs: path to exported trials for FC
FC_pathBeg = 'G:\dox\ms_switch_EI\data\v18_fc_trials\switchin_EI_IE_car';
FC_outDir = [FC_pathBeg filesep 'DTF_outputFiles'];
if ~exist(FC_outDir, 'dir')
    mkdir(FC_outDir);
end    
    
%% selected NN
FC_roi = {...
    'Default',          'Yeo7';
    'Dorsal Attention', 'Yeo7';
    };
FC_freqBands = {...
    'delta',    [0.5, 3]; ...
    'theta',    [4, 7]; ...
    'alpha',    [8, 12]; ...
    'beta',     [13 20]; ...
    'beta',     [21, 30]; ...
    'loGamma',  [31, 40]; ...
    'loGamma',  [41, 50]; ...
    };
hifbs = cat(2, 1+[50:10:110]', 10+[50:10:110]');
for n = 1:size(hifbs,1)
    FC_freqBands{end+1,1} = 'hiGamma';
    FC_freqBands{end,2} = hifbs(n,:);
end

%% =============================LOAD DATA==================================
%% ch2roi for Yeo7 -> GRUP 
if load_ch2roi
    clear GRUP;
    groupInfo = struct(...
        'brainAtlas', {{brainAtlas}}, ... % brain atlases: 'Yeo7','Yeo7_addLat','Yeo17','Yeo17_addLat','Mars','Mars_addLat'
        'freqBands', {{'delta','theta','alpha','beta','loGamma','hiGamma'}}, ...  % freq. bands (FB)
        'signifOfActiv', {{'was_sgnf'}}, ...   % significance of activation
        'significanceLevel', P_level, ...      % significance level for ROI (e.g. DMN)
        'significanceTimeThr', t_siglen, ...    % in [s], significance time threshold (ch is signif. if its significance > threshold) 
        'GW_selection', {{'any'}}, ...      % grey / white matter: choices: 'gm','wm','csf','bone', 'any'
        'GW_threshold', 60, ...    % probability, in [%]
        'plotStyle', {{'allFB_allROI'}}, ... 'oneFB_allROI','allFB_allROI','spectra_SNR'
        'maxDistToROI', 10, ...   % in [mm], maximum distance (tolerance) to region of interest
        'minSubjInROI', 1, ...  % minimum of subjects (patients) in a ROI (if a ROI has less patients, it is not considered in plots)
        'subtractAvgResponse', false, ...    % subtracts mean response over classes (conditions) from each clz
        'dirSuffix', dirSuffix, ... % adds suffix to 'ch2roi' directory name -> ch2roi_outDirSuffix
        'figSuffix', '', ... % adds suffix to figure name
        'time2plot', t_sel ... % in [s], w.r.t. cutting point (= 0)
        );
    ch2roi_wrapper;
end

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
params.connectivity.selectedChnls = {...    % each row = tuple: atlasName, areaName (must match names in anatomicalAreas_getList.m) 
    'all','all'; ...
    };

%% compute DTF for each subj
for subj = 1:size(params.storage.subjList,1)
    subjTag = params.storage.subjList{subj,1};
    disp(['FC - DTF: processing, subj = ' subjTag]);
    
    % --- cache file: FB activation
    cacheFile = [params.storage.dir_results filesep params.storage.outName filesep subjTag filesep 'cacheFile.mat'];
    assert(exist(cacheFile,'file') == 2);
    clear H selCh_H_resp;
    load(cacheFile, 'H', 'selCh_H_resp');   
    selCh_H_FBA = selCh_H_resp;         

    % --- cache file: exported trials for FC 
    cacheFile = [FC_pathBeg filesep subjTag filesep 'cacheFile.mat'];
    assert(exist(cacheFile,'file') == 2);        
    clear H selCh_H_resp;
    load(cacheFile, 'H', 'selCh_H_resp');   
    selCh_H_FC = selCh_H_resp;       
    assert(size(selCh_H_FC,2) == size(selCh_H_FBA,2));  % !!! MUST be the same !!!
    params.storage.cacheFile = cacheFile;               % used in connectivityAnalysis.m
    
    % load raw trials
    clear trials
    load(cacheFile, 'trials');  
    
    % -----------DTF: FB = 0 -> DTF on raw data -> M matrix (saved)--------
    % channel selection (based on alpha)
    [tf, i_freq] = ismember('alpha', groupInfo.list_FB(:,1)); assert(tf);
    i_chnls = [];
    FC_chnls = cell(1, size(FC_roi,1));
    chNotFound = false;
    for aa = 1:size(FC_roi,1)
        [tf, i_aa] = ismember(FC_roi{aa}, groupInfo.list_AA(:,1)); assert(tf);
        i_subj = find(GRUP.ind_chSubj{i_aa,i_freq,1}(:,2) == subj); % select subj ID
        if isempty(i_subj)
            disp(['FC - DTF: no channels found, subj = ' subjTag ', AA = ' FC_roi{aa}]);
            chNotFound = true;
            break;
        end
        FC_chnls{aa} = GRUP.ind_chSubj{i_aa,i_freq,1}(i_subj,1);     % select channels of this subj
        i_chnls = cat(1, i_chnls, FC_chnls{aa});
    end
    if chNotFound
        continue;
    end
        
    % compute FC: DTF 
    FC_fileName = [FC_outDir filesep 'M_DTF_FB0_subj' num2str(subj) '_' subjTag '.mat'];
    trials_dtf = trials;
    trials_dtf.data = trials.data(:,i_chnls,:);     % channel selection
    trials_dtf.selCh_all = i_chnls;
    trials_dtf.selCh_groups = FC_chnls;
    connectivityAnalysis(params, trials_dtf, FC_fileName);
    
    % -----------DTF: for FB -> M matrix (saved)--------
    for freq = 1:size(FC_freqBands,1)
        disp([' - computing DTF: freq band = ' num2str(freq) '/' num2str(size(FC_freqBands,1))]);
        
        % channel selection (based on freq. sub-band)
        [tf, i_freq] = ismember(FC_freqBands{freq,1}, groupInfo.list_FB(:,1)); assert(tf);
        i_chnls = [];
        FC_chnls = cell(1, size(FC_roi,1));
        chNotFound = false;
        for aa = 1:size(FC_roi,1)
            [tf, i_aa] = ismember(FC_roi{aa}, groupInfo.list_AA(:,1)); assert(tf);
            i_subj = find(GRUP.ind_chSubj{i_aa,i_freq,1}(:,2) == subj); % select subj ID
            if isempty(i_subj)
                disp(['FC - DTF: no channels found, subj = ' subjTag ', AA = ' FC_roi{aa}]);
                chNotFound = true;
                break;
            end
            FC_chnls{aa} = GRUP.ind_chSubj{i_aa,i_freq,1}(i_subj,1);     % select channels of this subj
            i_chnls = cat(1, i_chnls, FC_chnls{aa});
        end        
        if chNotFound
            break;
        end
        
        % channel selection
        trials_dtf = trials;
        trials_dtf.data = trials.data(:,i_chnls,:);     % channel selection
        trials_dtf.selCh_all = i_chnls;
        trials_dtf.selCh_groups = FC_chnls;        
        
        % filter data: bp, HT, ampEnv -> trials_filt
        params.triggering.trialsFiltering = {'bandPass', 'hilbAmp'};
        params.bp_freq.freqBand = FC_freqBands{freq,2};
        trials_dtf = trialsFiltering(params, trials_dtf);
        
        % DTF: FB = freq -> DTF on raw data -> M matrix (saved)
%         params.connectivity.MVAR_order = 10;
        params.connectivity.subBands = [1, FC_freqBands{freq,2}(2)];    % takes all freq until cutoff (spectrum of hilbAmpEnv can be rich also below the low-freq cutoff)
        FC_fileName = [FC_outDir filesep 'M_DTF_FB' num2str(freq) '_subj' num2str(subj) '_' subjTag '.mat'];
        connectivityAnalysis(params, trials_dtf, FC_fileName);
        
    end % of freq
        
end % of subj
