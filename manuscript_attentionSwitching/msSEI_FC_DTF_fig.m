%% ms SEI: DTF fig - connectivity between 2 selected AAs

% (c) Jiri, Oct23

%% settings: v1 - relaxed conditions (p-level, length of significance, ...)
FC_outputFiles = 'DTF_v2_outputFiles';
% FC_outputFiles = 'DTF_outputFiles';
% FC_outputFiles = 'DTF_vRJ_outputFiles';
t_sel = [-1.0, 2.0];    % the true sel time will be cut by 1 s on edges (= template length): [-2, 2] -> [-1, 1]
P_level = 0.05;         % significance level
t_siglen = 0.10;        % in [s], significance time threshold (ch is signif. if its significance > threshold) 

%% get job-related settings (params): jobExecutor_SEI.m, v11 -> params
% --- specify analysis to be performed -> runner
runner = struct;
runner.newJobs = false; % set to 'false' if the results are already computed
runner.newFDR = false;  % true = updates the FDR correction over all patients (and freq. bands, ...) 

% --- job settings: specific (define spectral method) -> params
spectralMethod = 'v24_stft_baseRS_car';         % ~ v24 (CAR), STFT with RS baseline
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
% params.connectivity.selectedChnls = {...    % each row = tuple: atlasName, areaName (must match names in anatomicalAreas_getList.m) 
%     'all','all'; ...
%     };
% params.connectivity.selectedChnls = {...    % each row = tuple: atlasName, areaName (must match names in anatomicalAreas_getList.m) 
%     'Yeo7',  'Default',           'was_sgnf', {'delta','theta','alpha','beta','loGamma','hiGamma'}; ...
%     'Yeo7',  'Dorsal Attention',  'was_sgnf', {'delta','theta','alpha','beta','loGamma','hiGamma'}; ...  
%     'Yeo17', 'Salience',          'was_sgnf', {'delta','theta','alpha','beta','loGamma','hiGamma'}; ...  
%     };  

% figDirSuffix = 'onlyAlpha';
% params.connectivity.selectedChnls = {...    % each row = tuple: atlasName, areaName (must match names in anatomicalAreas_getList.m) 
%     'Yeo7',  'Default',           'was_sgnf', {'delta','theta','alpha','beta','loGamma','hiGamma'}; ...
%     'Yeo7',  'Dorsal Attention',  'was_sgnf', {'delta','theta','alpha','beta','loGamma','hiGamma'}; ...  
%     'Yeo17', 'Salience',          'was_sgnf', {'delta','theta','alpha','beta','loGamma','hiGamma'}; ...  
%     };  

% figDirSuffix = 'onlyAlpha';
% params.connectivity.selectedChnls = {...    % each row = tuple: atlasName, areaName (must match names in anatomicalAreas_getList.m) 
%     'Yeo7',  'Default',           'was_sgnf', {'alpha'}; ...
%     'Yeo7',  'Dorsal Attention',  'was_sgnf', {'alpha'}; ...  
%     'Yeo17', 'Salience',          'was_sgnf', {'alpha'}; ...  
%     };  

figDirSuffix = 'onlyAlpha_DMN_DAN';
params.connectivity.selectedChnls = {...    % each row = tuple: atlasName, areaName (must match names in anatomicalAreas_getList.m) 
    'Yeo7',  'Default',           'was_sgnf', {'alpha'}; ...
    'Yeo7',  'Dorsal Attention',  'was_sgnf', {'alpha'}; ...  
    }; 

%% ============================= LOAD DTF on raw data =================================
do_zscore = true;
was_init_M = false;
for subj = 1:size(params.storage.subjList,1)
    subjTag = params.storage.subjList{subj,1};
    
    % get selected channels -> chGroups
    params.P_level = P_level;
    params.t_siglen = t_siglen;
    chGroups = ch2roi_getSelChnls(params, subjTag);

    % load & plot DTF matrix
    FC_outDir = [params.storage.dir_results filesep params.storage.outName filesep FC_outputFiles];
    assert(exist(FC_outDir,'dir')==7);
    FC_fileName = [FC_outDir filesep 'M_DTF_mvar' num2str(params.connectivity.MVAR_order) '_FB0_subj' num2str(subj) '_' subjTag '.mat'];
    if exist(FC_fileName,'file') == 2
        
        %% load DTF matrix -> 5D matrix: M(ch_trg,ch_src,freq,time,clz)
        clear M EC_info;
        load(FC_fileName, 'M', 'EC_info');
        
        %% zscore over time
        if do_zscore
            M = zscore(M, 0, 4);    % over 4th dim = time !
        end
        
        %% plot FC -> M_subj = 5D: grp_trg x grp_src x freq x time x clz x groups
        EC_info.subjTag = subjTag;
        EC_info.outDir = [FC_outDir filesep 'DTF_fig_zscore' num2str(do_zscore) '_' figDirSuffix];
        EC_info.figName = ['DTF_raw_subj_' EC_info.subjTag];
        M_subj = plotFC_spectral(M, EC_info, chGroups);
        
        %% cat to all
        if ~was_init_M
            M_all = nan(size(M_subj,1),size(M_subj,2),size(M_subj,3),size(M_subj,4),size(M_subj,5),size(params.storage.subjList,1));
            was_init_M = true;
        end
        M_all(:,:,:,:,:,subj) = M_subj;
    end
end

%% plot: mean over subj
M = nanmean(M_all,6);
EC_info.subjTag = 'ALL';
EC_info.outDir = [FC_outDir filesep 'DTF_fig_zscore' num2str(do_zscore) '_' figDirSuffix];
EC_info.figName = ['DTF_raw_subj_' EC_info.subjTag];
for nn = 1:size(chGroups.ch_inds_trials,1)
    chGroups_NN.ch_inds_trials{nn,1} = nn;
    chGroups_NN.nn_name{nn,1} = chGroups.nn_name{nn};
    chGroups_NN.ch_names{nn,1} = params.storage.subjList';
end
M_subj = plotFC_spectral(M, EC_info, chGroups_NN);
        










return;

%% old, MSU? (now=Oct23)
%% selected NN
FC_roi = {...
    'Default',          'DMN', 'Yeo7';
    'Dorsal Attention', 'DAN', 'Yeo7';
    };
FC_roiInteraction = {...
    'DMN -> DAN';
    'DAN -> DMN';
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

%% settigs: path to exported trials for FC
FC_pathBeg = 'G:\dox\ms_switch_EI\data\v18_fc_trials\switchin_EI_IE_car';
FC_outDir = [FC_pathBeg filesep 'DTF_outputFiles'];
if ~exist(FC_outDir, 'dir')
    mkdir(FC_outDir);
end    
    
%% settings: v1 - relaxed conditions (p-level, length of significance, ...)
dirSuffix = '_msFig_AA_v1_relaxCond';
t_sel = [-1.0, 2.0];    % the true sel time will be cut by 1 s on edges (= template length): [-2, 2] -> [-1, 1]
P_level = 0.05;         % significance level
t_siglen = 0.10;        % in [s], significance time threshold (ch is signif. if its significance > threshold) 
brainAtlas = 'Yeo7';

%% settings: v3
% figureSuffix = '_v3_strictCond';
% t_sel = [-1.0, 2.0];    % the true sel time will be cut by 1 s on edges (= template length): [-2, 2] -> [-1, 1]
% P_level = 0.01;         % significance level
% t_siglen = 0.25;        % in [s], significance time threshold (ch is signif. if its significance > threshold) 

%% settings
load_ch2roi = false;             % load again the GRUP struct ?

%% run: jobExecutor_SEI.m, v11 -> params
% specify analysis to be performed -> runner
runner = struct;
runner.newJobs = false; % set to 'false' if the results are already computed
runner.newFDR = false;  % true = updates the FDR correction over all patients (and freq. bands, ...) 

% job settings: common (data preprocessing & analysis) -> params
params = SEI_job_settingsCommon();

% job settings: specific (define spectral method) -> params
% spectralMethod = 'stft_sessions';       % ~ v12
spectralMethod = 'stft_baseRS';         % ~ v11  (CAR)
% spectralMethod = 'stft_baseRS_bip';     % ~ v11b (BIP)
params = SEI_job_spectralMethod(params, spectralMethod);

% job settings: define the classes for the triggers -> params
clzDefName = 'switchin_EI_IE';
params = SEI_job_classDefinitions(params, clzDefName);

% >>>>>>>>>>>>>>>>>>>>>>>> run jobs <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
if runner.newJobs
    for subj = 1:size(params.storage.subjList,1)
        main_execute(params, subj);
%         main_execute;
    end
    disp('Done.');
    disp(['Time: ' datestr(now, 'dd-mmm-yyyy HH:MM:SS')]);
end
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

%% FDR correction
params.storage.subjList = subjList_updateAnalyzed(params);
if runner.newFDR
    significanceLevel = P_level;
    fields_signifVals = {'pVals','hVals'};
    poolSubjects_fdrCorrection(params, significanceLevel, fields_signifVals);
end

%% ============================= LOAD DTF on raw data =================================
M_all = [];
do_zscore = true;
for subj = 1:size(params.storage.subjList,1)
    subjTag = params.storage.subjList{subj,1};
    
    % load & plot DTF matrix
    FC_fileName = [FC_outDir filesep 'M_DTF_FB0_subj' num2str(subj) '_' subjTag '.mat'];
    if exist(FC_fileName,'file') == 2
        
        %% load DTF matrix
        clear M EC_info;
        load(FC_fileName, 'M', 'EC_info');
        nGroups = size(EC_info.selCh_groups,2);
        nClz = size(EC_info.info.clzNames,1);
        i_ch_g1 = EC_info.selCh_groups{1};
        i_ch_g2 = EC_info.selCh_groups{2};
        FC_g1_g2 = squeeze(mean(mean(M(i_ch_g1,i_ch_g2,:,:,:),1),2));   % 3D: freq x t x clz
        FC_g2_g1 = squeeze(mean(mean(M(i_ch_g2,i_ch_g1,:,:,:),1),2));   % 3D: freq x t x clz
        M_subj = cat(4, FC_g1_g2, FC_g2_g1);                            % 4D: freq x t x clz x grp
        if do_zscore
            M_subj = zscore(M_subj, 0, 2);
            cLims = [0 0.2];
        end
        
        %% interaction labels (hardcoded for i_ch_g1)
        g_src = 1;
        g_trg = 2;
        FC_roiInteraction{1,1} = [FC_roi{g_src,2} '(' num2str(size(i_ch_g1,2)) ') -> ' FC_roi{g_trg,2} '(' num2str(size(i_ch_g2,2)) ')'];
        g_src = 2;
        g_trg = 1;        
        FC_roiInteraction{2,1} = [FC_roi{g_src,2} '(' num2str(size(i_ch_g2,2)) ') -> ' FC_roi{g_trg,2} '(' num2str(size(i_ch_g1,2)) ')'];
        
        %% plot as spectra for 2 channels
        f = fig_make;
        nCols = 2;
        nRows = 2;
        nPlot = 1;
        for clz = 1:nClz
            for grp = 1:nGroups
                % axes
                subplot(nRows, nCols, nPlot);
                hold on;
                
                % plot: imagesc
                h = imagesc(EC_info.time, EC_info.freq, M_subj(:,:,clz,grp));
%                 h = imagesc(EC_info.time, EC_info.freq, M_subj(:,:,clz,grp), cLims);
                colormap(gca,brewermap(256,'*RdBu')); 
                axis tight
                colorbar;
                
                % time t = 0
                plot([0, 0], ylim, '--k');
                
                % labels
                xlabel('time (s)');
                ylabel('freq (Hz)');
                title(['DTF: subj = ' subjTag ', ROI = ' FC_roiInteraction{grp} ', clz = ' EC_info.info.clzNames{clz}]);
                nPlot = nPlot+1;
            end
        end               
        
        % save
        figName = ['DTF_raw_subj_' subjTag];
        outDir = [FC_outDir filesep 'DTF_fig'];
        fig_save(f, figName, outDir, 'format','png', 'res',600);
%         close(f);         
        
        %% cat to all
        M_all = cat(5, M_all, M_subj);
    end
end

% plot all
g_src = 1;
g_trg = 2;
FC_roiInteraction{1,1} = [FC_roi{g_src,2} ' -> ' FC_roi{g_trg,2}];
g_src = 2;
g_trg = 1;        
FC_roiInteraction{2,1} = [FC_roi{g_src,2} ' -> ' FC_roi{g_trg,2}];
M_subj = nanmean(M_all, 5);

% --- plot
f = fig_make;
nCols = 2;
nRows = 2;
nPlot = 1;
for clz = 1:nClz
    for grp = 1:nGroups
        % axes
        subplot(nRows, nCols, nPlot);
        hold on;

        % plot: imagesc
        h = imagesc(EC_info.time, EC_info.freq, M_subj(:,:,clz,grp));
%                 h = imagesc(EC_info.time, EC_info.freq, M_subj(:,:,clz,grp), cLims);
        colormap(gca,brewermap(256,'*RdBu')); 
        axis tight
        colorbar;

        % time t = 0
        plot([0, 0], ylim, '--k');

        % labels
        xlabel('time (s)');
        ylabel('freq (Hz)');
        title(['DTF: subj = ALL, ROI = ' FC_roiInteraction{grp} ', clz = ' EC_info.info.clzNames{clz}]);
        nPlot = nPlot+1;
    end
end               

% save
figName = ['DTF_raw_subj_ALL'];
outDir = [FC_outDir filesep 'DTF_fig'];
fig_save(f, figName, outDir, 'format','png', 'res',600);
%         close(f);         
