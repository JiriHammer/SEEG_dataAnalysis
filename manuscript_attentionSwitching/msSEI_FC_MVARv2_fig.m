%% FC figures
% job = msSEI_FC_MVARv2_job.m

% (c) Jiri, Dec23

% %% outside loop
% outDirs = {
%     'FC_v10_outputFiles_mvar41', 41;...
%     'FC_v10_outputFiles_mvar119', 119; ...
%     'FC_v11_outputFiles_RJ', 17; ...
% };
% for zsc = 1:2
%     if zsc == 1
%         do_zscore = true;
%     else
%         do_zscore = false;
%     end
%     for d = 1:size(outDirs,1)
%         FC_outDirName = outDirs{d,1};   
%         MVAR_order = outDirs{d,2};    
%         for k = 1:2
%             if k == 1
%                 n_fc = 1;
%                 fcVars = {'PDC','COH','pCOH','DTF','dDTF','ffDTF','GGC'};
%             else
%                 n_fc = 2;
%                 fcVars = {'S','h','pCOH2', 'PDCF', 'coh','Af','GPDC'};
%             end
% 
%             msSEI_FC_MVARv2_fig;
%             clear EC_all M_all;
%         end
%     end
% end

%% select FC vars to plot
if ~exist('fcVars','var')
    n_fc = 1;
    fcVars = {'PDC','COH','pCOH','DTF','dDTF','ffDTF','GGC'};
    % fcVars = {'S','h','pCOH2', 'PDCF', 'coh','Af','GPDC'};
end

% if ~exist('fcVars','var')
%     n_fc = 2;
%     % fcVars = {'PDC','COH','pCOH','DTF','dDTF','ffDTF','GGC'};
%     fcVars = {'S','h','pCOH2', 'PDCF', 'coh','Af','GPDC'};
% end

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

%% uncomment for outside job
% params.connectivity.MVAR_order = MVAR_order;

%% comment-out for outside job
do_zscore = false;

% params.connectivity = struct(...
%     'computeConnectivity', true, ...
%     'connectivityMeasure', {fcVars}, ...   % choice: 'GGC','DTF','ffDTF','PDC','COH'
%     'MVAR_order', 21, ...
%     'subBands', [1 120], ... settings by RJ: [2 12;13 25;26 48;52 75; 76 98; 102 148];
%     'windowSize', 0.5, ...      % in [s]
%     'timeStep', 0.125, ...      % in [s], corresponds to N samples at 512 Hz (if N = 64, tStep = 1/(512/64) = 0.125)
%     'selectedChnls', 'see below', ... % see below: params.coherence.selectedChnls
%     'maxDistToROI', 10, ...      % in [mm], maximum distance (tolerance) to region of interest
%     'plotEC_ch2ch', true, ...       % plot coherence: each channel to all channels
%     'plotEC_selGroups', true, ...   % plot coherence: selected groups of channels (each group to each group)
%     'plotEC_zscore', true, ...       % apply z-score in figs over time dimension
%     'time2plot', [-3,3] ...     % in [s]
%     );    

%% ================== settings: connectivity -> params ====================
params.connectivity.connectivityMeasure = fcVars;

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

% --- v11: MVAR_ord = 20, DTF as RJ (band-pass filtered) in sel_FB = 'loFreqs','loGamma','hiGamma'
FC_outDirName = 'FC_v11_outputFiles_RJ';   
params.connectivity.MVAR_order = 17;
params.connectivity.subBands = [0 30; 31 48; 52, 98; 102, 120];
params.connectivity.selectedChnls = {...    % each row = tuple: atlasName, areaName (must match names in anatomicalAreas_getList.m) 
    'Yeo7',  'Default',           'was_sgnf', {'loFreq0','loGamma','hiGamma'}; ...
    'Yeo7',  'Dorsal Attention',  'was_sgnf', {'loFreq0','loGamma','hiGamma'}; ...  
    }; 

%% output directory
% FC_outDir = [params.storage.dir_results filesep params.storage.outName filesep 'DTF_vRJ_outputFiles'];
FC_outDir = [params.storage.dir_results filesep params.storage.outName filesep FC_outDirName];
if ~exist(FC_outDir, 'dir')
    mkdir(FC_outDir);
end    
    
%% ============================= LOAD DTF on raw data =================================
was_init_M = false;
i_subj = [];
mvar_ord = [];
thr_nChPerGroup = 5;    % only if each group has at least thr_nChPerGroup channels 
EC_all = struct;
for subj = 1:size(params.storage.subjList,1)
    subjTag = params.storage.subjList{subj,1};
    disp(['FC: loading data from subj = ' subjTag ' ...']);
    
    % FC file name -> FC_fileName
    FC_outDir = [params.storage.dir_results filesep params.storage.outName filesep FC_outDirName];
    assert(exist(FC_outDir,'dir')==7);
    FC_fileName = [FC_outDir filesep 'M_mvar' num2str(params.connectivity.MVAR_order) '_FB0_' subjTag '.mat'];
    
    % load & plot FC (DTF/ICOH) matrix
    if exist(FC_fileName,'file') == 2
        
        %% load EC struct -> EC, EC_info
        clear EC EC_info;
        load(FC_fileName, 'EC', 'EC_info');            

        %% cat MVAR order
        if isfield(EC_info, 'mvar_ord')
            mvar_ord = cat(1, mvar_ord, unique(EC_info.mvar_ord(:)));
        end
        
        %% plot FC -> M_subj = 5D: grp_trg x grp_src x freq x time x clz x groups
        EC_info.subjTag = subjTag;
%         M_subj = plotFC_spectral_2NN(M, EC_info, EC_info.info.selCh_groups, thr_nChPerGroup);
        
        %% figure
        nGroups = size(EC_info.info.selCh_groups.ch_inds_trials,1);
        nClz = size(EC_info.info.nClz,2);         
        
        f = fig_make;
        EC_info.nCols = size(params.connectivity.connectivityMeasure,2);
        EC_info.nRows = nGroups*nClz;
        
        % plot each connectivity measure (fcVar)
        for fc = 1:size(params.connectivity.connectivityMeasure,2)
            fcName = params.connectivity.connectivityMeasure{fc};
            EC_info.fcName = fcName;
            
            % FC matrix -> M
            M = EC.(fcName);    % M = 5D matrix: M(ch_trg,ch_src,freq,time,clz)
            assert(size(M,5) == size(EC_info.info.nClz,2));
            
            % zscore over time?
            if do_zscore
                M = zscore(M, 0, 4);    % over 4th dim = time !
            end
        
            % plot
            EC_info.thisCol = fc;
            M_subj = plotFC2axes_5DM_spectral(M, EC_info, EC_info.info.selCh_groups, thr_nChPerGroup);
            
            % init & allocate EC_all = 6D: [aa_src, aa_trg, freq, time, clz, subj]
            if ~isfield(EC_all, fcName)
                EC_all.(fcName) = nan(size(M_subj,1),size(M_subj,2),size(M_subj,3),size(M_subj,4),size(M_subj,5),size(params.storage.subjList,1));
                if ~isempty(M_subj)
                    i_subj = cat(1,i_subj, subj);
                end
            end

            % cat to EC_all
            if ~isempty(M_subj)
                EC_all.(fcName)(:,:,:,:,:,subj) = M_subj; % M_all = 6D: [aa_src, aa_trg, freq, time, clz, subj]
            end            
        end
        
        %% save
        EC_info.outDir = [FC_outDir filesep 'fig_zscore' num2str(do_zscore)];
        EC_info.figName = ['FC_vars' num2str(n_fc) '_subj_' EC_info.subjTag];        
        fig_save(f, EC_info.figName, EC_info.outDir, 'format','png');
        close(f);             
        
    end
end

%% ====================== FIG - ALL: spectral FC ==========================
EC_info.subjTag = 'ALL';
chGroups = EC_info.info.selCh_groups;
for nn = 1:size(chGroups.ch_inds_trials,1)
    chGroups_NN.ch_inds_trials{nn,1} = nn;
    chGroups_NN.nn_name{nn,1} = chGroups.nn_name{nn};
    chGroups_NN.ch_names{nn,1} = params.storage.subjList(i_subj)';
end

% figure
f = fig_make;
EC_info.nCols = size(params.connectivity.connectivityMeasure,2);
EC_info.nRows = nGroups*nClz;

% plot each connectivity measure (fcVar)
for fc = 1:size(params.connectivity.connectivityMeasure,2)
    fcName = params.connectivity.connectivityMeasure{fc};
    EC_info.fcName = fcName;

    % FC matrix -> M
    M_all = EC_all.(fcName);    % M = 6D: [aa_src, aa_trg, freq, time, clz, subj]
    M = nanmean(M_all,6);
    assert(size(M,5) == size(EC_info.info.nClz,2));

    % plot
    EC_info.thisCol = fc;
    M_subj = plotFC2axes_5DM_spectral(M, EC_info, chGroups_NN);
end

% save
EC_info.outDir = [FC_outDir filesep 'fig_zscore' num2str(do_zscore)];
EC_info.figName = ['FC_vars' num2str(n_fc) '_subj_' EC_info.subjTag];        
fig_save(f, EC_info.figName, EC_info.outDir, 'format','png', 'res',600);
close(f);             

%% ============= FIG - plotband: mean over freq & mean over subj ==========
if size(params.connectivity.selectedChnls,2) < 5
    selFreq = [0, 120];
elseif size(params.connectivity.selectedChnls,2) == 5
    selFreq = params.connectivity.selectedChnls{1,5};
end
FB_selFreq = {...
    'wholeBand', [0, 120];
    'delta',[0, 2];
    'theta',[6, 7];
    'alpha',[8, 13];
    'beta',[14, 29];
    'loFreq0',[0, 29];
    'loFreqs',[4, 29];
    'loGamma',[30, 49];
    'hiGamma',[50, 120];
    };
selTime = [-2.1, 2.1];
i_f = closestval(EC_info.freq, selFreq(1)):closestval(EC_info.freq, selFreq(2));
i_t = closestval(EC_info.time, selTime(1)):closestval(EC_info.time, selTime(2));
EC_info.time = EC_info.time(i_t);
EC_info.P_level = P_level;
nGroups = size(M_all,1);
nClz = size(M_all,5);

figName = EC_info.figName;
for fb = 1:size(FB_selFreq,1)
    % selected freq. bins
    selFreq = FB_selFreq{fb,2};
    i_f = closestval(EC_info.freq, selFreq(1)):closestval(EC_info.freq, selFreq(2));
    
    % figure
    f = fig_make;
    nCols = size(params.connectivity.connectivityMeasure,2);
    nRows = 3;  % g1->g2 & g2->g1 & (g1->g2) - (g2->g1)

    marg_h = [0.08 0.08];   % margin from [bottom, top]
    marg_w = [0.03 0.03];   % margin from [L, R] side
    gap = [0.02, 0.02];     % between axes from [top, side]

    % plot all FC vars
    for fc = 1:size(params.connectivity.connectivityMeasure,2)
        fcName = params.connectivity.connectivityMeasure{fc};
        EC_info.fcName = fcName;
        nPlot = 1;

        % FC matrix -> M
        M_all = EC_all.(fcName);    % M = 6D: [aa_src, aa_trg, freq, time, clz, subj]

        % -----------A: plotband DTF: g1 -> g2 (src -> trg) --------------------------
        src = 1;
        trg = 2;
        FC = squeeze(nanmean(M_all(trg,src,i_f,i_t,:,:),3));   % 2D: t x clz, mean over freq

        % axes
        thisPlot = fc + (nPlot-1)*nCols;   % plot clz below each other
        ax = subtightplot(nRows, nCols, thisPlot, gap, marg_h, marg_w); 
        hold on;

        str_title_g1_g2 = plotFC_getTitle(chGroups_NN, src, trg, [], EC_info.fcName);
        EC_info.str_title = str_title_g1_g2;

        plotEC_plotband_2NN(FC, EC_info);
        nPlot = nPlot+1;

        % -----------B: plotband DTF: g2 -> g1 (src -> trg) --------------------------
        src = 2;
        trg = 1;
        FC = squeeze(nanmean(M_all(trg,src,i_f,i_t,:,:),3));   % 2D: t x clz, mean over freq

        % axes
        thisPlot = fc + (nPlot-1)*nCols;   % plot clz below each other
        ax = subtightplot(nRows, nCols, thisPlot, gap, marg_h, marg_w); 
        hold on;

        str_title_g2_g1 = plotFC_getTitle(chGroups_NN, src, trg, [], EC_info.fcName);
        EC_info.str_title = str_title_g2_g1;

        plotEC_plotband_2NN(FC, EC_info);
        nPlot = nPlot+1;

        % -----------C: plotband DTF: (g1 -> g2) - (g2 -> g1) --------------------------
        src = 1;
        trg = 2;
        M_g1_g2 = nanmean(M_all(trg,src,i_f,i_t,:,:),3);  % mean over freq, 6D: [1,1,1,t,clz,subj]
        M_g2_g1 = nanmean(M_all(src,trg,i_f,i_t,:,:),3);  % mean over freq, 6D: [1,1,1,t,clz,subj]
        FC = squeeze(M_g1_g2 - M_g2_g1);           % 3D: t x clz x subj
        % FC_sem = squeeze(sem(M_g1_g2 - M_g2_g1,6));     % 2D: t x clz, SEM over subj

        % axes
        thisPlot = fc + (nPlot-1)*nCols;   % plot clz below each other
        ax = subtightplot(nRows, nCols, thisPlot, gap, marg_h, marg_w); 
        hold on;

        EC_info.str_title = [str_title_g1_g2 ' - ' str_title_g2_g1];

        plotEC_plotband_2NN(FC, EC_info);
        nPlot = nPlot+1;
    end

    % save
%     EC_info.figName = [EC_info.figName '_plotband'];
    EC_info.figName = [figName '_plotband_' num2str(fb) '_' FB_selFreq{fb,1}];
    fig_save(f, EC_info.figName, EC_info.outDir, 'format','png', 'res',600);
    close(f);   
end

