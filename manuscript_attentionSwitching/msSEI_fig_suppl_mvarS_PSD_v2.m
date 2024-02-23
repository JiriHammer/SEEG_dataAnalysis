%% TO DO !!! ms SEI suppl. fig: PSD & MVAR spectra (S)
% TO DO: remove the discrepency !!!
% works for two NN: DMN & DAN
% loads EC results (for 17 subjects with at least 5 signif. channels in the FBs)
% loads PSD results for the same subjects & channels
% FIGURE
% cols = FBs (delta, ..., hiGamma)
% rows =
%   - PSD
%   - MVAR S
% jobs: 
%   - MVAR S = msSEI_FC_MVARv2_job.m (filtered in 3 FB: LFB, LGB, HGB)
%   - PSD = jobsExecutor_SEI.m

% (c) Jiri, Jan24
% based on: msSEI_fig_3D_FB_DTF_v4.m

%% settings: relaxed conditions (p-level, length of significance, ...)
dirSuffix = '_msSEI_fig_suppl_mvarS_PSD';
t_sel = [-1.0, 2.0];    % the true sel time will be cut by 1 s on edges (= template length): [-2, 2] -> [-1, 1]
P_level = 0.05;         % significance level
t_siglen = 0.10;        % in [s], significance time threshold (ch is signif. if its significance > threshold) 
brainAtlas = 'Yeo7';

load_ch2roi = true;
load_dtfResults = true;
load_plvResults = true;

plot_fig1_brain3D = true;
plot_fig2_PLV_DTF = true;

%% FC matrix normalization (sensitive parameter!)
% norm_fcMatrix = 'zscore_time';
norm_fcMatrix = 'zscore_time_clz';
% norm_dtfMatrix = 'RS';

%% settings: FB + DTF results (same dir for SEI & RS) -> label_FB
% label_FB = {
%     'loFreqs', [4, 30],   'LFB: 4-30 Hz',   'DTF_v12_outputFiles_LFB';
%     'loGamma', [31, 50],  'LGB: 31-50 Hz',  'DTF_v13_outputFiles_LGB';
%     'hiGamma', [51, 120], 'HGB: 51-120 Hz', 'DTF_v14_outputFiles_HGB';
%     };

% label_FB = {
%     'delta',   [1,3],    '\delta: 0-3 Hz',  '\delta ',  'DTF_v12_7FB_RJ';
%     'theta',   [4,7],    '\theta: 4-7 Hz',  '\theta ',  'DTF_v12_7FB_RJ';
%     'alpha',   [8,12],   '\alpha: 8-12 Hz', '\alpha ',  'DTF_v12_7FB_RJ';
%     'beta',    [13, 30], '\beta: 13-30 Hz', '\beta ',   'DTF_v12_7FB_RJ';
%     'loGamma', [31, 50], 'LGB: 31-50 Hz',   'LGB',      'DTF_v12_7FB_RJ';
%     'hiGamma', [51, 120],'HGB: 51-120 Hz',  'HGB',      'DTF_v12_7FB_RJ';
%     };
% MVAR_order = 20;

label_FB = {
    'delta',   [0,3],    '\delta: 0-3 Hz',  '\delta ',  'FC_v11_outputFiles_RJ';
    'theta',   [4,7],    '\theta: 4-7 Hz',  '\theta ',  'FC_v11_outputFiles_RJ';
    'alpha',   [8,12],   '\alpha: 8-12 Hz', '\alpha ',  'FC_v11_outputFiles_RJ';
    'beta',    [13, 30], '\beta: 13-30 Hz', '\beta ',   'FC_v11_outputFiles_RJ';
    'loGamma', [31, 50], 'LGB: 31-50 Hz',   'LGB',      'FC_v11_outputFiles_RJ';
    'hiGamma', [51, 120],'HGB: 51-120 Hz',  'HGB',      'FC_v11_outputFiles_RJ';
    };

%% get job-related settings -> params
% --- specify analysis to be performed -> runner
runner = struct;
runner.newJobs = false; % set to 'false' if the results are already computed
runner.newFDR = false;  % true = updates the FDR correction over all patients (and freq. bands, ...) 

% --- job settings: specific (define spectral method) -> params
spectralMethod = 'v31_stft_session_bip';           % ~ v31 (BIP), similar to v12 (STFT base = sess, no zscore, no rej, FB = loFreq0)
% spectralMethod = 'v24_stft_baseRS_car';         % ~ v24 (CAR), STFT with RS baseline
% spectralMethod = 'v24_stft_baseRS_bip';         % ~ v24 (BIP), STFT with RS baseline

params = msSEI_getJobParams(runner, spectralMethod, P_level);
SEI_dir_results = params.storage.dir_results;

%% DTF settings 
% --- v11: MVAR_ord = 17, DTF as RJ (band-pass filtered) in sel_FB = 'loFreq0','loGamma','hiGamma'
% FC_outDirName = 'FC_v11_outputFiles_RJ';   
% params.connectivity.MVAR_order = 17;
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

%% settings: AA -> label_AA
label_AA = msSEI_get_label_AA(lower(brainAtlas));

sel_AA = params.connectivity.selectedChnls(:,2);
% sel_AA = {
%     'Dorsal Attention';
%     'Default';
% };
    

%% load DTF results to update subj list -> params.storage.subjList
params.connectivity.normalize = 'zscore_time';
if load_dtfResults
    fc_measure = 'DTF';
    disp(['Loading: ' fc_measure ' ...']);
    [M_all, EC_info] = msSEI_getDTF_results(params, FC_outDirName, fc_measure);
    EC_time = EC_info.time;   % time axis for DTF
    params.storage.subjList = EC_info.subjList; % updates subj (those having more than 3 channels/ROI)
    subjList = EC_info.subjList;
end
disp(['DTF results loaded. nSubj = ' num2str(size(EC_info.subjList,1))]);
    
%% ========================= FIGURE: RBP + MVAR-S ==========================
nCols = size(label_FB,1);   % ~ freq
nRows = 4;  % PSD, MVAR-S: DMN & DAN
nPlot = 1;

fig_W = 24;         % in [cm]
fig_H = 4*nRows;    % in [cm], 3 cm per row
f = fig_make('fig_W',fig_W, 'fig_H',fig_H);

marg_h = [0.10 0.10];   % margin from [bottom, top]
marg_w = [0.08 0.02];   % margin from [L, R] side
gap = [0.02, 0.02];     % between axes from [top, side] 

fontSize = 10;
time2plot = [-2.3, 2.3];

for nn = 1:size(sel_AA,1)
    [aa_found, aa] = ismember(sel_AA{nn}, label_AA(:,1));   % aa = index to GRUP !!!
    assert(aa_found);
    assert(strcmp(label_AA{aa,1}, sel_AA{nn}));        
    
    %% =========================== row 1: PLOT PSD - per FB ==========================
    for fb = 1:size(label_FB,1)
        selFreq = label_FB{fb,2}; 
        disp([' - plotting: ' label_FB{fb,1}]);

        % ------------------------ load: RBP ----------------------------------
        nCh = 0;       
        y = [];
        for subj = 1:size(params.storage.subjList,1)
            subjTag = params.storage.subjList{subj,1};
            
            % load -> H, selCh_H, trialsData
            cacheFile = [params.storage.dir_results filesep params.storage.outName filesep subjTag filesep 'cacheFile.mat'];
            assert(exist(cacheFile,'file') == 2);
            trialsData_name = ['trialsData_' label_FB{fb,1}];
            clear H selCh_H_resp;
            clear(trialsData_name);
            load(cacheFile, 'H', 'selCh_H_resp', trialsData_name);   
            selCh_H = selCh_H_resp; 
            trialsData = eval(trialsData_name);
            assert(size(selCh_H,2) == size(trialsData.yVals,2));
            
            % select channels
            i_ch = [];
            selCh_names = EC_info.ch_names{nn,subj};    % from DTF
            nCh = nCh + size(selCh_names,2);
            for ch1 = 1:size(selCh_names,2)
                for ch2 = 1:size(selCh_H,2)
                    if strcmp(selCh_names{ch1}, H.channels(selCh_H(ch2)).name) && strcmp(EC_info.subj_names{1,subj}, params.storage.subjList{subj})
                        i_ch = cat(2, i_ch, ch2);
                    end
                end
            end 
            y = cat(2, y, trialsData.yVals(:,i_ch,:));
        end
        assert(nCh == size(y,2));
        
        
        % axes
%         nPlot = fb + (nn-1)*nCols + nCols;
        ax = subtightplot(nRows, nCols, nPlot, gap, marg_h, marg_w);
        hold on;            
        set(gca, 'FontSize',fontSize);  

        % ------------------------ plot: RBP ----------------------------------
        % select the same channels as for DTF 
%         i_ch = [];
%         nCh = 0;
%         for subj = 1:size(params.storage.subjList)
%             selCh_names = EC_all_info.info.selCh_groups.ch_names{nn,subj+1};    % +1 as there is a bug in the catting (msSEI_FC_MVARv2_fig_S.m)
%             nCh = nCh + size(selCh_names,2);
%             for ch1 = 1:size(selCh_names,2)
%                 for ch2 = 1:size(GRUP.chSubj_name{aa,fb,1},1)
%                     if strcmp(GRUP.chSubj_name{aa,fb,1}{ch2,1}, selCh_names{ch1}) && strcmp(GRUP.chSubj_name{aa,fb,1}{ch2,2}, params.storage.subjList{subj})
%                         i_ch = cat(2, i_ch, ch2);
%                     end
%                 end
%             end
%         end
%         assert(nCh == size(i_ch,2));
        
        % data to plot
        i_t = closestval(trialsData.xVals,time2plot(1)):closestval(trialsData.xVals,time2plot(2));
        x = trialsData.xVals{aa,fb}(i_t);
        y_avg = mean(y(i_t,:,:),2);
        y_sem =  sem(y(i_t,:,:),2);

        % >>> plot grouped channels activity in GRUP.trials{aa,freq} <<<
        for clz = 1:size(y,3)
%             plotband(x, y_avg(:,:,clz), y_sem(:,:,clz), clrs(clz,:));
            msSEI_plotband(x, y_avg(:,:,clz), y_sem(:,:,clz), clz);
        end

        % plot the AVG (non-specific) response (mean over clz)
        y_nonSpec = mean(y_avg,3);  % mean over clz
        plot(x, y_nonSpec, ':k', 'LineWidth',1);

        % plot t=0 line
        axis tight;
        xlim(time2plot);
        plot([0 0], ylim, ':k', 'LineWidth',1);

        % plot significance of grouped activations
%         if isfield(GRUP, 'H_trials')
%             plot2axes_signif_filled(GRUP.time{aa,fb}(i_t), GRUP.H_trials{aa,fb}(i_t), GRUP.P_trials{aa,fb}(i_t));
%         end

        % labels
        if nPlot > nCols*(nRows-1)
            xlabel('time (s)', 'FontSize',fontSize);
        else
            set(gca, 'XTickLabel',{''});
        end
        if mod(nPlot,nCols) == 1
%             ylabel({label_AA{aa,2},['N_C=' num2str(N_chnl(3),3)], ['N_P=' num2str(N_subj(3),2)]}, 'FontSize',fontSize, 'FontWeight','bold');   % , 'Rotation',0, 'HorizontalAlignment','right', 'VerticalAlignment','center'
            ylabel('RBP (z-score)', 'FontSize',fontSize);
        end
%         title([label_AA{aa,2} ': ' label_FB{freq,3}], 'FontSize',fontSize);
        title(['RBP: ' label_AA{aa,2}], 'FontSize',fontSize);
        
        % Tick Labels
        set(gca, 'XTick',[-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2]);
        if nPlot > nCols*(nRows-1)
            set(gca, 'XTickLabel',{'-2', '', '-1', '', '0', '', '1', '', '2'});        
        else
            set(gca, 'XTickLabel',{''});        
        end
        set(gca, 'YTick',[-1:0.025:1]);
        set(gca, 'YTickLabel',{''});

        box on;
        grid on;

        nPlot = nPlot+1;
    end
    
    %% ========================= row 2: MVAR-S - pre FB ==========================
    for fb = 1:size(label_FB,1)
        selFreq = label_FB{fb,2}; 
        disp([' - plotting: ' label_FB{fb,1}]);

        % axes
%         nPlot = fb + (nn-1)*nCols + nCols;
        ax = subtightplot(nRows, nCols, nPlot, gap, marg_h, marg_w);
        hold on;            
        set(gca, 'FontSize',fontSize);  
        
        % ------------------------ plot: RBP ----------------------------------
        % data to plot
%         i_t = closestval(EC_info.time,time2plot(1)):closestval(EC_info.time,time2plot(2));
%         x = EC_info.time(i_t);
        i_t = closestval(EC_time,time2plot(1)):closestval(EC_time,time2plot(2));
        x = EC_time(i_t);        
        i_f = closestval(EC_info.freq,selFreq(1)):closestval(EC_info.freq,selFreq(2));
        S_subj = squeeze(mean(EC_all.S(nn,nn,i_f,i_t,:,:),3));  % 3D: t x clz x subj, mean over freq  
        S_avg = nanmean(S_subj,3);  % 2D: t x clz
        S_sem = sem(S_subj,3);      % 2D: t x clz

        % >>> plot MVAR RBP <<<
        for clz = 1:size(S_avg,2)
            msSEI_plotband(x, S_avg(:,clz), S_sem(:,clz), clz);
        end

        % plot the AVG (non-specific) response (mean over clz)
        y_nonSpec = mean(S_avg,2);  % mean over clz
        plot(x, y_nonSpec, ':k', 'LineWidth',1);

        % plot t=0 line
        axis tight;
        xlim(time2plot);
        plot([0 0], ylim, ':k', 'LineWidth',1);

        % plot significance of grouped activations
%         if isfield(GRUP, 'H_trials')
%             plot2axes_signif_filled(GRUP.time{aa,fb}(i_t), GRUP.H_trials{aa,fb}(i_t), GRUP.P_trials{aa,fb}(i_t));
%         end

        % significance
        FC = S_subj;   % 3D: t x clz x subj
        vals = reshape(FC, [size(FC,1), 1, size(FC,2)*size(FC,3)]); % 3D: [t, 1, subj-clz], where 1 = "ch", subj-clz = "trials"
        labels = repmat([1,2], [1, size(FC,3)]);
        [P_vals, H_vals] = getSignificance(vals, labels, 0.05);
        plot2axes_signif_filled(x, H_vals, P_vals);

        % labels
        if nPlot > nCols*(nRows-1)
            xlabel('time (s)', 'FontSize',fontSize);
        else
            set(gca, 'XTickLabel',{''});
        end
        if mod(nPlot,nCols) == 1
%             ylabel({label_AA{aa,2},['N_C=' num2str(N_chnl(3),3)], ['N_P=' num2str(N_subj(3),2)]}, 'FontSize',fontSize, 'FontWeight','bold');   % , 'Rotation',0, 'HorizontalAlignment','right', 'VerticalAlignment','center'
            ylabel('RBP (z-score)', 'FontSize',fontSize);
        end
%         title([label_AA{aa,2} ': ' label_FB{freq,3}], 'FontSize',fontSize);
        title(['RBP: ' label_AA{aa,2}], 'FontSize',fontSize);
        
        % Tick Labels
        set(gca, 'XTick',[-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2]);
        if nPlot > nCols*(nRows-1)
            set(gca, 'XTickLabel',{'-2', '', '-1', '', '0', '', '1', '', '2'});        
        else
            set(gca, 'XTickLabel',{''});        
        end
        set(gca, 'YTick',[-1:0.025:1]);
        set(gca, 'YTickLabel',{''});

        box on;
        grid on;

        nPlot = nPlot+1;
    end % of fb
    
end % of nn

%% save
figName = [FC_outDirName '_RBP' '_' params.connectivity.normalize '_mvarS_v2'];
outDir = [params.storage.dir_results filesep params.storage.outName filesep 'msSEI_fig_suppl_mvarS_PSD'];
fig_save(f, figName, outDir, 'format','png', 'res',600);
% close(f);

return;


%% LOAD & PLOT EC.S -> EC_all_info & EC_all.S = 6D: nn x nn x freq x time x clz x subj
% !!! make sure the selection of FC_outDirName (etc) is uncommented !!!
msSEI_FC_MVARv2_fig_S;
assert(size(EC_all.S,6) == size(subjList,1));

%% ch2roi for Yeo7 -> GRUP ===> all chnls MNI
if load_ch2roi
    params.storage.subjList = subjList;
    clear GRUP;
    groupInfo = struct(...
        'brainAtlas', {{brainAtlas}}, ... % brain atlases: 'Yeo7','Yeo7_addLat','Yeo17','Yeo17_addLat','Mars','Mars_addLat'
        'freqBands', {label_FB(:,1)'}, ...  % freq. bands (FB)
        'signifOfActiv', {{'any_sgnf'}}, ...   % significance of activation
        'significanceLevel', P_level, ...      % significance level for ROI (e.g. DMN)
        'significanceTimeThr', t_siglen, ...    % in [s], significance time threshold (ch is signif. if its significance > threshold) 
        'GW_selection', {{'any'}}, ...      % grey / white matter: choices: 'gm','wm','csf','bone', 'any'
        'GW_threshold', 60, ...    % probability, in [%]
        'plotStyle', {{'allFB_allROI_subjMean','allFB_allROI'}}, ... 'oneFB_allROI','allFB_allROI','spectra_SNR'
        'maxDistToROI', 10, ...   % in [mm], maximum distance (tolerance) to region of interest
        'minSubjInROI', 1, ...  % minimum of subjects (patients) in a ROI (if a ROI has less patients, it is not considered in plots)
        'subtractAvgResponse', false, ...    % subtracts mean response over classes (conditions) from each clz
        'dirSuffix', dirSuffix, ... % adds suffix to 'ch2roi' directory name -> ch2roi_outDirSuffix
        'figSuffix', '', ... % adds suffix to figure name
        'time2plot', t_sel ... % in [s], w.r.t. cutting point (= 0)
        );
    ch2roi_wrapper;
end
assert(strcmp(label_FB{fb,1}, groupInfo.list_FB{freq}));


    % chnls MNI
%     tmp_ch = EC_info.ch_MNI{nn};               % sel chnls MNI
%     tmp_ch = cat(2, tmp_ch, nn*ones(size(tmp_ch,1),1)); % add value = nn (was significant)
%     chnls_MNI_VAL = cat(1, chnls_MNI_VAL, tmp_ch);      % add (again) the significant channels, (note that they are duplicate)
% %     params.plot_brain3D.chVals_colorMap = cat(1, params.plot_brain3D.chVals_colorMap, msSwitch_get_clr4aa(sel_AA{nn}));
%     params.plot_brain3D.chVals_colorMap = cat(1, params.plot_brain3D.chVals_colorMap, chnls_clrs_nn(nn,:));

%% ========================= FIGURE: RBP + MVAR-S ==========================
nCols = size(label_FB,1);   % ~ freq
nRows = 4;  % PSD, MVAR-S: DMN & DAN
nPlot = 1;

fig_W = 24;         % in [cm]
fig_H = 4*nRows;    % in [cm], 3 cm per row
f = fig_make('fig_W',fig_W, 'fig_H',fig_H);

marg_h = [0.10 0.10];   % margin from [bottom, top]
marg_w = [0.08 0.02];   % margin from [L, R] side
gap = [0.02, 0.02];     % between axes from [top, side] 

fontSize = 10;
time2plot = [-2.3, 2.3];

for nn = 1:size(sel_AA,1)
    [aa_found, aa] = ismember(sel_AA{nn}, label_AA(:,1));   % aa = index to GRUP !!!
    assert(aa_found);
    assert(strcmp(label_AA{aa,1}, sel_AA{nn}));        
    
    %% =========================== row 1: PLOT PSD - per FB ==========================
    for fb = 1:size(label_FB,1)
        selFreq = label_FB{fb,2}; 
        disp([' - plotting: ' label_FB{fb,1}]);

        % axes
%         nPlot = fb + (nn-1)*nCols + nCols;
        ax = subtightplot(nRows, nCols, nPlot, gap, marg_h, marg_w);
        hold on;            
        set(gca, 'FontSize',fontSize);  
        
        % ------------------------ plot: RBP ----------------------------------
        % select the same channels as for DTF 
        i_ch = [];
        nCh = 0;
        for subj = 1:size(params.storage.subjList)
            selCh_names = EC_all_info.info.selCh_groups.ch_names{nn,subj+1};    % +1 as there is a bug in the catting (msSEI_FC_MVARv2_fig_S.m)
            nCh = nCh + size(selCh_names,2);
            for ch1 = 1:size(selCh_names,2)
                for ch2 = 1:size(GRUP.chSubj_name{aa,fb,1},1)
                    if strcmp(GRUP.chSubj_name{aa,fb,1}{ch2,1}, selCh_names{ch1}) && strcmp(GRUP.chSubj_name{aa,fb,1}{ch2,2}, params.storage.subjList{subj})
                        i_ch = cat(2, i_ch, ch2);
                    end
                end
            end
        end
        assert(nCh == size(i_ch,2));
        
        % data to plot
        i_t = closestval(GRUP.time{aa,fb},time2plot(1)):closestval(GRUP.time{aa,fb},time2plot(2));
        x = GRUP.time{aa,fb}(i_t);
        y_avg = mean(GRUP.trials{aa,fb}(i_t,i_ch,:),2);
        y_sem =  sem(GRUP.trials{aa,fb}(i_t,i_ch,:),2);

        % >>> plot grouped channels activity in GRUP.trials{aa,freq} <<<
        clrs = GRUP.info.clzColor;
        for clz = 1:size(GRUP.trials{aa,fb},3)
%             plotband(x, y_avg(:,:,clz), y_sem(:,:,clz), clrs(clz,:));
            msSEI_plotband(x, y_avg(:,:,clz), y_sem(:,:,clz), clz);
        end

        % plot the AVG (non-specific) response (mean over clz)
        y_nonSpec = mean(y_avg,3);  % mean over clz
        plot(x, y_nonSpec, ':k', 'LineWidth',1);

        % plot t=0 line
        axis tight;
        xlim(time2plot);
        plot([0 0], ylim, ':k', 'LineWidth',1);

        % plot significance of grouped activations
        if isfield(GRUP, 'H_trials')
            plot2axes_signif_filled(GRUP.time{aa,fb}(i_t), GRUP.H_trials{aa,fb}(i_t), GRUP.P_trials{aa,fb}(i_t));
        end

        % labels
        if nPlot > nCols*(nRows-1)
            xlabel('time (s)', 'FontSize',fontSize);
        else
            set(gca, 'XTickLabel',{''});
        end
        if mod(nPlot,nCols) == 1
%             ylabel({label_AA{aa,2},['N_C=' num2str(N_chnl(3),3)], ['N_P=' num2str(N_subj(3),2)]}, 'FontSize',fontSize, 'FontWeight','bold');   % , 'Rotation',0, 'HorizontalAlignment','right', 'VerticalAlignment','center'
            ylabel('RBP (z-score)', 'FontSize',fontSize);
        end
%         title([label_AA{aa,2} ': ' label_FB{freq,3}], 'FontSize',fontSize);
        title(['RBP: ' label_AA{aa,2}], 'FontSize',fontSize);
        
        % Tick Labels
        set(gca, 'XTick',[-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2]);
        if nPlot > nCols*(nRows-1)
            set(gca, 'XTickLabel',{'-2', '', '-1', '', '0', '', '1', '', '2'});        
        else
            set(gca, 'XTickLabel',{''});        
        end
        set(gca, 'YTick',[-1:0.025:1]);
        set(gca, 'YTickLabel',{''});

        box on;
        grid on;

        nPlot = nPlot+1;
    end
    
    %% ========================= row 2: MVAR-S - pre FB ==========================
    for fb = 1:size(label_FB,1)
        selFreq = label_FB{fb,2}; 
        disp([' - plotting: ' label_FB{fb,1}]);

        % axes
%         nPlot = fb + (nn-1)*nCols + nCols;
        ax = subtightplot(nRows, nCols, nPlot, gap, marg_h, marg_w);
        hold on;            
        set(gca, 'FontSize',fontSize);  
        
        % ------------------------ plot: RBP ----------------------------------
        % data to plot
%         i_t = closestval(EC_info.time,time2plot(1)):closestval(EC_info.time,time2plot(2));
%         x = EC_info.time(i_t);
        i_t = closestval(EC_time,time2plot(1)):closestval(EC_time,time2plot(2));
        x = EC_time(i_t);        
        i_f = closestval(EC_info.freq,selFreq(1)):closestval(EC_info.freq,selFreq(2));
        S_subj = squeeze(mean(EC_all.S(nn,nn,i_f,i_t,:,:),3));  % 3D: t x clz x subj, mean over freq  
        S_avg = nanmean(S_subj,3);  % 2D: t x clz
        S_sem = sem(S_subj,3);      % 2D: t x clz

        % >>> plot MVAR RBP <<<
        for clz = 1:size(S_avg,2)
            msSEI_plotband(x, S_avg(:,clz), S_sem(:,clz), clz);
        end

        % plot the AVG (non-specific) response (mean over clz)
        y_nonSpec = mean(S_avg,2);  % mean over clz
        plot(x, y_nonSpec, ':k', 'LineWidth',1);

        % plot t=0 line
        axis tight;
        xlim(time2plot);
        plot([0 0], ylim, ':k', 'LineWidth',1);

        % plot significance of grouped activations
%         if isfield(GRUP, 'H_trials')
%             plot2axes_signif_filled(GRUP.time{aa,fb}(i_t), GRUP.H_trials{aa,fb}(i_t), GRUP.P_trials{aa,fb}(i_t));
%         end

        % significance
        FC = S_subj;   % 3D: t x clz x subj
        vals = reshape(FC, [size(FC,1), 1, size(FC,2)*size(FC,3)]); % 3D: [t, 1, subj-clz], where 1 = "ch", subj-clz = "trials"
        labels = repmat([1,2], [1, size(FC,3)]);
        [P_vals, H_vals] = getSignificance(vals, labels, 0.05);
        plot2axes_signif_filled(x, H_vals, P_vals);

        % labels
        if nPlot > nCols*(nRows-1)
            xlabel('time (s)', 'FontSize',fontSize);
        else
            set(gca, 'XTickLabel',{''});
        end
        if mod(nPlot,nCols) == 1
%             ylabel({label_AA{aa,2},['N_C=' num2str(N_chnl(3),3)], ['N_P=' num2str(N_subj(3),2)]}, 'FontSize',fontSize, 'FontWeight','bold');   % , 'Rotation',0, 'HorizontalAlignment','right', 'VerticalAlignment','center'
            ylabel('RBP (z-score)', 'FontSize',fontSize);
        end
%         title([label_AA{aa,2} ': ' label_FB{freq,3}], 'FontSize',fontSize);
        title(['RBP: ' label_AA{aa,2}], 'FontSize',fontSize);
        
        % Tick Labels
        set(gca, 'XTick',[-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2]);
        if nPlot > nCols*(nRows-1)
            set(gca, 'XTickLabel',{'-2', '', '-1', '', '0', '', '1', '', '2'});        
        else
            set(gca, 'XTickLabel',{''});        
        end
        set(gca, 'YTick',[-1:0.025:1]);
        set(gca, 'YTickLabel',{''});

        box on;
        grid on;

        nPlot = nPlot+1;
    end % of fb
    
end % of nn

%% save
figName = [FC_outDirName '_RBP' '_' params.connectivity.normalize '_mvarS_v1'];
outDir = [params.storage.dir_results filesep params.storage.outName filesep 'msSEI_fig_suppl_mvarS_PSD'];
fig_save(f, figName, outDir, 'format','png', 'res',600);
% close(f); 
