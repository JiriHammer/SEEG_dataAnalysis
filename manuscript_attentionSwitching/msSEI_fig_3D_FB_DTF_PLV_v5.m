%% ms SEI fig: funct6ional & effective connectivity - 3D brains + PLV + DTF
% "final" figure, v5  :o0
% works for two NN: DMN & DAN
% FIGURE 1
% - 3D brains - sel. chnls
% FIGURE 2
% cols = FBs (delta, ..., hiGamma)
% rows =
%   - PLV
%   - DTF: DMN -> DAN
%   - DTF: DAN -> DMN
% jobs: 
%   - DTF = msSEI_FC_MVARv2_job (filtered in 3 FB: LFB, LGB, HGB)
%   - PLV = msSEI_FC_PLV_job

% (c) Jiri, Jan24
% based on: msSEI_fig_3D_FB_DTF_v4.m

%% settings: relaxed conditions (p-level, length of significance, ...)
dirSuffix = '_msFig_3D_FB_DTF_PLV_v5';
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
MVAR_order = 17;

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

%% 3D brain params & view angles -> sel_viewAngles
sel_viewAngles =  {
    'sagittal L', [-90, 0];          % sagittal R = 100,10 (90, 0)
    'axial',      [0, 90];          % axial = 0,90    
    'sagittal R', [90, 0];         % sagittal L = -100,10 (-90, 0)
    };
    
%% ========================= FIGURE 1: 3D brains ==========================
if plot_fig1_brain3D
%% load DTF results -> M_all, EC_info ===> used chnls MNI
params.connectivity.normalize = 'zscore_time';
if load_dtfResults
    fc_measure = 'DTF';
    [M_all, EC_info] = msSEI_getDTF_results(params, FC_outDirName, fc_measure);
    EC_time = EC_info.time;   % time axis for DTF
    params.storage.subjList = EC_info.subjList; % updates subj (those having more than 3 channels/ROI)
end
disp(['DTF loaded. nSubj = ' num2str(size(EC_info.subjList,1))]);
    
%% TBD???  ch2roi for Yeo7 -> GRUP ===> all chnls MNI
% if load_ch2roi
%     clear GRUP;
%     groupInfo = struct(...
%         'brainAtlas', {{brainAtlas}}, ... % brain atlases: 'Yeo7','Yeo7_addLat','Yeo17','Yeo17_addLat','Mars','Mars_addLat'
%         'freqBands', {{label_FB{fb,1}}}, ...  % freq. bands (FB)
%         'signifOfActiv', {{'any_sgnf'}}, ...   % significance of activation
%         'significanceLevel', P_level, ...      % significance level for ROI (e.g. DMN)
%         'significanceTimeThr', t_siglen, ...    % in [s], significance time threshold (ch is signif. if its significance > threshold) 
%         'GW_selection', {{'any'}}, ...      % grey / white matter: choices: 'gm','wm','csf','bone', 'any'
%         'GW_threshold', 60, ...    % probability, in [%]
%         'plotStyle', {{'allFB_allROI_subjMean','allFB_allROI'}}, ... 'oneFB_allROI','allFB_allROI','spectra_SNR'
%         'maxDistToROI', 10, ...   % in [mm], maximum distance (tolerance) to region of interest
%         'minSubjInROI', 1, ...  % minimum of subjects (patients) in a ROI (if a ROI has less patients, it is not considered in plots)
%         'subtractAvgResponse', false, ...    % subtracts mean response over classes (conditions) from each clz
%         'dirSuffix', dirSuffix, ... % adds suffix to 'ch2roi' directory name -> ch2roi_outDirSuffix
%         'figSuffix', '', ... % adds suffix to figure name
%         'time2plot', t_sel ... % in [s], w.r.t. cutting point (= 0)
%         );
%     ch2roi_wrapper;
% end
% assert(strcmp(label_FB{fb,1}, groupInfo.list_FB{freq}));

%% figure
nCols = size(sel_viewAngles,1);   % ~ freq
nRows = 1;
fig_W = 24;         % in [cm]
fig_H = 10;    % in [cm], 3 cm per row
f = fig_make('fig_W',fig_W, 'fig_H',fig_H);
marg_h = [0.05 0.05];   % margin from [bottom, top]
marg_w = [0.05 0.10];   % margin from [L, R] side
gap = [-0.01, -0.01];     % between axes from [top, side]  

%  ------------------------ PLOT: 3D brain for both NN ---------------
%     marg_w = [0.02 0.10];   % margin from [L, R] side
% settings
params.plot_brain3D = struct(...        % user interface structure: holds most (but not all!) of the user settings, see also inside brain3D_plot.m
    'plot_brain3D', true, ...           % 3D brain model
    'volumes2plot', {{...               % volumes to plot (coregistered .nii files)
        'colin27_t1_tal_lin masked', [1 1 1], 0.20;  ...     % 1x3 cell: name, color, transparency (other choices for norm brain = wc1T1_colin27)
        }}, ...          
    'chVals_asCircles', true, ...       % plots channel values as color circles
    'chVals_asPatches', false, ...       % each channel value is mapped onto surface (but may be overwritten by latter points!), see also distanceThreshold for focality inside brain3D_plot.m
    'chVals_colorMap', [0 0 0], ...    % colormap of channel values
    'chVals_patchTransp', 0.2, ...      % transparency of channel values
    'circleSizeLims', [5 5], ...       % in points, maps chnl values to circle size
    'size_interpolate', 1.0, ...        % in [mm], voxel size to which the brain is interpolated
    'backgroundColor', 'w', ...         % figure / axes background color, good choices: 'w', 'k'
    'text_LR', false, ...                % marks L / R side
    'text_tag', '', ...    % string
    'viewAngle', [0, 90], ...           % axial = 0,90; sagittal = 100,10; coronal front = 0,180; coronal back = 0,0
    'visible_axis', false ...           % if axis & their labels are visible
); 

% 3D brain: include Yeo7 NNs
chnls_MNI_VAL = [];
chnls_clrs_nn = [
    [139, 0, 0]./255;  % DMN = dark red
    [0, 80, 0]./255;  % DAN = dark green
];
for nn = 1:size(sel_AA,1) 
    [aa_found, aa] = ismember(sel_AA{nn}, label_AA(:,1));   % aa = index to GRUP !!!
    assert(aa_found);
    assert(strcmp(label_AA{aa,1}, sel_AA{nn}));

    % include Yeo7 NNs
    params.plot_brain3D.volumes2plot{end+1,1} = sel_AA{nn};
    params.plot_brain3D.volumes2plot{end,2} = msSwitch_get_clr4aa(sel_AA{nn});
    params.plot_brain3D.volumes2plot{end,3} = 0.2;

    % chnls MNI
    tmp_ch = EC_info.ch_MNI{nn};               % sel chnls MNI
    tmp_ch = cat(2, tmp_ch, nn*ones(size(tmp_ch,1),1)); % add value = nn (was significant)
    chnls_MNI_VAL = cat(1, chnls_MNI_VAL, tmp_ch);      % add (again) the significant channels, (note that they are duplicate)
%     params.plot_brain3D.chVals_colorMap = cat(1, params.plot_brain3D.chVals_colorMap, msSwitch_get_clr4aa(sel_AA{nn}));
    params.plot_brain3D.chVals_colorMap = cat(1, params.plot_brain3D.chVals_colorMap, chnls_clrs_nn(nn,:));

end % of nn

% brain volumes -> surfaces
[brainVols, clrmap, inds_vol] = brainVolumes_fixed(params, params.plot_brain3D);

% plot: 3D brains: different view angles
for v = 1:size(sel_viewAngles,1)   

    % project channels MNI
    chnls_MNI_VAL_toPlot = chnls_MNI_VAL;
    if strcmp(sel_viewAngles{v,1}, 'sagittal L')
        chnls_MNI_VAL_toPlot(chnls_MNI_VAL(:,1)>=0,:) = [];  % only channels on L side
        chnls_MNI_VAL_toPlot(:,1) = min(brainVols{1}.xi);   % all channels on L side - up front!
    elseif strcmp(sel_viewAngles{v,1}, 'axial')
        chnls_MNI_VAL_toPlot(:,3) = max(brainVols{1}.zi);  % all channels on top - up front!
    elseif strcmp(sel_viewAngles{v,1}, 'sagittal R')
        chnls_MNI_VAL_toPlot(chnls_MNI_VAL(:,1)<0,:) = [];  % only channels on R side
        chnls_MNI_VAL_toPlot(:,1) = max(brainVols{1}.xi);   % all channels on R side - up front!
    end     

    % axes
%     nPlot = fb + (v-1)*nCols;
    nPlot = v;
    ax = subtightplot(nRows, nCols, nPlot, gap, marg_h, marg_w);
    hold on;

    % set: plot_brain3D & plotInfo
    plotInfo = params.plot_brain3D;
    plotInfo.axHandle = ax;
    plotInfo.figHandle = f;
    plotInfo.chVals_lims = [0, 2];
    
    % PLOT 3D brain vols
    for vol = 1:size(brainVols,2)
        if brainVols{vol}.loaded
            plotInfo.thisVolume = vol;
            % plot all chnl values at once: uncomment line 137 in brain3D_plot.m
            h = brain3D_plot(brainVols{vol}, chnls_MNI_VAL_toPlot, plotInfo);   
        end
    end    

    % title: FB
%     if v == 1
%         title(label_FB{fb,3}, 'FontSize',fontSize);  % , 'FontWeight','normal'
%     end

    % set view angle
    view(sel_viewAngles{v,2});
end % view angles

%% save
figName = [FC_outDirName '_brain3D_v3'];
outDir = [params.storage.dir_results filesep params.storage.outName filesep 'msSEI_fig_3D_FB_DTF_PLV_v5'];
fig_save(f, figName, outDir, 'format','png', 'res',600);
% close(f); 
end % of plot_fig1_brain3D

%% ========================= FIGURE 2: PLV + DTF ==========================
nCols = size(label_FB,1);   % ~ freq
nRows = 3;  % PLV, DTF: DMN->DAN, DTF: DAN->DMN
nPlot = 1;

fig_W = 18;         % in [cm]
fig_H = 4*nRows;    % in [cm], 3 cm per row
f = fig_make('fig_W',fig_W, 'fig_H',fig_H);

marg_h = [0.13 0.10];   % margin from [bottom, top]
marg_w = [0.08 0.02];   % margin from [L, R] side
gap = [0.02, 0.05];     % between axes from [top, side] 

fontSize = 10;
time2plot = [-2.3, 2.3];

% ========================= row 1: PLV ==========================

% load PLV results -> M_all, EC_info
fc_measure = 'PLV';
PLV_outDirName = 'PLV_v3';
params.connectivity.normalize = norm_fcMatrix;
if load_plvResults
    disp(['Loading: ' fc_measure ' ...']);
    [PLV_all, PLV_info] = msSEI_getDTF_results(params, PLV_outDirName, fc_measure);
    PLV_time = PLV_info.time;   % time axis for DTF
%     params.storage.subjList = PLV_info.subjList; % updates subj (those having more than 3 channels/ROI)
end
EC_info = PLV_info;
EC_time = PLV_time;
disp([fc_measure ' loaded. nSubj = ' num2str(size(EC_info.subjList,1))]);

% =========================== PLOT PLV - per FB ==========================
yLims = [-0.08, 0.18];
for fb = 1:size(label_FB,1)
    disp([' - plotting: ' label_FB{fb,1}]);
    
    % ------------------------ plot: PLV ----------------------------------
    nn = 1; 
    src = 1;    % DMN
    trg = 2;    % DAN    
    [aa_found, aa] = ismember(sel_AA{nn}, label_AA(:,1));   % aa = index to GRUP !!!
    assert(aa_found);
    assert(strcmp(label_AA{aa,1}, sel_AA{nn}));        
    selFreq = label_FB{fb,2};
    i_f = closestval(EC_info.freq, selFreq(1)):closestval(EC_info.freq, selFreq(2))
    if fb == 1
        i_f = 1
%     elseif fb == 2
%         i_f = 2
%     elseif fb == 3
%         i_f = 3    
    end
    i_t = closestval(EC_time, time2plot(1)):closestval(EC_time, time2plot(2));  
    EC_info.time = EC_time(i_t);
    EC_info.P_level = P_level;    

    % PLV data -> FC
    FC = squeeze(nanmean(PLV_all(trg,src,i_f,i_t,:,:),3));   % 3D: t x clz x subj, mean over freq

    % axes
%     nPlot = fb + (nn-1)*nCols + 4*nCols;
    ax = subtightplot(nRows, nCols, nPlot, gap, marg_h, marg_w);
    hold on;            
    set(gca, 'FontSize',fontSize);  
%         ax.YAxisLocation = 'right';
%     yyaxis right;

    % title
    [tf, i_src] = ismember(sel_AA{src}, label_AA(:,1)); assert(tf);
    [tf, i_trg] = ismember(sel_AA{trg}, label_AA(:,1)); assert(tf);    
    EC_info.str_title = ['PLV: ' label_AA{i_src,2} ' - ' label_AA{i_trg,2}];

    % >>> plot <<<
%     plotEC_plotband_2NN(FC, EC_info);
    for clz = 1:size(FC,2)
%             plotband(EC_info.time, nanmean(FC(:,clz,:),3), sem(FC(:,clz,:),3), EC_info.info.colors(clz,:), clz);
        msSEI_plotband(EC_info.time, nanmean(FC(:,clz,:),3), sem(FC(:,clz,:),3), clz);
    end
    axis tight;
    if fb >=2
        ylim(yLims);
    end
    
    plot([0 0], ylim, ':k', 'LineWidth',1);

    % significance
    vals = reshape(FC, [size(FC,1), 1, size(FC,2)*size(FC,3)]); % 3D: [t, 1, subj-clz], where 1 = "ch", subj-clz = "trials"
    labels = repmat([1,2], [1, size(FC,3)]);
    [P_vals, H_vals] = getSignificance(vals, labels, 0.001);
%     [P_vals, H_vals] = getSignificance(vals, labels, EC_info.P_level);
    plot2axes_signif_filled(EC_info.time, H_vals, P_vals);

    % labels
    if nPlot > nCols*(nRows-1)
        xlabel('Time (s)', 'FontSize',fontSize);
    else
        set(gca, 'XTickLabel',{''});
    end
    if mod(nPlot,nCols) == 1
%         ylabel('PLV (z-score)');
    else
%         set(gca, 'YTickLabel',{''});
    end
%     title(EC_info.str_title);
    title(label_FB{fb,3}, 'FontSize',fontSize);

    set(gca, 'XTick',[-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2]);
    if nPlot > nCols*(nRows-1)
        set(gca, 'XTickLabel',{'-2', '', '-1', '', '0', '', '1', '', '2'});        
    else
        set(gca, 'XTickLabel',{''});        
    end

    grid on;
    box on;

    nPlot = nPlot+1;
end
disp('Plotting PLV done.');

% ========================= row 2: DTF ==========================
% load DTF results -> M_all, EC_info ===> used chnls MNI
if load_dtfResults
    fc_measure = 'DTF';
    disp(['Loading: ' fc_measure ' ...']);
    [DTF_all, DTF_info] = msSEI_getDTF_results(params, FC_outDirName, fc_measure);
    DTF_time = DTF_info.time;   % time axis for DTF
%     params.storage.subjList = EC_info.subjList; % updates subj (those having more than 3 channels/ROI)
end
EC_info = DTF_info;
EC_time = DTF_time;
disp(['DTF loaded. nSubj = ' num2str(size(EC_info.subjList,1))]);

% =========================== PLOT DTF - per FB ==========================
yLims = [-0.25, 0.25];
for fb = 1:size(label_FB,1)
    disp([' - plotting: ' label_FB{fb,1}]);
    
    % ------------------------ plot: DTF ----------------------------------
    for nn = 1:size(sel_AA,1) 
        [aa_found, aa] = ismember(sel_AA{nn}, label_AA(:,1));   % aa = index to GRUP !!!
        assert(aa_found);
        assert(strcmp(label_AA{aa,1}, sel_AA{nn}));        
        selFreq = label_FB{fb,2};
        i_f = closestval(EC_info.freq, selFreq(1)):closestval(EC_info.freq, selFreq(2))
        i_t = closestval(EC_time, time2plot(1)):closestval(EC_time, time2plot(2));  
        EC_info.time = EC_time(i_t);
        EC_info.P_level = P_level;    

        % define source (src) & target (trg)
        if nn == 1
            src = 1;    % DMN
            trg = 2;    % DAN
        else
            src = 2;    % DAN
            trg = 1;    % DMN
        end

        % DTF data -> FC
        FC = squeeze(nanmean(DTF_all(trg,src,i_f,i_t,:,:),3));   % 3D: t x clz x subj, mean over freq

        % axes
        nPlot = fb + (nn-1)*nCols + nCols;
        ax = subtightplot(nRows, nCols, nPlot, gap, marg_h, marg_w);
        hold on;            
        set(gca, 'FontSize',fontSize);  
%         ax.YAxisLocation = 'right';
    %     yyaxis right;

        % title
        [tf, i_src] = ismember(sel_AA{src}, label_AA(:,1)); assert(tf);
        [tf, i_trg] = ismember(sel_AA{trg}, label_AA(:,1)); assert(tf);    
    %     EC_info.str_title = ['DTF: ' sel_AA{src} ' -> ' sel_AA{trg}];
        EC_info.str_title = ['DTF: ' label_AA{i_src,2} ' -> ' label_AA{i_trg,2}];

        % >>> plot DTF <<<
    %     plotEC_plotband_2NN(FC, EC_info);
        for clz = 1:size(FC,2)
%             plotband(EC_info.time, nanmean(FC(:,clz,:),3), sem(FC(:,clz,:),3), EC_info.info.colors(clz,:), clz);
            msSEI_plotband(EC_info.time, nanmean(FC(:,clz,:),3), sem(FC(:,clz,:),3), clz);
        end
        axis tight;
        ylim(yLims);
        
        plot([0 0], ylim, ':k', 'LineWidth',1);

        % significance
        vals = reshape(FC, [size(FC,1), 1, size(FC,2)*size(FC,3)]); % 3D: [t, 1, subj-clz], where 1 = "ch", subj-clz = "trials"
        labels = repmat([1,2], [1, size(FC,3)]);
        [P_vals, H_vals] = getSignificance(vals, labels, 0.001);
%         [P_vals, H_vals] = getSignificance(vals, labels, EC_info.P_level);
        plot2axes_signif_filled(EC_info.time, H_vals, P_vals);
        
        % labels
        if nPlot > nCols*(nRows-1)
            xlabel('Time (s)', 'FontSize',fontSize);
        else
            set(gca, 'XTickLabel',{''});
        end
        if mod(nPlot,nCols) == 1
%             ylabel('DTF (z-score)');
        else
%             set(gca, 'YTickLabel',{''});
        end
%         title(EC_info.str_title);

        set(gca, 'XTick',[-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2]);
        if nPlot > nCols*(nRows-1)
            set(gca, 'XTickLabel',{'-2', '', '-1', '', '0', '', '1', '', '2'});        
        else
            set(gca, 'XTickLabel',{''});        
        end
%         ylim([-0.25, 0.25]);

        grid on;
        box on;
        
        nPlot = nPlot+1;
    end  % of nn
end

%% save
figName = [FC_outDirName '_PLV_v4' '_' params.connectivity.normalize '_DTF_v5'];
outDir = [params.storage.dir_results filesep params.storage.outName filesep 'msSEI_fig_3D_FB_DTF_PLV_v5'];
% fig_save(f, figName, outDir, 'format','png', 'res',600);
fig_save(f, figName, outDir, 'res',600);
% close(f); 

%% ms txt: infor about number of channels
disp('number of channels:');
for nn = 1:size(sel_AA,1)
    nCh = sum(EC_info.nCh_group(:,nn),1);
    nSubj = size(EC_info.nCh_group,1);
    nAvg = mean(EC_info.nCh_group(:,nn),1);
    nSem = sem(EC_info.nCh_group(:,nn),1);
    disp([sel_AA{nn} ': nCh = ' num2str(nCh), ', nSubj = ' num2str(nSubj)]);
    disp([sel_AA{nn} ': nCh/subj, mean +/- SEM = ' num2str(nAvg), ' +/- ' num2str(nSem)]);
end
