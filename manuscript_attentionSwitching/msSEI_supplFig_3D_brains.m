%% ms SEI: suppl. fig - 3D brains (highlight DMN & DAN)
% 25 subj - v31

% (c) Jiri, Jan24

%% settings: v1 - relaxed conditions (p-level, length of significance, ...)
dirSuffix = '_msSEI_supplFig_3D_brains_v1';
t_sel = [-1.0, 2.0];    % the true sel time will be cut by 1 s on edges (= template length): [-2, 2] -> [-1, 1]
P_level = 0.05;         % significance level
t_siglen = 0.10;        % in [s], significance time threshold (ch is signif. if its significance > threshold) 
brainAtlas = 'Yeo7';
load_ch2roi = true;

%% get job-related settings -> params
% --- specify analysis to be performed -> runner
runner = struct;
runner.newJobs = false; % set to 'false' if the results are already computed
runner.newFDR = false;  % true = updates the FDR correction over all patients (and freq. bands, ...) 

% --- job settings: specific (define spectral method) -> params
spectralMethod = 'v31_stft_session_bip';           % ~ v31 (BIP), similar to v12 (STFT base = sess, no zscore, no rej, FB = loFreq0)
% spectralMethod = 'v30_stft_baseRS_bip';         % ~ v30 (BIP), STFT with RS baseline
% spectralMethod = 'v29_stft_baseRS_bip';         % ~ v29 (BIP), STFT with RS baseline
% spectralMethod = 'v24_stft_baseRS_car';         % ~ v24 (CAR), STFT with RS baseline
% spectralMethod = 'v24_stft_baseRS_bip';         % ~ v24 (BIP), STFT with RS baseline
% spectralMethod = 'stft_sessions';       % ~ v12
% spectralMethod = 'stft_baseRS';         % ~ v11  (CAR)
% spectralMethod = 'stft_baseRS_bip';     % ~ v11b (BIP)

params = msSEI_getJobParams(runner, spectralMethod, P_level);


%% subj list
    subjList = {
        '20_PR3', 'P1';
        '20_PR4', 'P2';
        '20_PR5', 'P3';
        '20_PR6', 'P4';
        '20_PR8', 'P5';
        '20_PR9', 'P6';
        '20_PR10', 'P7';
        '20_PR11', 'P8';
        '20_PR13', 'P9';
        '20_PR12', 'P10';
        '20_PR14', 'P11';
        '21_PR1', 'P12';
        '21_PR2', 'P13';
        '21_PR4', 'P14';
        '21_PR5', 'P15';
        '21_PR6', 'P16';
        '21_PR7', 'P17';
        '22_PR1', 'P18';
        '22_PR3', 'P19';
        '22_PR4', 'P20';
        '22_PR5', 'P21';
        '22_PR8', 'P22';
        '23_PR2', 'P23';
        '23_PR3', 'P24';
        '23_PR4', 'P25';
    };
% params.storage.subjList = subjList;

%% settings: AA -> label_AA
label_AA = msSEI_get_label_AA(lower(brainAtlas));

sel_AA = {
    'Default';
    'Dorsal Attention';
};

%% 3D brain params & view angles -> sel_viewAngles
sel_viewAngles =  {
%     'sagittal L', [-90, 0];          % sagittal R = 100,10 (90, 0)
    'axial',      [0, 90];          % axial = 0,90    
%     'sagittal R', [90, 0];         % sagittal L = -100,10 (-90, 0)
    };
    
%% ========================= FIGURE =======================================
nCols = 5; 
nRows = 5;  
assert(size(params.storage.subjList,1) <= nCols*nRows);

fig_W = 24;     % in [cm]
fig_H = 24;     % in [cm]
f = fig_make('fig_W',fig_W, 'fig_H',fig_H);
marg_h = [0.02 0.02];   % margin from [bottom, top]
marg_w = [0.02 0.02];   % margin from [L, R] side
gap = [-0.02, -0.02];     % between axes from [top, side]  
fontSize = 10;

% chnls_clrs_nn = [
%     [139, 0, 0]./255;  % DMN = dark red
%     [0, 80, 0]./255;  % DAN = dark green
% ];
% chnls_clrs_nn = [
%     [255, 0, 0]./255;  % DMN = red
%     [255, 0, 255]./255;  % DAN = magenta
% ];
chnls_clrs_nn = [
    [0, 0, 255]./255;  % DMN = blue
%     [3 213 255]./255;  % DMN = light blue
    [255, 0, 0]./255;  % DAN = red
];

% =============================== PLOT ===================================
nPlot = 1;
for subj = 1:size(subjList,1)
    params.storage.subjList = [];
    params.storage.subjList = subjList(subj,:);
    subjTag = subjList{subj};
    
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
        'circleSizeLims', [8 8], ...       % in points, maps chnl values to circle size
        'size_interpolate', 1.0, ...        % in [mm], voxel size to which the brain is interpolated
        'backgroundColor', 'w', ...         % figure / axes background color, good choices: 'w', 'k'
        'text_LR', false, ...                % marks L / R side
        'text_tag', '', ...    % string
        'viewAngle', [0, 90], ...           % axial = 0,90; sagittal = 100,10; coronal front = 0,180; coronal back = 0,0
        'visible_axis', false ...           % if axis & their labels are visible
    ); 
    
    %% ch2roi for singleSubj Yeo7 -> GRUP 
    clear GRUP;
    groupInfo = struct(...
        'brainAtlas', {{brainAtlas}}, ... % brain atlases: 'Yeo7','Yeo7_addLat','Yeo17','Yeo17_addLat','Mars','Mars_addLat'
        'freqBands', {{'hiGamma'}}, ...  % freq. bands (FB)
        'signifOfActiv', {{'any_sgnf'}}, ...   % significance of activation
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

    %% load cacheFile -> H, selCh_H_resp
    cacheFile = [params.storage.dir_results filesep params.storage.outName filesep subjTag filesep 'cacheFile.mat'];
    assert(exist(cacheFile,'file') == 2);
    clear H selCh_H_resp;
    load(cacheFile, 'H', 'selCh_H_resp');   
    selCh_H = selCh_H_resp; 
    
    %% all chnls_MNI_VAL -> [x,y,z,0];    
    chnls_MNI_VAL = [];    
    for ch = 1:size(selCh_H,2)
        tmp_ch = [H.channels(selCh_H(ch)).MNI_x, H.channels(selCh_H(ch)).MNI_y, H.channels(selCh_H(ch)).MNI_z, 0]; % 0 -> black dot
        chnls_MNI_VAL = cat(1, chnls_MNI_VAL, tmp_ch);
    end
    params.plot_brain3D.chVals_colorMap = [0 0 0];  % black dots = all
        
    %% sel_AA -> chnls_MNI_VAL + nn (sel)
    for nn = 1:size(sel_AA,1) % = rows in the figure
        [aa_found, aa] = ismember(sel_AA{nn}, label_AA(:,1));   % aa = index to GRUP !!!
        assert(aa_found);
        assert(strcmp(label_AA{aa,1}, groupInfo.list_AA{aa,1}));

        % channels: MNI + VALS, cumulate over FB -> chnls_MNI_VAL
        if ~isempty(GRUP.chnls_MNI_all{aa,1})
            tmp_ch = GRUP.chnls_MNI_all{aa,1};               % all chnls MNI
            tmp_ch = cat(2, tmp_ch, nn*ones(size(tmp_ch,1),1));   % add value = 0 (not significant)
            chnls_MNI_VAL = cat(1, chnls_MNI_VAL, tmp_ch);      % all (non-significant) channels  
        end
        params.plot_brain3D.chVals_colorMap = cat(1, params.plot_brain3D.chVals_colorMap, chnls_clrs_nn(nn,:));
        
        % include Yeo7 NNs
%         if subj == 1
%             params.plot_brain3D.volumes2plot{end+1,1} = sel_AA{nn};
%             params.plot_brain3D.volumes2plot{end,2} = msSwitch_get_clr4aa(sel_AA{nn});
%             params.plot_brain3D.volumes2plot{end,3} = 0.4;
%             params.plot_brain3D.text_LR = true;
%         end
    end

    %% brain volumes -> surfaces
    [brainVols, clrmap, inds_vol] = brainVolumes_fixed(params, params.plot_brain3D);

    %% plot: 3D brains: different view angles
    for v = 1:size(sel_viewAngles,1)   

        % project channels MNI
        if strcmp(sel_viewAngles{v,1}, 'sagittal L')
            chnls_MNI_VAL(chnls_MNI_VAL(:,1)>=0,:) = [];  % only channels on L side
            chnls_MNI_VAL(:,1) = min(brainVols{1}.xi);   % all channels on L side - up front!
        elseif strcmp(sel_viewAngles{v,1}, 'axial')
            chnls_MNI_VAL(:,3) = max(brainVols{1}.zi);  % all channels on top - up front!
        elseif strcmp(sel_viewAngles{v,1}, 'sagittal R')
            chnls_MNI_VAL(chnls_MNI_VAL(:,1)<0,:) = [];  % only channels on R side
            chnls_MNI_VAL(:,1) = max(brainVols{1}.xi);   % all channels on R side - up front!
        end 

        % axes
        ax = subtightplot(nRows, nCols, nPlot, gap, marg_h, marg_w);
        hold on;

        % set: plot_brain3D & plotInfo
        plotInfo = params.plot_brain3D;
        plotInfo.axHandle = ax;
        plotInfo.figHandle = f;
        
        % PLOT 3D brain vols
        plotInfo.text_tag = subjList{subj,2};        
        plotInfo.chVals_lims = [0 2];
        for vol = 1:size(brainVols,2)
            if brainVols{vol}.loaded
                plotInfo.thisVolume = vol;
                % plot all chnl values at once: uncopmment line 137 in brain3D_plot.m
                h = brain3D_plot(brainVols{vol}, chnls_MNI_VAL, plotInfo);   
            end
        end    

        % title: FB
%         if v == 1
%             title(label_FB{freq,2}, 'FontSize',fontSize);  % , 'FontWeight','normal'
%         end

        % set view angle
        view(sel_viewAngles{v,2});

        nPlot = nPlot + 1;
    end % view angles

end

%% save
figName = ['subjBrains_v6_res0'];
outDir = [params.storage.dir_results filesep params.storage.outName filesep 'msSEI_supplFig_3D_brains'];
fig_save(f, figName, outDir, 'format','png', 'res',600);
% fig_save(f, figName, outDir, 'format','png');
% close(f); 
