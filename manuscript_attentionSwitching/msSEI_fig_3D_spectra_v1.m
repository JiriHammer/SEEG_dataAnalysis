%% ms SEI fig - 3D brains + spectra (E-I, I-E)
% rows = NN: DMN & DAN
% cols:
%   - 3D brains 
%   - spectra (LFB, HBG)

% (c) Jiri, Jan24

%% settings: v1 - relaxed conditions (p-level, length of significance, ...)
dirSuffix = '_msSEI_fig_3D_spectra_v1';
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

%% ch2roi for Yeo7 -> GRUP 
if load_ch2roi
    clear GRUP;
    groupInfo = struct(...
        'brainAtlas', {{brainAtlas}}, ... % brain atlases: 'Yeo7','Yeo7_addLat','Yeo17','Yeo17_addLat','Mars','Mars_addLat'
        'freqBands', {{'hiGamma'}}, ...  % freq. bands (FB)
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

%% settings: AA -> label_AA
label_AA = msSEI_get_label_AA(lower(brainAtlas));

sel_AA = {
    'Default';
    'Dorsal Attention';
};

%% 3D brain params & view angles -> sel_viewAngles
sel_viewAngles =  {
    'sagittal L', [-90, 0];          % sagittal R = 100,10 (90, 0)
    'axial',      [0, 90];          % axial = 0,90    
    'sagittal R', [90, 0];         % sagittal L = -100,10 (-90, 0)
    };
    
%% ========================= FIGURE =======================================
nCols = size(sel_viewAngles,1) + 2;   % views + 2 spectra
nRows = size(sel_AA,1);       % 2 for brains, for each AA: 1 = woSubtr, 2 = wiSubtr

fig_W = 18;         % in [cm]
fig_H = 5*nRows;    % in [cm], 3 cm per row
f = fig_make('fig_W',fig_W, 'fig_H',fig_H);
marg_h = [0.15 0.10];   % margin from [bottom, top]
marg_w = [0.02 0.10];   % margin from [L, R] side
gap = [-0.01, -0.01];     % between axes from [top, side]  

fontSize = 10;
time2plot = [-2.3, 2.3];
freq2plot = [0, 120];

label_clz = {
    'E-I','External - Internal';
    'I-E','Internal - External';
};
cLims_nn = [...
    -0.05, 0.05;
    -0.08, 0.08;
    ];

% =============================== PLOT ===================================
nPlot = 1;
for nn = 1:size(sel_AA,1) % = rows in the figure
    [aa_found, aa] = ismember(sel_AA{nn}, label_AA(:,1));   % aa = index to GRUP !!!
    assert(aa_found);
    assert(strcmp(label_AA{aa,1}, groupInfo.list_AA{aa,1}));

    % =========================== PLOT: 3D brain (all chnls) =================
    % settings
    params.plot_brain3D = struct(...        % user interface structure: holds most (but not all!) of the user settings, see also inside brain3D_plot.m
        'plot_brain3D', true, ...           % 3D brain model
        'volumes2plot', {{...               % volumes to plot (coregistered .nii files)
            'colin27_t1_tal_lin masked', [1 1 1], 0.20;  ...     % 1x3 cell: name, color, transparency (other choices for norm brain = wc1T1_colin27)
            }}, ...          
        'chVals_asCircles', true, ...       % plots channel values as color circles
        'chVals_asPatches', false, ...       % each channel value is mapped onto surface (but may be overwritten by latter points!), see also distanceThreshold for focality inside brain3D_plot.m
        'chVals_colorMap', [0 0 0; 0 0 0], ...    % colormap of channel values
        'chVals_patchTransp', 0.2, ...      % transparency of channel values
        'circleSizeLims', [3 3], ...       % in points, maps chnl values to circle size
        'size_interpolate', 1.0, ...        % in [mm], voxel size to which the brain is interpolated
        'backgroundColor', 'w', ...         % figure / axes background color, good choices: 'w', 'k'
        'text_LR', false, ...                % marks L / R side
        'text_tag', '', ...    % string
        'viewAngle', [0, 90], ...           % axial = 0,90; sagittal = 100,10; coronal front = 0,180; coronal back = 0,0
        'visible_axis', false ...           % if axis & their labels are visible
    ); 

    % axes settings
    gap = [-0.01, 0.00];     % between axes from [top, side]  
    marg_w = [0.02 0.12];   % margin from [L, R] side
%     yyaxis left;
    
    % include Yeo7 NNs
    params.plot_brain3D.volumes2plot{end+1,1} = sel_AA{nn};
    params.plot_brain3D.volumes2plot{end,2} = msSwitch_get_clr4aa(sel_AA{nn});
    params.plot_brain3D.volumes2plot{end,3} = 0.4;

    % brain volumes -> surfaces
    [brainVols, clrmap, inds_vol] = brainVolumes_fixed(params, params.plot_brain3D);

    % plot: 3D brains: different view angles
    for v = 1:size(sel_viewAngles,1)   

        % channels: MNI + VALS, cumulate over FB -> chnls_MNI_VAL
        chnls_MNI_VAL = [];    
        tmp_ch = GRUP.chnls_MNI_all{aa,1};               % all chnls MNI
        tmp_ch = cat(2, tmp_ch, zeros(size(tmp_ch,1),1));   % add value = 0 (not significant)
        chnls_MNI_VAL = cat(1, chnls_MNI_VAL, tmp_ch);      % all (non-significant) channels     
%         for freq = 1:size(GRUP.trials,2)
%             tmp_ch = GRUP.chnls_MNI_sel{aa,freq};           % sel chnls MNI
%             tmp_ch = cat(2, tmp_ch, ones(size(tmp_ch,1),1));% add value = 1 (was significant)
%             chnls_MNI_VAL = cat(1, chnls_MNI_VAL, tmp_ch);  % add (again) the significant channels, (note that they are duplicate)
%         end
%         if strcmp(sel_viewAngles{v,1}, 'sagittal L')
% %                 chnls_MNI_VAL(:,1) = -abs(chnls_MNI_VAL(:,1));  % all channels on L side
%             chnls_MNI_VAL(:,1) = min(brainVols{1}.xi);  % all channels on L side - up front!
%         elseif strcmp(sel_viewAngles{v,1}, 'axial')
%             chnls_MNI_VAL(:,3) = max(brainVols{1}.zi);  % all channels on top - up front!
%         end

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
        plotInfo.chVals_lims = [0 1];
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

    % =========================== PLOT: spectra ===========================
    % axes settings
    gap = [0.05, 0.02];     % between axes from [top, side] 
    marg_w = [0.11 0.10];   % margin from [L, R] side
    for clz = 1:size(GRUP.spectra_data{aa},3)

        % axes
        ax = subtightplot(nRows, nCols, nPlot, gap, marg_h, marg_w);
        hold on;            
        set(gca, 'FontSize',fontSize);        

        % data to plot
        i_t = closestval(GRUP.spectra_time,time2plot(1)):closestval(GRUP.spectra_time,time2plot(2));
        i_f = closestval(GRUP.spectra_freq,freq2plot(1)):closestval(GRUP.spectra_freq,freq2plot(2));
        xVals = GRUP.spectra_time(i_t);
        yVals = GRUP.spectra_freq(i_f);
        cVals = GRUP.spectra_data{aa}(i_f,i_t,clz);
        cLims = cLims_nn(nn,:);
        
        % >>> plot grouped channels activity in GRUP.trials{aa,freq} <<<
        h = imagesc(xVals, yVals, cVals, cLims);
        colormap(gca,brewermap(256,'*RdBu'));   % set colormap for this axes (independent of figure colormap)

        % axes props
        axis tight;

        % plot t=0 line
        axis tight;
        xlim(time2plot);
        plot([0 0], ylim, ':k', 'LineWidth',2);
        freqVals = [13, 30, 50, 100];
        for fb = freqVals
            plot(xlim, [fb fb], '--k', 'LineWidth',0.5);
        end
        timeVals = [-2:0.5:-0.5, 0.5:0.5:2];
        for tt = timeVals
            plot([tt tt], ylim, '--k', 'LineWidth',0.5);
        end
        
        % text
%         % increase y-axes to fit the text
%         yLims = get(gca,'ylim');
%         set(gca, 'ylim', [yLims(1), yLims(1) + 1.2*diff(yLims)]);   % increase by 10 %
%         xLims = get(gca,'xlim');
%         yLims = get(gca,'ylim');
%         txt_x = xLims(1) + 0.20*diff(xLims);
%         txt_y = yLims(1) + 0.95*diff(yLims);
% 
%         % number of channels, subjects, ratio
        freq = 1;
        N_chnl = [size(GRUP.trials{aa,freq},2), size(GRUP.trials{aa,freq},2)/GRUP.nChnls_roi(aa,freq), GRUP.nChnls_roi(aa,freq)];
        N_subj = [size(unique(GRUP.chSubj_inds{aa,freq,1}(:,2)),1), size(unique(GRUP.chSubj_inds{aa,freq,1}(:,2)),1)/size(unique(GRUP.chSubj_inds{aa,freq,2}(:,2)),1), size(unique(GRUP.chSubj_inds{aa,freq,2}(:,2)),1)]; 
        disp([sel_AA{nn} ': nCh = ' num2str(N_chnl(1)) ', nSubj = ' num2str(N_subj(1))]);
%         str_chInfo = ['R_C=' num2str(round(N_chnl(2)*100)) '%' , ... 
%                     ', R_P=' num2str(round(N_subj(2)*100)) '%'];
%         text(txt_x, txt_y, str_chInfo, 'FontSize',fontSize-1);

        % labels
        if nPlot > nCols*(nRows-1)
            xlabel('Time (s)', 'FontSize',fontSize);
        else
            set(gca, 'XTickLabel',{''});
        end
        set(gca, 'YTick',[0, freqVals, 120]);
        if mod(nPlot,nCols) == size(sel_viewAngles,1)+1
%             ylabel({label_AA{aa,2},['N_C=' num2str(N_chnl(3),3)], ['N_P=' num2str(N_subj(3),2)]}, 'FontSize',fontSize, 'FontWeight','bold');   % , 'Rotation',0, 'HorizontalAlignment','right', 'VerticalAlignment','center'
            ylabel('Frequency (Hz)', 'FontSize',fontSize);
            set(gca, 'YTickLabel',cat(2, {0}, num2cell(freqVals), {120}));
        else
            set(gca, 'YTickLabel',{''});
        end
%         title(label_clz{clz,2}, 'FontSize',fontSize);

        % Tick Labels
        set(gca, 'XTick',[-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2]);
        if nPlot > nCols*(nRows-1)
            set(gca, 'XTickLabel',{'-2', '', '-1', '', '0', '', '1', '', '2'});        
        else
            set(gca, 'XTickLabel',{''});        
        end
        
        % colorbar
        if mod(nPlot,nCols) == 0 
            hCb = colorbar;
            ylabel(hCb, 'PSD (z-score)', 'FontSize',10, 'FontWeight','normal');
            set(hCb, 'position', [ax.Position(1)+ax.Position(3)+0.01, ax.Position(2), 0.015, ax.Position(4)]);
        end

        nPlot = nPlot+1;
    end  % of clz   
    
end  % of nn

%% save
figName = ['anySignif_v4'];
outDir = [params.storage.dir_results filesep params.storage.outName filesep 'msSEI_fig_3D_spectra_v1'];
% fig_save(f, figName, outDir, 'format','png', 'res',600);
fig_save(f, figName, outDir, 'res',600);
% close(f); 
