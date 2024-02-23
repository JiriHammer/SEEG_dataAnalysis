%% ms SEI fig - anatomy area (AA): 3D brains + spectra + FBs, Yeo7 + Yeo17

% (c) Jiri, Oct23
% small known bugs:
% 1) FDR should be done for all P values at once (now fragmented in ch2roi_load)
% 2) bip. channel B8 can be assigned to DAN (B7-B8) & DMN (B8-B9)

%% settings: v1 - relaxed conditions (p-level, length of significance, ...)
dirSuffix = '_msFig_AA_v1_relaxCond';
t_sel = [-2.3, 2.8];    % the true sel time will be cut by 1 s on edges (= template length): [-2, 2] -> [-1, 1]
P_level = 0.05;         % significance level
t_siglen = 0.10;        % in [s], significance time threshold (ch is signif. if its significance > threshold) 
brainAtlas = 'Yeo17';

%% settings: v2 - relaxed conditions (p-level, length of significance, ...)
% dirSuffix = '_msFig_AA_v2_longSignif';
% t_sel = [-2.3, 2.8];    % the true sel time will be cut by 1 s on edges (= template length): [-2, 2] -> [-1, 1]
% P_level = 0.05;         % significance level
% t_siglen = 0.50;        % in [s], significance time threshold (ch is signif. if its significance > threshold) 
% brainAtlas = 'Yeo7';

%% settings: v3
% dirSuffix = '_msFig_AA_v3_strictCond';
% t_sel = [-1.0, 2.0];    % the true sel time will be cut by 1 s on edges (= template length): [-2, 2] -> [-1, 1]
% P_level = 0.01;         % significance level
% t_siglen = 0.25;        % in [s], significance time threshold (ch is signif. if its significance > threshold) 
% brainAtlas = 'Yeo7';

%% settings
load_ch2roi = true;             % load again the GRUP struct ?

%% get job-related settings (params): jobExecutor_SEI.m, v11 -> params
% --- specify analysis to be performed -> runner
runner = struct;
runner.newJobs = false; % set to 'false' if the results are already computed
runner.newFDR = false;  % true = updates the FDR correction over all patients (and freq. bands, ...) 

% --- job settings: specific (define spectral method) -> params
spectralMethod = 'v31_stft_session_bip';           % ~ v31 (BIP), similar to v12 (STFT base = sess, no zscore, no rej, FB = loFreq0)
% spectralMethod = 'v30_stft_baseRS_bip';           % ~ v30 (BIP), similar to v24 (new= no zscore, mo rej, FB = loFreq0), STFT with RS baseline
% spectralMethod = 'v29_stft_baseRS_bip';         % z-scored, loFreq0
% spectralMethod = 'v28_stft_baseRS_bip';         % ~ v24 (BIP), sparse BIP
% spectralMethod = 'v27_stft_baseRS_bip';         % z-scored, no loFreq0
% spectralMethod = 'v24_stft_baseRS_car';         % ~ v24 (CAR), STFT with RS baseline
% spectralMethod = 'v24_stft_baseRS_bip';         % ~ v24 (BIP), STFT with RS baseline
% spectralMethod = 'stft_sessions';       % ~ v12
% spectralMethod = 'stft_baseRS';         % ~ v11  (CAR)
% spectralMethod = 'stft_baseRS_bip';     % ~ v11b (BIP)

params = msSEI_getJobParams(runner, spectralMethod, P_level);

%% (un)comment: redefine subjList
% params.storage.subjList = msSEI_getSubjList('SEI_DTF_v31_5ch6fb');
% dirSuffix = '_msFig_AA_v2_subj_DTF_v31_5ch6fb';

%% =============================LOAD DATA==================================
%% ch2roi for Yeo7 -> GRUP 
if load_ch2roi
    clear GRUP;
    groupInfo = struct(...
        'brainAtlas', {{brainAtlas}}, ... % brain atlases: 'Yeo7','Yeo7_addLat','Yeo17','Yeo17_addLat','Mars','Mars_addLat'
        'freqBands', {{'delta','theta','alpha','beta','loGamma','hiGamma'}}, ...  % freq. bands (FB)
        'signifOfActiv', {{'was_sgnf'}}, ...   % significance of activation
        'significanceLevel', 0.05, ...      % significance level for ROI (e.g. DMN)
        'significanceTimeThr', t_siglen, ...    % in [s], significance time threshold (ch is signif. if its significance > threshold) 
        'GW_selection', {{'any'}}, ...      % grey / white matter: choices: 'gm','wm','csf','bone', 'any'
        'GW_threshold', 60, ...    % probability, in [%]
        'plotStyle', {{'allFB_allROI_subjMean','allFB_allROI','spectra_SNR'}}, ... 'oneFB_allROI','allFB_allROI','spectra_SNR'
        'maxDistToROI', 10, ...   % in [mm], maximum distance (tolerance) to region of interest
        'minSubjInROI', 1, ...  % minimum of subjects (patients) in a ROI (if a ROI has less patients, it is not considered in plots)
        'subtractAvgResponse', false, ...    % subtracts mean response over classes (conditions) from each clz
        'dirSuffix', dirSuffix, ... % adds suffix to 'ch2roi' directory name -> ch2roi_outDirSuffix
        'figSuffix', '', ... % adds suffix to figure name
        'time2plot', t_sel ... % in [s], w.r.t. cutting point (= 0)
        );
    ch2roi_wrapper;     % not that P = 0.001 is hardcoded in ch2roi_load.m (line: 192)
end

%% labels
label_FB = {
    'delta', '\delta: 0-3 Hz';
    'theta', '\theta: 4-7 Hz';
    'alpha', '\alpha: 8-12 Hz';
    'beta', '\beta: 13-30 Hz';
    'loGamma', 'LGB: 31-50 Hz';
    'hiGamma', 'HGB: 51-120 Hz';
    };

% label_AA_all = struct(...   % 1st = Yeo7NN, other = Yeo17 NNs
%     'VIS', {{'Visual', 'Visual peripheral', 'Visual central'}}, ...
%     'DAN', {{'Dorsal Attention', 'Dorsal attention A', 'Dorsal attention B'}}, ...
%     'SMT', {{'Somatomotor', 'Somatomotor A', 'Somatomotor B'}}, ...
%     'VAN', {{'Ventral Attention', 'Ventral attention', 'Salience'}}, ...
%     'FPN', {{'Frontoparietal', 'Control C', 'Control A', 'Control B'}}, ...
%     'LIM', {{'Limbic', 'Limbic N9', 'Limbic N10'}}, ...
%     'DMN', {{'Default', 'Default D (Auditory)', 'Default C', 'Default A', 'Default B'}} ...
%     );

label_AA_yeo7 = {
    'Visual',           'VIS';
    'Dorsal Attention', 'DAN';
    'Somatomotor',      'SMT';
    'Ventral Attention','VAN';
    'Frontoparietal',   'FPN';
    'Limbic',           'LIM';
    'Default',          'DMN';
    };

% Yeo17: all
label_AA_yeo17 = {
    'Visual peripheral',    'VIS P';
    'Visual central',       'VIS C';
    'Dorsal attention A',   'DAN A';
    'Dorsal attention B',   'DAN B';
    'Somatomotor A',        'SMT A';
    'Somatomotor B',        'SMT B';
    'Ventral attention',    'VAN';
    'Salience',             'SAL';
    'Control C',            'FPN C';
    'Control A',            'FPN A';
    'Control B',            'FPN B';
    'Limbic N9',            'LIM N9';
    'Limbic N10',           'LIM N10';
    'Default D (Auditory)', 'DMN D';
    'Default C',            'DMN C';
    'Default A',            'DMN A';
    'Default B',            'DMN B';
    };

label_AA = eval(['label_AA_' lower(brainAtlas)]);

%% ============================= plot each NN =============================
for nn = 1:size(label_AA,1)
%     sel_AA{1} = 'Default';
    sel_AA{1} = label_AA{nn,1};
    
    % view angles
    sel_viewAngles =  {
        'sagittal L', [-90, 0];          % sagittal R = 100,10 (90, 0)
        'axial',      [0, 90];          % axial = 0,90    
    %     'sagittal R', [100, 10];         % sagittal L = -100,10 (-90, 0)
        };

    % FIGURE
    nCols = size(GRUP.trials,2);   % FB
    nRows = size(sel_viewAngles,1) + 2*size(sel_AA,2);       % 2 for brains, for each AA: 1 = woSubtr, 2 = wiSubtr

    fig_W = 18;         % in [cm]
    fig_H = 3*nRows;    % in [cm], 3 cm per row
    f = fig_make('fig_W',25, 'fig_H',20);       % for Yeo7
    marg_h = [0.10 0.10];   % margin from [bottom, top]
    marg_w = [0.05 0.02];   % margin from [L, R] side
    gap = [-0.01, -0.01];     % between axes from [top, side]  

    fontSize = 10;

    % 3D brains of Yeo17 NNs
    params.plot_brain3D = struct(...        % user interface structure: holds most (but not all!) of the user settings, see also inside brain3D_plot.m
        'plot_brain3D', true, ...           % 3D brain model
        'volumes2plot', {{...               % volumes to plot (coregistered .nii files)
%             'colin27_t1_tal_lin masked left_h', [1 1 1], 0.20;  ...     % 1x3 cell: name, color, transparency (other choices for norm brain = wc1T1_colin27)
            'colin27_t1_tal_lin masked', [1 1 1], 0.20;  ...     % 1x3 cell: name, color, transparency (other choices for norm brain = wc1T1_colin27)
            }}, ...          
        'chVals_asCircles', true, ...       % plots channel values as color circles
        'chVals_asPatches', false, ...       % each channel value is mapped onto surface (but may be overwritten by latter points!), see also distanceThreshold for focality inside brain3D_plot.m
        'chVals_colorMap', [0 0 0; 1 0 0], ...    % colormap of channel values
        'chVals_patchTransp', 0.2, ...      % transparency of channel values
        'circleSizeLims', [2 4], ...       % in points, maps chnl values to circle size
        'size_interpolate', 1.0, ...        % in [mm], voxel size to which the brain is interpolated
        'backgroundColor', 'w', ...         % figure / axes background color, good choices: 'w', 'k'
        'text_LR', false, ...                % marks L / R side
        'text_tag', '', ...    % string
        'viewAngle', [0, 90], ...           % axial = 0,90; sagittal = 100,10; coronal front = 0,180; coronal back = 0,0
        'visible_axis', false ...           % if axis & their labels are visible
    ); 

    % include Yeo7 NNs
    for aa = 1:size(sel_AA,2)
        params.plot_brain3D.volumes2plot{end+1,1} = sel_AA{aa};
        params.plot_brain3D.volumes2plot{end,2} = msSwitch_get_clr4aa(sel_AA{aa});
        params.plot_brain3D.volumes2plot{end,3} = 0.4;
    end

    % ------------------ brain volumes -> surfaces ---------------------------
    [brainVols, clrmap, inds_vol] = brainVolumes_fixed(params, params.plot_brain3D);

    % plot: 3D brains
    nPlot = 1;
    for freq = 1:size(GRUP.trials,2)
        % plot different view angles
        for v = 1:size(sel_viewAngles,1)    
%             if strcmp(sel_viewAngles{v,1},'axial')
%                 params.plot_brain3D.volumes2plot{1,3} = 0.8;
%             else
%                 params.plot_brain3D.volumes2plot{1,3} = 0.2;
%             end

            % -----------------------channels: MNI + VALS--------------------------
            chnls_MNI_VAL = [];
            [aa_found, aa] = ismember(sel_AA{1}, label_AA(:,1));
            assert(aa_found);
            assert(strcmp(label_AA{nn,1}, groupInfo.list_AA{aa,1}));
            tmp_ch = GRUP.chnls_MNI_all{aa,freq};        % all chnls MNI
            tmp_ch = cat(2, tmp_ch, zeros(size(tmp_ch,1),1));   % add value = 0 (not significant)
            chnls_MNI_VAL = cat(1, chnls_MNI_VAL, tmp_ch);
            tmp_ch = GRUP.chnls_MNI_sel{aa,freq};        % sel chnls MNI
            tmp_ch = cat(2, tmp_ch, ones(size(tmp_ch,1),1));   % add value = 1 (was significant)
            if ~isempty(tmp_ch)
                chnls_MNI_VAL = cat(1, chnls_MNI_VAL, tmp_ch);    
            end
            if strcmp(sel_viewAngles{v,1}, 'sagittal L')
%                 chnls_MNI_VAL(:,1) = -abs(chnls_MNI_VAL(:,1));  % all channels on L side
                chnls_MNI_VAL(:,1) = min(brainVols{1}.xi);  % all channels on L side - up front!
            elseif strcmp(sel_viewAngles{v,1}, 'axial')
                chnls_MNI_VAL(:,3) = max(brainVols{1}.zi);  % all channels on top - up front!
            end

            % ---------------------------axes--------------------------------------
            ax = subtightplot(nRows, nCols, nPlot+(v-1)*nCols, gap, marg_h, marg_w);
            hold on;

            % ------------------------ set: plot_brain3D & plotInfo ---------------
            plotInfo = params.plot_brain3D;
            plotInfo.axHandle = ax;
            plotInfo.figHandle = f;

            % --------------------- PLOT 3D brain vols ------------------------
            for vol = 1:size(brainVols,2)
                if brainVols{vol}.loaded
                    plotInfo.thisVolume = vol;
                    % plot all chnl values at once: uncopmment line 137 in brain3D_plot.m
                    h = brain3D_plot(brainVols{vol}, chnls_MNI_VAL, plotInfo);   
                end
            end    

            % title: FB
            if v == 1
                title(label_FB{freq,2}, 'FontSize',fontSize);  % , 'FontWeight','normal'
            end

            % set view angle
            view(sel_viewAngles{v,2});

        end     % view angles
        nPlot = nPlot + 1;
    end     

    % ----------------------------------- plot FBs ----------------------------
    nPlot = size(sel_viewAngles,1)*nCols+1;
    clrs = GRUP.info.clzColor;

%     time2plot = [-1.3, 2.3];
    time2plot = t_sel;
    gap = [0.02, 0.01];     % between axes from [top, side]  
    for freq = 1:size(GRUP.trials,2)
        assert(strcmp(label_FB{freq,1}, groupInfo.list_FB{freq}));

        for subt = 1:2  % 1 = with subtraction, 2 = wo subtr.
            thisPlot = nPlot+(subt-1)*nCols;    % subplot below

            % axes
            subtightplot(nRows, nCols, thisPlot, gap, marg_h, marg_w);
            hold on;            
            set(gca, 'FontSize',fontSize);        

            % data to plot
            i_t = closestval(GRUP.time{aa,freq},time2plot(1)):closestval(GRUP.time{aa,freq},time2plot(2));
            x = GRUP.time{aa,freq}(i_t);
            y_avg = mean(GRUP.trials{aa,freq}(i_t,:,:),2);
            y_sem = sem(GRUP.trials{aa,freq}(i_t,:,:),2);
            if subt == 2
                y_avg = y_avg - mean(y_avg,3);  % subtract AVG (non-specific) response (mean over clz)
            end

            % >>> plot grouped channels activity in GRUP.trials{aa,freq} <<<
            for clz = 1:size(GRUP.trials{aa,freq},3)
%                 plotband(x, y_avg(:,:,clz), y_sem(:,:,clz), clrs(clz,:));
                msSEI_plotband(x, y_avg(:,:,clz), y_sem(:,:,clz), clz);
            end

            % plot the AVG (non-specific) response (mean over clz)
            if subt == 1
                y_nonSpec = mean(y_avg,3);  % mean over clz
                plot(x, y_nonSpec, '--k', 'LineWidth',1);
            end

            % plot baseline values
            if isfield(GRUP, 'BASE_vals') && subt == 1
                if ~isempty(GRUP.BASE_vals{aa,freq})
                    B_avg = GRUP.BASE_vals{aa,freq}(1);
                    B_sem = GRUP.BASE_vals{aa,freq}(2);
    %                     plotband(xlim, [B_avg,B_avg], [B_sem,B_sem], 'g');
                    plot(xlim, [B_avg,B_avg], 'Color',[0 0.5 0], 'LineWidth',2);
                end
            end

            % plot t=0 line
            axis tight;
            xlim(time2plot);
            plot([0 0], ylim, ':k', 'LineWidth',1);

            % plot significance of grouped activations
            if isfield(GRUP, 'H_trials')
                plot2axes_signif_filled(GRUP.time{aa,freq}(i_t), GRUP.H_trials{aa,freq}(i_t), GRUP.P_trials{aa,freq}(i_t));
            end

            % plot grouped paradigm times (RT, ...)
            if isfield(GRUP, 'paraTimes')
                plot2axes_paraTimes(GRUP.paraTimes{aa,freq}, clrs);
            end

            % ------------------------text, labels-----------------------------
            if subt == 1
                % increase y-axes to fit the text
                yLims = get(gca,'ylim');
                set(gca, 'ylim', [yLims(1), yLims(1) + 1.4*diff(yLims)]);   % increase by 10 %
                xLims = get(gca,'xlim');
                yLims = get(gca,'ylim');
                txt_x = xLims(1) + 0.10*diff(xLims);
                txt_y = yLims(1) + 0.85*diff(yLims);

                % number of channels, subjects, ratio
                N_chnl = [size(GRUP.trials{aa,freq},2), size(GRUP.trials{aa,freq},2)/GRUP.nChnls_roi(aa,freq), GRUP.nChnls_roi(aa,freq)];
                N_subj = [size(unique(GRUP.nSubj_roi{aa,freq,1}),1), size(unique(GRUP.nSubj_roi{aa,freq,1}),1)/size(unique(GRUP.nSubj_roi{aa,freq,2}(:,1)),1), size(unique(GRUP.nSubj_roi{aa,freq,2}(:,1)),1)]; 
                str_chInfo = ['R_C=' num2str(round(N_chnl(2)*100)) '%' , ... 
                            ', R_S=' num2str(round(N_subj(2)*100)) '%'];
                text(txt_x, txt_y, str_chInfo, 'FontSize',fontSize-1);
            end

            % xlabel
            if thisPlot > nCols*(nRows-1)
                xlabel('Time (s)', 'FontSize',fontSize);
            else
                set(gca, 'XTickLabel',{''});
            end
            if mod(thisPlot,nCols) == 1
    %             ylabel({label_AA{aa,2},['N_C=' num2str(N_chnl(3),3)], ['N_P=' num2str(N_subj(3),2)]}, 'FontSize',fontSize, 'FontWeight','bold');   % , 'Rotation',0, 'HorizontalAlignment','right', 'VerticalAlignment','center'
                ylabel('RBP (z-score)', 'FontSize',fontSize);
            end

            % Tick Labels
            set(gca, 'XTick',[-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5]);
            if thisPlot > nCols*(nRows-1)
                set(gca, 'XTickLabel',{'-2', '', '-1', '', '0', '', '1', '', '2', ''});        
            else
                set(gca, 'XTickLabel',{''});        
            end
            set(gca, 'YTick',[-1:0.025:1]);
            set(gca, 'YTickLabel',{''});

            box on;
            grid on;

        end
        nPlot = nPlot+1;
    end

    % text on upper part of the figure
    tx = axes('visible','off', 'position',[0 0 1 1]);
    mytitle = ['neural network: ' label_AA{nn,2} ' (' brainAtlas '), N_C = ' num2str(N_chnl(3)) ', N_P = ' num2str(N_subj(3))];
    % mytitle = strrep(mytitle, '_','\_');
    text(0.05, 0.93, mytitle, 'fontsize', fontSize, 'fontw', 'bold');
    
    % save
    figName = ['v13_P001_aa' num2str(nn) '_' label_AA{aa,2}];
%     figName = ['v8_P05_aa' num2str(nn) '_' label_AA{nn,2}];
    outDir = [params.storage.dir_results filesep params.storage.outName filesep 'ch2roi' groupInfo.dirSuffix filesep anatAtlas_name filesep 'msSEI_fig_AA'];
%     fig_save(f, figName, outDir, 'format','png', 'res',600);
    fig_save(f, figName, outDir, 'res',600);
    close(f); 
    
end