%% ms SEI - switching times as x-over time (XT) points on channels (avg over trials)
% low-pass filteres trials data
% determines "switching times" as cross-over points = t_switch
% plots histogram of t_switch

% (c) Jiri, Sep23
% based on: msSEI_switchingTimes_xover.m

%% ============================ settings ==================================
load_ch2roi = true;
lowPassTrials = true;
plotSingleChannels = false;

%% settings: v1 - relaxed conditions (p-level, length of significance, ...)
dirSuffix = '_msFig_XT_v1_relaxCond';
time4snr = [-2.3, 2.3];% time over which to compute SNR
t_sel = [-0.5, 2.0];    % time interval for X-over time point
P_level = 0.05;         % significance level
t_siglen = 0.10;        % in [s], significance time threshold (ch is signif. if its significance > threshold) 
brainAtlas = 'Yeo7';

%% settings: v2 - relaxed conditions + long time window 4 switching = not nice results
% dirSuffix = '_msFig_XT_v2_longWindow4XT';
% time4snr = [-2.3, 2.3];% time over which to compute SNR
% t_sel = [-2.0, 2.0];    % time interval for X-over time point
% P_level = 0.05;         % significance level
% t_siglen = 0.10;        % in [s], significance time threshold (ch is signif. if its significance > threshold) 
% brainAtlas = 'Yeo7';

%% settings: v3
% dirSuffix = '_msFig_AA_v3_strictCond';
% t_sel = [-1.0, 2.0];    % the true sel time will be cut by 1 s on edges (= template length): [-2, 2] -> [-1, 1]
% P_level = 0.01;         % significance level
% t_siglen = 0.25;        % in [s], significance time threshold (ch is signif. if its significance > threshold) 
% brainAtlas = 'Yeo7';

%% selected NN
sel_AA = {
    'Dorsal Attention';
    'Default';
};

%% freq. bands
label_FB = {
    'delta', '\delta: 0-3 Hz',  '\delta ';
    'theta', '\theta: 4-7 Hz',  '\theta ';
    'alpha', '\alpha: 8-12 Hz', '\alpha ';
    'beta', '\beta: 13-30 Hz',  '\beta ';
    'loGamma', 'LGB: 31-50 Hz', 'LGB';
    'hiGamma', 'HGB: 51-120 Hz','HGB';
    };

% label_FB = {
%     'loFreqs', 'LFB: 3-30 Hz',  'LFB';
%     'loGamma', 'LGB: 31-50 Hz', 'LGB';
%     'hiGamma', 'HGB: 51-120 Hz','HGB';
%     };

% label_FB = {
%     'loFreq0', 'LFB: 0-30 Hz',  'LFB';
%     'loGamma', 'LGB: 31-50 Hz', 'LGB';
%     'hiGamma', 'HGB: 51-120 Hz','HGB';
%     };

%% settings: AA -> label_AA
label_AA = {
    'Visual',           'VIS';
    'Dorsal Attention', 'DAN';
    'Somatomotor',      'SMT';
    'Ventral Attention','VAN';
    'Frontoparietal',   'FPN';
    'Limbic',           'LIM';
    'Default',          'DMN';
    };

%% get job-related settings -> params
% --- specify analysis to be performed -> runner
runner = struct;
runner.newJobs = false; % set to 'false' if the results are already computed
runner.newFDR = true;  % true = updates the FDR correction over all patients (and freq. bands, ...) 

% --- job settings: specific (define spectral method) -> params
spectralMethod = 'v31_stft_session_bip';           % ~ v31 (BIP), similar to v12 (STFT base = sess, no zscore, no rej, FB = loFreq0)
% spectralMethod = 'v30_stft_baseRS_bip';         % ~ v30 (BIP), STFT with RS baseline
% spectralMethod = 'v24_stft_baseRS_car';         % ~ v24 (CAR), STFT with RS baseline
% spectralMethod = 'v24_stft_baseRS_bip';         % ~ v24 (BIP), STFT with RS baseline
% spectralMethod = 'stft_sessions';       % ~ v12
% spectralMethod = 'stft_baseRS';         % ~ v11  (CAR)
% spectralMethod = 'stft_baseRS_bip';     % ~ v11b (BIP)

params = msSEI_getJobParams(runner, spectralMethod, P_level);

%% ch2roi: get only significant activations for Yeo7 & all bands -> GRUP 
if load_ch2roi
    groupInfo = struct(...
        'brainAtlas', {{brainAtlas}}, ... % brain atlases: 'Yeo7','Yeo7_addLat','Yeo17','Yeo17_addLat','Mars','Mars_addLat'
        'freqBands', {label_FB(:,1)'}, ...  % freq. bands (FB)
        'signifOfActiv', {{'was_sgnf'}}, ...   % significance of activation
        'significanceLevel', P_level, ...      % significance level for ROI (e.g. DMN)
        'significanceTimeThr', t_siglen, ...    % in [s], significance time threshold (ch is signif. if its significance > threshold) 
        'GW_selection', {{'any'}}, ...      % grey / white matter: choices: 'gm','wm','csf','bone', 'any'
        'GW_threshold', 60, ...    % probability, in [%]
        'plotStyle', {{'oneFB_allROI','allFB_allROI','spectra_SNR','allFB_allROI_SNR'}}, ... 'oneFB_allROI','allFB_allROI','spectra_SNR'
        'maxDistToROI', 10, ...   % in [mm], maximum distance (tolerance) to region of interest
        'minSubjInROI', 1, ...  % minimum of subjects (patients) in a ROI (if a ROI has less patients, it is not considered in plots)
        'subtractAvgResponse', true, ...    % subtracts mean response over classes (conditions) from each clz
        'dirSuffix', dirSuffix, ... % adds suffix to 'ch2roi' directory name -> ch2roi_outDirSuffix
        'figSuffix', '', ... % adds suffix to figure name
        'time2plot', time4snr ... % in [s], w.r.t. cutting point (= 0)
        );
    groupInfo.file_brainClusters = 'C:\Users\Jiri\Documents\MATLAB\brainClusters_tables\LabelsAK_1.xlsx';
    ch2roi_wrapper;
end

%% >>>>>>> SNR data: 2 groups (= AA: DAN & DMN) of 6 params (= FB) <<<<<<<<
nFreq = size(GRUP.trials,2);
SNR_avg = nan(size(sel_AA,1),nFreq);    % 2D: nn x freq
SNR_sem = nan(size(sel_AA,1),nFreq);    % 2D: nn x freq
P_snr = nan(nFreq,nFreq,size(sel_AA,1));
for nn = 1:size(sel_AA,1) % = rows in the figure
    [aa_found, aa] = ismember(sel_AA{nn}, label_AA(:,1));   % aa = index to GRUP !!!
    assert(aa_found);
    assert(strcmp(label_AA{aa,1}, groupInfo.list_AA{aa,1}));
    
    for freq = 1:nFreq
        assert(strcmp(label_FB{freq,1}, groupInfo.list_FB{freq}));
        i_t = closestval(GRUP.time{aa,freq},time4snr(1)):closestval(GRUP.time{aa,freq},time4snr(2));
        
        % AVG +/- SEM
        SNR_avg(nn,freq) = mean(sum(GRUP.SNR_trials{aa,freq}(i_t,:),1)./diff(time4snr),2);
        SNR_sem(nn,freq) = sem(sum(GRUP.SNR_trials{aa,freq}(i_t,:),1)./diff(time4snr),2);
        
        % significance (over chnls in each aa)
        for fb2 = 1:nFreq
            x_fb1 = sum(GRUP.SNR_trials{aa,freq}(i_t,:),1)./diff(time4snr);
            x_fb2 = sum(GRUP.SNR_trials{aa,fb2}(i_t,:),1)./diff(time4snr);
            P_snr(freq,fb2,nn) = ranksum(x_fb1,x_fb2, 'alpha',P_level);
        end
    end
end
            
% FDR correction -> H_fdr -> H_matrix
P_vals = P_snr(:);
% [N_fdr, i_fdr] = fdr(P_vals, P_level, 'general');
[N_fdr, i_fdr] = fdr(P_vals, 0.001, 'general');
H_fdr = zeros(size(P_vals,1),1);    % default = 1 -> not significant
if N_fdr > 1
    H_fdr(i_fdr) = 1;
end
H_snr = reshape(H_fdr, [size(P_snr,1), size(P_snr,2), size(P_snr,3)]);

%% >>>>>>>>>> TX: switching times from xover points -> SWITCH_lags{aa,freq} <<<<<<<<<<
timeAxis = GRUP.time{1,1};
i_t = closestval(timeAxis, t_sel(1)):closestval(timeAxis, t_sel(2));
nAA = size(sel_AA,1);
SWITCH_lags = cell(nAA, size(GRUP.trials,2));  % 2D cell: aa x freq, times of switching
SWITCH_subj = cell(nAA, size(GRUP.trials,2));  % 2D cell: aa x freq, subject tags (how many different subject were in each ROI)
for nn = 1:nAA
    [tf, aa] = ismember(sel_AA{nn}, groupInfo.list_AA(:,1));
    for freq = 1:size(GRUP.trials,2)
        
        % low-pass filter data: GRUP.trials{aa,freq} -> filtData
        if lowPassTrials
            fs = 1/mean(diff(GRUP.time{aa,freq}));      % in [Hz], computed from time vector
            freqNyquist = fs/2;                         % in [Hz]
            loF = 0;                                    % in [Hz]
            hiF = 1;                                    % in [Hz]
            Wn = hiF/freqNyquist;                       % normalized bandpass frequencies
            n = 3;                                      % butterworth order
            [b,a] = butter(n, Wn, 'low');               % returns polynoms of Butterw. filter
            filtData = filtfilt(b, a, GRUP.trials{aa,freq});
        else
            filtData = GRUP.trials{aa,freq};
        end
        
        % subtract mean over clz (non-spec. response) -> D_trials
        D_trials = filtData - mean(filtData,3);     
        GRUP.trials{aa,freq} = GRUP.trials{aa,freq} - mean(GRUP.trials{aa,freq},3);  
        
        % time of switching
        for ch = 1:size(GRUP.trials{aa,freq},2)
            
            % ---------------time of switching -> t_switch-----------------
            t_switch = msSEI_getSwitchingTime(D_trials(i_t,ch,1), timeAxis(i_t));
            if ~isempty(t_switch)
                SWITCH_lags{nn,freq} = cat(2, SWITCH_lags{nn,freq}, t_switch);        % append    
                SWITCH_subj{nn,freq} = cat(2, SWITCH_subj{nn,freq}, GRUP.nSubj_roi{aa,freq,1}(ch));        % append    
            end
            
            % --------------------------- plot ----------------------------
            if plotSingleChannels 
                f = fig_make;
                hold on;

                for clz = 1:size(D_trials,3)
                    plot(timeAxis, GRUP.trials{aa,freq}(:,ch,clz), 'Color',GRUP.info.clzColor(clz,:), 'LineStyle','-');
                    plot(timeAxis, D_trials(:,ch,clz), 'Color',GRUP.info.clzColor(clz,:), 'LineWidth',2, 'LineStyle','-', 'Marker','.');
                end

                % mark significances
                for s = 1:size(GRUP.sel_hVals{aa,freq},1)
                    if GRUP.sel_hVals{aa,freq}(s,ch) == 1
                        plot(timeAxis(s), 0, 'Marker','*', 'Color','k');
                    end
                end        
                if ~isempty(t_switch)
                    plot([t_switch, t_switch], ylim, '--m', 'LineWidth',2);
                end
                
                % plot sel time
                plot([t_sel(1), t_sel(1)], ylim, '--k');
                plot([t_sel(2), t_sel(2)], ylim, '--k');
                grid on;

                title(['ROI = ' list_AA{aa} ', FB = ' groupInfo.list_FB{freq}]);
                xlabel('time (s)');
                ylabel([groupInfo.list_FB{freq} ' zscore']);
                
                % text on upper part of the figure
                if isfield(GRUP.info, 'txtOnFig')
                    tx = axes('visible','off', 'position',[0 0 1 1]);
                    subj = GRUP.nSubj_roi{aa,freq,1}(ch);
                    assert(strcmp(GRUP.chSubj_name{aa,freq,1}{ch,2}, params.storage.subjList{subj,1}));
%                     mytitle = ['subj = ' params.storage.subjList{subj,1} ', ch = ' num2str(GRUP.sel_chNums{aa,freq}(ch)) ', ' GRUP.info.txtOnFig];
                    mytitle = ['subj = ' params.storage.subjList{subj,1} ', ch = ' GRUP.chSubj_name{aa,freq,1}{ch,1} ', ' GRUP.info.txtOnFig];
                    mytitle = strrep(mytitle, '_','\_');
                    text(0.016, 0.98, mytitle, 'fontsize', 14, 'fontw', 'bold');
                end

                % save
                figName = ['ROI_' num2str(aa) '_FB_' num2str(freq) '_ch' num2str(ch) '_' list_AA{aa} '_' groupInfo.list_FB{freq}];
                outDir = [params.storage.dir_results filesep params.storage.outName filesep 'ch2roi' groupInfo.dirSuffix filesep anatAtlas_name filesep 'xover_singleCh' filesep list_AA{aa} filesep groupInfo.list_FB{freq}];
                fig_save(f, figName, outDir, 'format','png');
                close(f); 
            end
            % -------------------------------------------------------------
            
        end
        disp(['Switching points: aa = ' num2str(aa) ', freq = ' num2str(freq) ' - done!'])
    end
end

%% TX: significance of difference in distribution: AA - to - AA for each FB
% compute significances -> P_vals vector of all pairs of comparisons
nAA = size(SWITCH_lags,1);
nFreq = size(SWITCH_lags,2);
P_matrix = nan(nAA,nAA,nFreq); % 3D: aa x aa x freq

% ranksum test: all-2-all for significance
for freq = 1:nFreq
    labels_aa = [];
    for aa1 = 1:size(SWITCH_lags,1)
        for aa2 = 1:size(SWITCH_lags,1) 
            if ~isempty(SWITCH_lags{aa1,freq}) && ~isempty(SWITCH_lags{aa2,freq})
                P_matrix(aa1,aa2,freq) = ranksum(SWITCH_lags{aa1,freq},SWITCH_lags{aa2,freq}, 'alpha',P_level);
            else
                P_matrix(aa1,aa2,freq) = 1;
            end
        end
%         labels_aa{aa1} = [list_AA{aa1,1}];
    end
end

% FDR correction -> H_fdr -> H_matrix
P_vals = P_matrix(:);
% [N_fdr, i_fdr] = fdr(P_vals, P_level, 'general');
[N_fdr, i_fdr] = fdr(P_vals, 0.001, 'general');
H_fdr = zeros(size(P_vals,1),1);    % default = 1 -> not significant
if N_fdr > 1
    H_fdr(i_fdr) = 1;
end
H_matrix = reshape(H_fdr, [size(P_matrix,1), size(P_matrix,2), size(P_matrix,3)]);

%% =========================== FIGURE =====================================
fig_W = 18;         % in [cm]
fig_H = 12;    % in [cm], 3 cm per row
f = fig_make('fig_W',fig_W, 'fig_H',fig_H);
marg_h = [0.10 0.10];   % margin from [bottom, top]
marg_w = [0.10 0.02];   % margin from [L, R] side
gap = [0.05, 0.08];     % between axes from [top, side]  

nRows = 1;
nCols = 2;
nPlot = 1;

fontSize = 10;

% ============ A: barplot - mean SNR over *channels in GRUP.SNR_trials ========
% axes
ax = subtightplot(nRows, nCols, nPlot, gap, marg_h, marg_w);
hold on;

[hBar, hErrorBar] = barwitherr(SNR_sem, SNR_avg);

set(gca,'XTick',[1,2]);
set(gca,'XTickLabel',{'DAN','DMN'}, 'FontSize',fontSize);
legend(hBar, label_FB(:,2), 'location','NorthEast', 'FontSize',fontSize);

ylabel('E-I vs. I-E difference strength (SNR)', 'FontSize',fontSize);
    
% ============ B: boxplot - crossing time (TX) ========
% axes
nPlot = nPlot+1;
ax = subtightplot(nRows, nCols, nPlot, gap, marg_h, marg_w);
hold on;

x = [1:size(SWITCH_lags,2); 1:size(SWITCH_lags,2)];
x(1,:) = x(1,:)-0.15;
x(2,:) = x(2,:)+0.15;
xtick_label = cell(size(SWITCH_lags,2),1);
for freq = 1:size(SWITCH_lags,2)
    for aa = 1:size(SWITCH_lags,1)
        clr_aa = msSwitch_get_clr4aa(sel_AA{aa});
%         h = boxplot(SWITCH_lags{aa,freq}, 'Positions',x(aa,freq), 'Symbol','w', 'Widths',0.15, 'Colors',clr_aa); 
%         h = boxplot(SWITCH_lags{aa,freq}, 'Positions',x(aa,freq), 'Symbol','w', 'Widths',0.15, 'whisker',1.5);  % whisker = full range of data
        h = boxplot(SWITCH_lags{aa,freq}, 'Positions',x(aa,freq), 'Symbol','w', 'Widths',0.15, 'whisker',0.4);  % whisker = 0.4 => 5th & 95th prctile range of data
        set(h(5), 'Color', clr_aa);
        set(h(5), 'LineWidth', 1);
        xtick_label{freq} = [label_FB{freq,3}, ' (N_C = ' num2str(size(SWITCH_lags{aa,freq},2)) ', N_P = ' num2str(size(unique(SWITCH_subj{aa,freq}),2)) ')'];
    end
    
    % significance
    if H_matrix(1,2,freq) == 1
        P = 0;
    else
        P = 1;
    end
    CNN_makeSignificanceBar([x(1,freq), x(2,freq)], 1.5, P, P_level);
end
set(gca,'xlim',[0 x(1,end)+1.0]); 
set(gca,'XTick',1:size(SWITCH_lags,2));
ax.TickLabelInterpreter = 'tex';
set(gca,'XTickLabel',label_FB(:,3), 'FontSize',fontSize);
% set(gca,'XTickLabel',xtick_label, 'FontSize',fontSize);
% set(gca,'XTickLabelRotation',35);
ylim([-0.55, 1.7]);
ylabel('Switching time (s) between E-I & I-E', 'FontSize',fontSize);
plot(xlim, [0 0], ':k');

grid on;
box on;

%% save
% figName = 'TX_DMN_DAN_v3_loFreq0';
figName = 'TX_DMN_DAN_v6_7FB';
outDir = [params.storage.dir_results filesep params.storage.outName filesep 'ch2roi' groupInfo.dirSuffix filesep anatAtlas_name filesep 'msSEI_fig_SNR_TX'];
% fig_save(f, figName, outDir, 'format','png', 'res',600);
fig_save(f, figName, outDir, 'res',600);
% close(f); 
