%% ms SEI - switching times as x-over time (XT) points on channels (avg over trials)
% low-pass filters trials data
% determines "switching times" as cross-over points = t_switch
% plots histogram of t_switch

% (c) Jiri, Sep23
% based on: msSEI_switchingTimes_xover.m

%% ============================ settings ==================================
load_ch2roi = true;
lowPassTrials = true;
plotSingleChannels = true;

%% settings: v1 - relaxed conditions (p-level, length of significance, ...)
dirSuffix = '_msFig_suppl_XT_sCh';
time4snr = [-2.3, 2.3];% time over which to compute SNR
t_sel = [-0.5, 2.0];    % time interval for X-over time point
P_level = 0.05;         % significance level
t_siglen = 0.10;        % in [s], significance time threshold (ch is signif. if its significance > threshold) 
brainAtlas = 'Yeo7';

%% selected NN
sel_AA = {
%     'Dorsal Attention';
    'Default';
};

%% freq. bands
label_FB = {
%     'delta', '\delta: 0-3 Hz',  '\delta ';
%     'theta', '\theta: 4-7 Hz',  '\theta ';
%     'alpha', '\alpha: 8-12 Hz', '\alpha ';
%     'beta', '\beta: 13-30 Hz',  '\beta ';
%     'loGamma', 'LGB: 31-50 Hz', 'LGB';
    'hiGamma', 'HGB: 51-120 Hz','HGB';
    };

%% settings: AA -> label_AA
label_AA = {
%     'Visual',           'VIS';
%     'Dorsal Attention', 'DAN';
%     'Somatomotor',      'SMT';
%     'Ventral Attention','VAN';
%     'Frontoparietal',   'FPN';
%     'Limbic',           'LIM';
    'Default',          'DMN';
    };

%% get job-related settings -> params
% --- specify analysis to be performed -> runner
runner = struct;
runner.newJobs = false; % set to 'false' if the results are already computed
runner.newFDR = false;  % true = updates the FDR correction over all patients (and freq. bands, ...) 

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
                %% figure
                fig_W = 12;         % in [cm]
                fig_H = 10;         % in [cm]
                f = fig_make('fig_W',fig_W, 'fig_H',fig_H);                
%                 f = fig_make;
                hold on;
                
                % plot sCh trial-mean & lowPass
                clrs = {
                    [0 0 255]./255, [255 0 255]./255;   % dark b -> m
                    [255 0 0]./255, [81,174,255]./255;  % r -> light b
                };
                i_zero = closestval(timeAxis, 0.0);    % index at t = 0
                for clz = 1:size(D_trials,3)
                    y_avg = GRUP.trials{aa,freq}(:,ch,clz);
                    plot(timeAxis(1:i_zero),    y_avg(1:i_zero),   'Color',clrs{clz,1});    % plot from beg -> 0 (= task switch)
                    plot(timeAxis(i_zero:end),  y_avg(i_zero:end), 'Color',clrs{clz,2});     % plot from 0 -> end
                    
                    y_avg = D_trials(:,ch,clz);
                    plot(timeAxis(1:i_zero),    y_avg(1:i_zero),   'Color',clrs{clz,1}, 'LineWidth',2);    % plot from beg -> 0 (= task switch)
                    plot(timeAxis(i_zero:end),  y_avg(i_zero:end), 'Color',clrs{clz,2}, 'LineWidth',2);     % plot from 0 -> end                    
%                     plot(timeAxis, GRUP.trials{aa,freq}(:,ch,clz), 'Color',GRUP.info.clzColor(clz,:), 'LineStyle','-');
%                     plot(timeAxis, D_trials(:,ch,clz), 'Color',GRUP.info.clzColor(clz,:), 'LineWidth',2, 'LineStyle','-', 'Marker','.');
                end
                plot([0, 0], ylim, ':k', 'LineWidth',2);
                
                % mark significances
                pVals = 0.0005*ones(size(GRUP.sel_hVals{aa,freq},1),1); % fake pVals
                plot2axes_signif_filled(timeAxis, GRUP.sel_hVals{aa,freq}(:,ch), pVals)
%                 for s = 1:size(GRUP.sel_hVals{aa,freq},1)
%                     if GRUP.sel_hVals{aa,freq}(s,ch) == 1
%                         plot(timeAxis(s), 0, 'Marker','*', 'Color','k');
%                     end
%                 end        
                if ~isempty(t_switch)
                    plot([t_switch, t_switch], ylim, '--m', 'LineWidth',2);
                end
                
                % plot sel time
                plot([t_sel(1), t_sel(1)], ylim, '--k');
                plot([t_sel(2), t_sel(2)], ylim, '--k');
                grid on;

%                 title(['ROI = ' list_AA{aa} ', FB = ' groupInfo.list_FB{freq}]);
                xlabel('Time (s)');
                ylabel([label_FB{freq,2} ' (z-score)']);
                xlim([-2.3, 2.3]);
                
                % text on upper part of the figure
%                 if isfield(GRUP.info, 'txtOnFig')
%                     tx = axes('visible','off', 'position',[0 0 1 1]);
%                     subj = GRUP.nSubj_roi{aa,freq,1}(ch);
%                     assert(strcmp(GRUP.chSubj_name{aa,freq,1}{ch,2}, params.storage.subjList{subj,1}));
% %                     mytitle = ['subj = ' params.storage.subjList{subj,1} ', ch = ' num2str(GRUP.sel_chNums{aa,freq}(ch)) ', ' GRUP.info.txtOnFig];
%                     mytitle = ['subj = ' params.storage.subjList{subj,1} ', ch = ' GRUP.chSubj_name{aa,freq,1}{ch,1} ', ' GRUP.info.txtOnFig];
%                     mytitle = strrep(mytitle, '_','\_');
%                     text(0.016, 0.98, mytitle, 'fontsize', 10, 'fontw', 'bold');
%                 end

                %% save
                figName = ['ROI_' num2str(aa) '_FB_' num2str(freq) '_ch' num2str(ch) '_' list_AA{aa} '_' groupInfo.list_FB{freq}];
                outDir = [params.storage.dir_results filesep params.storage.outName filesep 'ch2roi' groupInfo.dirSuffix filesep anatAtlas_name filesep 'xover_singleCh_v2' filesep list_AA{aa} filesep groupInfo.list_FB{freq}];
                fig_save(f, figName, outDir, 'format','png', 'res',600);
                close(f); 
                
                if ch == 50
                    why;
                end
            end
            % -------------------------------------------------------------
            
        end
        disp(['Switching points: aa = ' num2str(aa) ', freq = ' num2str(freq) ' - done!'])
    end
end
