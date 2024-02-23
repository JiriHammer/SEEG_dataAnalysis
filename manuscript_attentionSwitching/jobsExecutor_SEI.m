%% trigger analysis: SWITCH-EI
% class 1: switching: E -> I
% class 2: switching: I -> E
% (c) Jiri, May20

%% specify analysis to be performed -> runner
runner = struct;
runner.newJobs = false; % set to 'false' if the results are already computed
runner.newFDR = false;  % true = updates the FDR correction over all patients (and freq. bands, ...) 

%% job settings: common (data preprocessing & analysis) -> params
% params = SEI_job_commonSettings();    % TO DO: del fcn
params = SEI_job_settingsCommon();

%% job settings: specific (define spectral method) -> params
spectralMethod = 'v31_stft_session_bip';           % ~ v31 (BIP), similar to v12 (STFT base = sess, no zscore, no rej, FB = loFreq0)
% spectralMethod = 'v30_stft_baseRS_bip';           % ~ v30 (BIP), similar to v24 (new= no zscore, no rej, FB = loFreq0), STFT with RS baseline
% spectralMethod = 'v30_stft_baseRS_car';           % ~ v30 (CAR), similar to v24 (new= no zscore, no rej, FB = loFreq0), STFT with RS baseline
% spectralMethod = 'v29_stft_baseRS_bip';           % ~ v29 (BIP), similar to v27 (new=zscore, rej = 10 SD, FB = loFreq0), STFT with RS baseline
% spectralMethod = 'v28_stft_baseRS_bip';         % ~ v28 (BIP), new = bip: 12-34, STFT with RS baseline
% spectralMethod = 'v27_stft_baseRS_bip';         % ~ v27 (BIP), new=zscore, STFT with RS baseline
% spectralMethod = 'v26_rawRS_bip';         % ~ v26 (BIP), just load & preprocess RS, for RS connectivity
% spectralMethod = 'v25_stft_baseRS_diff_car';         % ~ v24 (CAR), STFT with RS baseline
% spectralMethod = 'v25_stft_baseRS_diff_bip';         % ~ v24 (BIP), STFT with RS baseline
% spectralMethod = 'v24_stft_baseRS_car';         % ~ v24 (CAR), STFT with RS baseline
% spectralMethod = 'v24_stft_baseRS_bip';         % ~ v24 (BIP), STFT with RS baseline
% spectralMethod = 'bipDiffZ_trials4EC';    % ~ v23
% spectralMethod = 'diffZ_trials4GC';    % ~ v21 & v22
% spectralMethod = 'surrog_trials4GC';      % ~ v20
% spectralMethod = 'alphaBPHT_trials4GC'; % ~ v19
% spectralMethod = 'export_trials4FC';    % TBD: 'export_trials4GC';    % ~ v18
% spectralMethod = 'dtf_trials';          % ~ v17
% spectralMethod = 'coh_trials';          % ~ v16
% spectralMethod = 'cwt_trials';          % ~ v15
% spectralMethod = 'psd_trials';          % ~ v14 (STFT implementation by R. Janca)
% spectralMethod = 'stft_lowPassDiff';    % ~ v13
% spectralMethod = 'stft_sessions';       % ~ v12
% spectralMethod = 'stft_trials';
% spectralMethod = 'stft_baseRS';         % ~ v11 (CAR) - manually change to 'car' in SEI_job_settingsCommon
% spectralMethod = 'stft_baseRS_bip';     % ~ v11b (BIP) - manually change to 'bip' in SEI_job_settingsCommon
% spectralMethod = 'mtft';                % ~ v6
params = SEI_job_spectralMethod(params, spectralMethod);

%% job settings: define the classes for the triggers -> params
clzDefName = 'switchin_EI_IE';
% clzDefName = 'switchin_IE_EE';
% clzDefName = 'switchin_EE_II';
% clzDefName = 'switchin_EI_IE_EE_II';
% clzDefName = 'task_E_I';
params = SEI_job_classDefinitions(params, clzDefName);

%% logging to log_commandWindow
setLogFile_commandWindow(params);

%% >>>>>>>>>>>>>>>>>>>>>>>> run jobs <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
if runner.newJobs
    for subj = 1:size(params.storage.subjList,1)
        main_execute(params, subj);
%         main_execute;
    end
    disp('Done.');
    disp(['Time: ' datestr(now, 'dd-mmm-yyyy HH:MM:SS')]);
end
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

%% change paths to output directory of results (uncomment)
% dir_outputResults = 'G:\dox\proj_switching_EI';
% params = changeDir_outputResults(params, dir_outputResults);

%% FDR correction
% params.storage.subjList = subjList_updateAnalyzed(params);
% if runner.newFDR
%     significanceLevel = 0.05;
%     fields_signifVals = {'pVals','hVals'};
%     poolSubjects_fdrCorrection(params, significanceLevel, fields_signifVals);
% 
%     fields_signifVals = {'pVals_base','hVals_base'};
%     poolSubjects_fdrCorrection(params, significanceLevel, fields_signifVals);
% end

%% ---------------- PLOT channels: slices & activations -------------------
% sort channels into folders based on groupInfo.anatomy_signif_list
sliceInfo = struct;
sliceInfo.time2plot = [-2.1,2.1]; ... [-3.7,3.7] ...[-1,4] ... % in [s], w.r.t. cutting point (= 0), currently only used for significance, not for plotting!!!
sliceInfo.significanceTimeThr = 0.25;   % in [s], significance time threshold (ch is signif. if its significance > threshold) 
sliceInfo.subtractAvgResponse = true;    % subtracts mean response over classes (conditions) from each clz
sliceInfo.maxDistToROI = 10;     % in [mm], maximum distance (tolerance) to region of interest
sliceInfo.file_brainClusters = 'C:\Users\Jiri\Documents\MATLAB\brainClusters_tables\LabelsAK_1.xlsx';
sliceInfo.anatomy_signif_list = {...    % anat_atlas, freqBand, signifType
%     'Yeo7_addLat', 'alpha', 'was_sgnf'; ...   % 'not_sgnf', 'all'
%     'Yeo7_addLat', 'beta', 'was_sgnf'; ...   % 'not_sgnf', 'all'
    'Yeo7', 'hiGamma', 'was_sgnf'; ...   % 'not_sgnf', 'all'
    'Yeo7', 'hiGamma', 'not_sgnf'; ...   % 'not_sgnf', 'all'
    'Yeo17', 'hiGamma', 'was_sgnf'; ...   % 'not_sgnf', 'all'
    'Yeo17', 'hiGamma', 'not_sgnf'; ...   % 'not_sgnf', 'all'  
    'Yeo7', 'loGamma', 'was_sgnf'; ...   % 'not_sgnf', 'all'
    'Yeo7', 'loGamma', 'not_sgnf'; ...   % 'not_sgnf', 'all'
    'Yeo17', 'loGamma', 'was_sgnf'; ...   % 'not_sgnf', 'all'
    'Yeo17', 'loGamma', 'not_sgnf'; ...   % 'not_sgnf', 'all'      
    'Yeo7', 'loFreq0', 'was_sgnf'; ...   % 'not_sgnf', 'all'
    'Yeo7', 'loFreq0', 'not_sgnf'; ...   % 'not_sgnf', 'all'
    'Yeo17', 'loFreq0', 'was_sgnf'; ...   % 'not_sgnf', 'all'
    'Yeo17', 'loFreq0', 'not_sgnf'; ...   % 'not_sgnf', 'all'  
    };

params.storage.subjList = {
%     '20_PR11'; ...
    '20_PR14'; ...
%     '22_PR5'; ...
%     '23_PR2'; ...
    };
params.plot_brainSlices = struct(...           % user interface structure: holds most (but not all!) of the user settings
    'plot_slices', true, ...                            % 2D brain slices
    'volumes2plot', {{...       % volumes to plot (coregistered .nii files)
        'wT1', gray(128), [1, 1];  ...     % name & colormap & transparency limits
%         'wrCT', copper(128), [0 0.8]; ...
%         'ch-spec_Yeo7', [0 1 0; 0 1 0], [0, 0.3]; ...  % prev. used colormap: autumn(128)
        'ch-spec_Yeo7', [3 213 255; 3 213 255]./255, [0, 0.3]; ...  % prev. used colormap: autumn(128)
%         'ch-spec_Yeo17', [1 0 0; 1 0 0], [0, 0.3]; ...
        }}, ...             
    'slicePlanes', {{'axial','sagittal','coronal'}},... % planes to plot, options: 'axial','sagittal','coronal'
    'size_interpolate', 0.7, ...1.0, ...                        % in [mm], voxel size to which the brain is interpolated
    'printFormats', {{'png'}}, ...
    'printResolution', 0 ...                           % choices: 0 (= screen resolution) or 600 (= dpi). Resolution of the figures.                          
); 
poolSubjects_slicesActivity(params, sliceInfo);

%% export contacts into table (.xlsx)
pool_chnls2table(params, sliceInfo);

%% pooling (grouping) analysis: all
% groupInfo = struct(...
%     'brainAtlas', {{'Yeo17', 'Yeo7', 'Yeo7_addLat', 'LabelsAK', 'LabelsAK_addLat', 'Yeo17_addLat', 'Destrie'}}, ... % brain atlases
%     'freqBands', {{'erp','delta','theta','alpha','beta','loGamma','hiGamma'}}, ...  % freq. bands (FB)
%     'signifOfActiv', {{'was_sgnf','not_sgnf','any_sgnf'}}, ...   % significance of activation
%     'plotStyle', {{'oneFB_allROI','allFB_allROI','allFB_oneROI_3D','stacked_FB','spectra_SNR'}}, ...
%     'maxDistToROI', 5, ...   % in [mm], maximum distance (tolerance) to region of interest
%     'minSubjInROI', 3, ...  % minimum of subjects (patients) in a ROI (if a ROI has less patients, it is not considered in plots)
%     'dirSuffix', '', ... % adds suffix to 'ch2roi' directory name -> ch2roi_outDirSuffix
%     'figSuffix', '', ... % adds suffix to figure name
%     'time2plot', [-2,2] ...
%     );
% groupInfo = struct(...
%     'brainAtlas', {{'Yeo7','Yeo17'}}, ... % brain atlases: 'Yeo7','Yeo7_addLat','Yeo17','Yeo17_addLat','Mars','Mars_addLat'
%     'freqBands', {{'theta','alpha','beta','loGamma','hiGamma'}}, ...  % freq. bands (FB)
%     'signifOfActiv', {{'was_sgnf','not_sgnf','any_sgnf'}}, ...   % significance of activation
%     'significanceLevel', 0.05, ...      % significance level for ROI (e.g. DMN)
%     'GW_selection', {{'any'}}, ...      % grey / white matter: choices: 'gm','wm','csf','bone', 'any'
%     'GW_threshold', 60, ...    % probability, in [%]
%     'plotStyle', {{'oneFB_allROI','allFB_allROI','stacked_FB','spectra_SNR','AUC_allFB_allROI'}}, ...
%     'maxDistToROI', 10, ...   % in [mm], maximum distance (tolerance) to region of interest
%     'minSubjInROI', 1, ...  % minimum of subjects (patients) in a ROI (if a ROI has less patients, it is not considered in plots)
%     'dirSuffix', '', ... % adds suffix to 'ch2roi' directory name -> ch2roi_outDirSuffix
%     'figSuffix', '', ... % adds suffix to figure name
%     'time2plot', [-3.7,3.7] ...[-1,4] ... 
%     );
% groupInfo.AUC_settings = struct(...
%     'baselineWindow', [-0.5 0], ...
%     'timeWindow', [0 0.5]  ...
%     );
% groupInfo = struct(...
%     'brainAtlas', {{'Mars_addLat'}}, ... % brain atlases
%     'freqBands', {{'hiGamma','beta','alpha'}}, ...  % freq. bands (FB)
%     'signifOfActiv', {{'was_sgnf','not_sgnf','any_sgnf'}}, ...   % significance of activation
%     'plotStyle', {{'oneFB_allROI'}}, ...
%     'maxDistToROI', 10, ...   % in [mm], maximum distance (tolerance) to region of interest
%     'minSubjInROI', 3, ...  % minimum of subjects (patients) in a ROI (if a ROI has less patients, it is not considered in plots)
%     'dirSuffix', '', ... % adds suffix to 'ch2roi' directory name -> ch2roi_outDirSuffix
%     'figSuffix', '', ... % adds suffix to figure name
%     'time2plot', [-3.7,3.7] ...
%     );
groupInfo = struct(...
    'brainAtlas', {{'Yeo7','Yeo17'}}, ... % brain atlases: 'Yeo7','Yeo7_addLat','Yeo17','Yeo17_addLat','Mars','Mars_addLat'
    'freqBands', {{'delta','theta','alpha','beta','loGamma','hiGamma'}}, ...  % freq. bands (FB)
    'signifOfActiv', {{'was_sgnf','not_sgnf','any_sgnf'}}, ...   % significance of activation
    'significanceLevel', 0.05, ...      % significance level for ROI (e.g. DMN)
    'significanceTimeThr', 0.10, ...    % in [s], significance time threshold (ch is signif. if its significance > threshold) 
    'GW_selection', {{'any'}}, ...      % grey / white matter: choices: 'gm','wm','csf','bone', 'any'
    'GW_threshold', 60, ...    % probability, in [%]
    'plotStyle', {{'allFB_allROI_subjMean','oneFB_allROI','allFB_allROI','spectra_SNR','allFB_allROI_SNR'}}, ... 'oneFB_allROI','allFB_allROI','spectra_SNR'
    'maxDistToROI', 10, ...   % in [mm], maximum distance (tolerance) to region of interest
    'minSubjInROI', 1, ...  % minimum of subjects (patients) in a ROI (if a ROI has less patients, it is not considered in plots)
    'subtractAvgResponse', true, ...    % subtracts mean response over classes (conditions) from each clz
    'dirSuffix', '_p0.05_tSigThr0.10_wiSubtr', ... % adds suffix to 'ch2roi' directory name -> ch2roi_outDirSuffix
    'figSuffix', '', ... % adds suffix to figure name
    'time2plot', [-2.1,2.1] ... [-3.7,3.7] ...[-1,4] ... % in [s], w.r.t. cutting point (= 0)
    );
groupInfo.file_brainClusters = 'C:\Users\Jiri\Documents\MATLAB\brainClusters_tables\LabelsAK_1.xlsx';
ch2roi_wrapper; % script requires groupInfo struct !!!

%% pooling (grouping) analysis: single-subj test
sList = params.storage.subjList;
for subj = 1:size(sList,1)      % start from 20_PR7
    params.storage.subjList = [];
    params.storage.subjList{1} = sList{subj};
    
%     groupInfo = struct(...
%         'brainAtlas', {{'Yeo7','Yeo17'}}, ... % brain atlases: 'Yeo7','Yeo7_addLat','Yeo17','Yeo17_addLat','Mars','Mars_addLat'
%         'freqBands', {{'delta','theta','alpha','beta','loGamma','hiGamma'}}, ...  % freq. bands (FB)
%         'signifOfActiv', {{'was_sgnf','not_sgnf','any_sgnf'}}, ...   % significance of activation
%         'significanceLevel', 0.05, ...      % significance level for ROI (e.g. DMN)
%         'GW_selection', {{'any'}}, ...      % grey / white matter: choices: 'gm','wm','csf','bone', 'any'
%         'GW_threshold', 60, ...    % probability, in [%]
%         'plotStyle', {{'oneFB_allROI','allFB_allROI'}}, ...
%         'maxDistToROI', 10, ...   % in [mm], maximum distance (tolerance) to region of interest
%         'minSubjInROI', 1, ...  % minimum of subjects (patients) in a ROI (if a ROI has less patients, it is not considered in plots)
%         'dirSuffix', '_singleSubj', ... % adds suffix to 'ch2roi' directory name -> ch2roi_outDirSuffix
%         'figSuffix', ['_' sList{subj}], ... % adds suffix to 'ch2roi' directory name -> ch2roi_outDirSuffix
%         'time2plot', [-1.7,1.7] ...[-1,4] ... 
%         );
    groupInfo = struct(...
        'brainAtlas', {{'Yeo7','Yeo17'}}, ... % brain atlases: 'Yeo7','Yeo7_addLat','Yeo17','Yeo17_addLat','Mars','Mars_addLat'
        'freqBands', {{'delta','theta','alpha','beta','loGamma','hiGamma'}}, ...  % freq. bands (FB)
        'signifOfActiv', {{'was_sgnf','not_sgnf','any_sgnf'}}, ...   % significance of activation
        'significanceLevel', 0.05, ...      % significance level for ROI (e.g. DMN)
        'significanceTimeThr', 0.10, ...    % in [s], significance time threshold (ch is signif. if its significance > threshold) 
        'GW_selection', {{'any'}}, ...      % grey / white matter: choices: 'gm','wm','csf','bone', 'any'
        'GW_threshold', 60, ...    % probability, in [%]
        'plotStyle', {{'oneFB_allROI','allFB_allROI','spectra_SNR','allFB_allROI_SNR'}}, ... 'oneFB_allROI','allFB_allROI','spectra_SNR'
        'maxDistToROI', 10, ...   % in [mm], maximum distance (tolerance) to region of interest
        'minSubjInROI', 1, ...  % minimum of subjects (patients) in a ROI (if a ROI has less patients, it is not considered in plots)
        'subtractAvgResponse', true, ...    % subtracts mean response over classes (conditions) from each clz
        'dirSuffix', ['_p0.05_tSigThr0.10_wiSubtr' filesep 'singleSubj'], ... % adds suffix to 'ch2roi' directory name -> ch2roi_outDirSuffix
        'figSuffix', ['_' sList{subj}], ... % adds suffix to figure name
        'time2plot', [-2.1,2.1] ... [-3.7,3.7] ...[-1,4] ... % in [s], w.r.t. cutting point (= 0)
        );    
    ch2roi_wrapper; % script requires groupInfo struct !!!
end

%% transparent 3D brains: DMN-DAN, vals + signif + SNR
groupInfo = struct;
groupInfo.time2plot = [-0.6:0.2:0.5];
groupInfo.vals2plot = {...
    'SNR_trials', 'prctile5', jet(256), 'circleSize_linScaled';
    'trials', 'prctile5', brewermap(256,'*RdBu'), 'circleSize_absScaled';
    'chnls_hVals', [0, 1], [0 0 0; 1 0 0], 'default';
    };
groupInfo.anatomy_signif_list = {...    % anat_atlas, freqBand, signifType
    'Yeo7_networks', 'beta', 'all', {'Default'; 'Dorsal Attention'}; ...   % 'not_sgnf', 'all'
    'Yeo7_networks', 'hiGamma', 'all', {'Default'; 'Dorsal Attention'};
    };
params.plot_brain3D = struct(...        % user interface structure: holds most (but not all!) of the user settings, see also inside brain3D_plot.m
    'plot_brain3D', true, ...           % 3D brain model
    'volumes2plot', {{...               % volumes to plot (coregistered .nii files)
        'icbm152', [1 1 1], 0.07;  ...     % 1x3 cell: name, color, transparency (other choices for norm brain = wc1T1_colin27)
        'wDefault', [0 1 0], 0.2; ...
        'wDorsal Attention', [165 113 78]./255, 0.2; ...
        }}, ...          
    'chVals_asCircles', true, ...       % plots channel values as color circles
    'chVals_asPatches', false, ...       % each channel value is mapped onto surface (but may be overwritten by latter points!), see also distanceThreshold for focality inside brain3D_plot.m
    'chVals_colorMap', [0 0 0; 1 0 0], ...    % colormap of channel values
    'chVals_patchTransp', 0.2, ...      % transparency of channel values
    'circleSizeLims', [5 60], ...       % in points, maps chnl values to circle size
    'size_interpolate', 1.0, ...        % in [mm], voxel size to which the brain is interpolated
    'backgroundColor', 'w', ...         % figure / axes background color, good choices: 'w', 'k'
    'text_LR', true, ...                % marks L / R side
    'text_tag', '', ...    % string
    'viewAngle', [0, 90], ...           % axial = 0,90; sagittal = 100,10; coronal front = 0,180; coronal back = 0,0
    'visible_axis', false ...           % if axis & their labels are visible
); 
% poolSubjects_brain3D(params, groupInfo);
poolSubjects_brain3D;

%% transparent 3D brains: vals of all chnls
groupInfo = struct;
groupInfo.time2plot = [-2.5:0.5:-1.0, -0.8:0.1:0.8, 1:0.5:2.5];
groupInfo.vals2plot = {...
    'SNR_trials', 'prctile5', jet(256), 'circleSize_linScaled';
...    'trials', 'prctile5', brewermap(256,'*RdBu'), 'circleSize_absScaled';
    };
groupInfo.anatomy_signif_list = {...    % anat_atlas, freqBand, signifType
    'Yeo7_networks', 'beta', 'all', 'Yeo7_networks'; ...   % 'not_sgnf', 'all'
    'Yeo7_networks', 'hiGamma', 'all', 'Yeo7_networks';
    };
params.plot_brain3D = struct(...        % user interface structure: holds most (but not all!) of the user settings, see also inside brain3D_plot.m
    'plot_brain3D', true, ...           % 3D brain model
    'volumes2plot', {{...               % volumes to plot (coregistered .nii files)
        'wc1T1_colin27', [1 1 1], 0.07;  ...     % 1x3 cell: name, color, transparency
        'wDefault', [0 1 0], 0.2; ...
        'wDorsal Attention', [165 113 78]./255, 0.2; ...
        }}, ...          
    'chVals_asCircles', true, ...       % plots channel values as color circles
    'chVals_asPatches', false, ...       % each channel value is mapped onto surface (but may be overwritten by latter points!), see also distanceThreshold for focality inside brain3D_plot.m
    'chVals_colorMap', brewermap(256,'*RdBu'), ...    % colormap of channel values
    'chVals_patchTransp', 0.2, ...      % transparency of channel values
    'circleSizeLims', [5 40], ...       % in points, maps chnl values to circle size
    'size_interpolate', 1.0, ...        % in [mm], voxel size to which the brain is interpolated
    'backgroundColor', 'w', ...         % figure / axes background color, good choices: 'w', 'k'
    'text_LR', true, ...                % marks L / R side
    'text_tag', '', ...    % string
    'viewAngle', [0, 90], ...           % axial = 0,90; sagittal = 100,10; coronal front = 0,180; coronal back = 0,0
    'visible_axis', false ...           % if axis & their labels are visible
); 
% poolSubjects_brain3D_activity(params, groupInfo);
poolSubjects_brain3D_activity;

%% compare across analysis types (in outNames)
pathBeg = 'G:\dox\ms_switch_EI\data';
dirsList = {...
    [pathBeg filesep 'v2_stft_sessions\switchin_EE_II_car']; ...
    [pathBeg filesep 'v2_stft_sessions\switchin_EI_IE_car']; ...
    [pathBeg filesep 'v2_stft_sessions\switchin_EI_IE_EE_II_car']; ...
    [pathBeg filesep 'v11_stft_baseRS\switchin_EI_IE_car']; ...
    [pathBeg filesep 'v12_stft_sessions_zscore\switchin_EI_IE_car']; ...
    };
figDir = [pathBeg filesep 'chnlsActivations_dirsList'];
pooling_outNames(params, dirsList, figDir);

%% ### close logging of command window input (if program crashed, do it manually!) ###
if strcmp(get(0,'Diary'),'on')
    diary off;
end
