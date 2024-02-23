%% ch2roi_wrapper: assignment of iEEG channel to ROI (region of interest)
% assignment based on:
%   - brain atlas: isarg_atlas from RJ (brain_atlas_v07_standalone.m)
%   - frequency band: in struct trialsData_freqBand (e.g. trialsData_hiGamma)
%   - significance of activation: was*, not*, any*
% coded in structure 'groupInfo', full example:
% groupInfo = struct(...
%     'brainAtlas', {{'Yeo7', 'Yeo7_addLat', 'Yeo17', 'Yeo17_addLat', 'Destrie', 'DestLat', 'Mars', 'Mars_addLat', 'LabelsAK', 'LabelsAK_addDir'}}, ... % brain atlases
%     'freqBands', {{'erp','delta','theta','alpha','beta','loGamma','hiGamma'}}, ...  % freq. bands
%     'signifOfActiv', {{'was_sgnf','not_sgnf','any_sgnf'}}, ...   % significance of activation
%     'GW_selection', {{'any','gm','wm'}}, ...      % grey / white matter: choices: 'gm','wm','csf','bone', 'any'
%     'GW_threshold', 60, ...    % probability, in [%]
%     'plotStyle', {{'tResFB_oneROI_3D','oneFB_allROI','allFB_allROI','allFB_oneROI_3D','stacked_FB','spectra_SNR'}}, ...
%     'maxDistToROI', 2, ...   % in [mm], maximum distance (tolerance) to region of interest
%     'minSubjInROI', 3, ...   % minimum of subjects (patients) in a ROI (if a ROI has less patients, it is not considered in plots)    
%     'dirSuffix', '', ... % adds suffix to 'ch2roi' directory name -> ch2roi_outDirSuffix
%     'figSuffix', '', ... % adds suffix to figure name
%     'time2plot', [-2,2] ...
%     );
% required variables: groupInfo, params

% (c) Jiri, Jul22

%% required / default fields
groupInfo.path_isarg_atlas = params.path2others.isarg_atlas;
if isempty(params.storage.subjList)
    disp('WARNING: no subjects selected.');
    return;
end
if ~isfield(groupInfo, 'GW_selection')
    groupInfo.GW_selection = {'any'};
end
if ~isfield(groupInfo, 'dirSuffix')
    groupInfo.dirSuffix = '';
end
if ~isfield(groupInfo, 'figSuffix')
    groupInfo.figSuffix = '';
end
groupInfo.clzNames = params.triggering.classes;

%% go through selected brain atlases
for atls = 1:size(groupInfo.brainAtlas,2)
    anatAtlas_name = groupInfo.brainAtlas{atls};
    list_AA = anatomicalAreas_getList(anatAtlas_name, groupInfo);   % define ROIs (anatom. areas) 
    disp(['Grouping analysis, atlas = : ' anatAtlas_name ' ...']);
    
    %% go trough different significance of activations (SOA)
    for sgnf = 1:size(groupInfo.signifOfActiv,2)
        sel_SOA = groupInfo.signifOfActiv{sgnf};
        
        for gw = 1:size(groupInfo.GW_selection,2)
            sel_GW = groupInfo.GW_selection{gw};        
            
            list_anatomy_signif = {anatAtlas_name, nan, sel_SOA, sel_GW}; % 1 x 4 cell: anatAtlas, freqBand, signifType, greyWhiteMatter; for example: 'Yeo7_networks', 'all', 'was_sgnf'; ...   % 'not_sgnf', 'all'

            %% >>> load & assign ch2ROI -> GRUP (grouped data) <<<
            groupInfo.list_AA = list_AA;
            groupInfo.list_FB = groupInfo.freqBands';
            groupInfo.list_anatomy_signif = list_anatomy_signif;
            GRUP = ch2roi_load(params, groupInfo);

            %% number of used ROIs
            if ~isfield(groupInfo, 'minSubjInROI')
                groupInfo.minSubj = 3;
            else
                groupInfo.minSubj = groupInfo.minSubjInROI;
            end
            [groupInfo.nROI, groupInfo.nSpectra]  = getUsedROIs(GRUP, groupInfo.minSubj);

            %% plot: AUC for all FB & for all ROIs
            if ismember('AUC_allFB_allROI', groupInfo.plotStyle)
                groupInfo.outDir = [params.storage.dir_results filesep params.storage.outName filesep 'ch2roi' groupInfo.dirSuffix filesep anatAtlas_name filesep 'AUC_allFB_allROI'];
                groupInfo.figName = [anatAtlas_name '_' sel_SOA '_' sel_GW groupInfo.figSuffix];
                ch2roi_plot_AUC_allFB_allROI(GRUP, groupInfo);
            end
            
            %% plot: sel freq. bands
            if ismember('oneFB_allROI', groupInfo.plotStyle)
                for freq = 1:size(groupInfo.list_FB,1)
                    groupInfo.freqBand = groupInfo.freqBands{freq};
                    [tf, groupInfo.i_freq] = ismember(groupInfo.freqBand, groupInfo.list_FB);
                    assert(tf);
                    groupInfo.outDir = [params.storage.dir_results filesep params.storage.outName filesep 'ch2roi' groupInfo.dirSuffix filesep anatAtlas_name filesep 'selFreq'];
                    groupInfo.figName = [anatAtlas_name '_' num2str(freq) '_' groupInfo.freqBand '_' sel_SOA '_' sel_GW groupInfo.figSuffix];
                    ch2roi_plot_selFreq(GRUP, groupInfo);
                end
            end

            %% plot: all freq. bands for all ROIs
            if ismember('allFB_allROI', groupInfo.plotStyle)
                groupInfo.outDir = [params.storage.dir_results filesep params.storage.outName filesep 'ch2roi' groupInfo.dirSuffix filesep anatAtlas_name filesep 'allFreq'];
                groupInfo.figName = [anatAtlas_name '_' sel_SOA '_' sel_GW groupInfo.figSuffix];
                ch2roi_plot_allFreq(GRUP, groupInfo);
            end

            %% plot mean over subj: all freq. bands for all ROIs
            if ismember('allFB_allROI_subjMean', groupInfo.plotStyle)
                groupInfo.outDir = [params.storage.dir_results filesep params.storage.outName filesep 'ch2roi' groupInfo.dirSuffix filesep anatAtlas_name filesep 'allFreq_subjMean'];
                groupInfo.figName = [anatAtlas_name '_' sel_SOA '_' sel_GW groupInfo.figSuffix];
                ch2roi_plot_allFreq_subjMean(GRUP, groupInfo);
            end
            
            %% plot: all freq. bands for all ROIs
            if ismember('allFB_allROI_SNR', groupInfo.plotStyle)
                groupInfo.outDir = [params.storage.dir_results filesep params.storage.outName filesep 'ch2roi' groupInfo.dirSuffix filesep anatAtlas_name filesep 'SNR_allFreq'];
                groupInfo.figName = ['SNR_' anatAtlas_name '_' sel_SOA '_' sel_GW groupInfo.figSuffix];
                ch2roi_plot_allFreq_SNR(GRUP, groupInfo);
            end
            
            %% plot: time resolved brain activations of selected FB: brain3D + all freq. bands for each ROI
            if ismember('tResFB_oneROI_3D', groupInfo.plotStyle)
                for roi = 1:size(GRUP.trials,1)
                    for freq  =1:size(GRUP.trials,2)
                        groupInfo.i_roi = roi;      % pointer to selected ROI
                        groupInfo.i_freq = freq;      % pointer to selected ROI
                        groupInfo.outDir = [params.storage.dir_results filesep params.storage.outName filesep 'ch2roi' groupInfo.dirSuffix filesep ...
                                            anatAtlas_name filesep 'tSeq_brain3D' filesep list_AA{roi,1} '_' sel_SOA '_' num2str(freq) '_' groupInfo.list_FB{freq}];
                        groupInfo.figName = [list_AA{roi,1} '_' sel_SOA '_' sel_GW groupInfo.figSuffix];
                        ch2roi_plot_tSeqFB_brain3D(GRUP, groupInfo);
                    end
                end
            end

            %% plot: brain3D + all freq. bands for each ROI
            if ismember('allFB_oneROI_3D', groupInfo.plotStyle)
                for roi = 1:size(GRUP.trials,1)
                    groupInfo.i_roi = roi;      % pointer to selected ROI
                    groupInfo.outDir = [params.storage.dir_results filesep params.storage.outName filesep 'ch2roi' groupInfo.dirSuffix filesep anatAtlas_name filesep 'allFB_brain3D'];
                    groupInfo.figName = [anatAtlas_name '_' list_AA{roi,1} '_' sel_SOA '_' sel_GW groupInfo.figSuffix];
                    ch2roi_plot_allFB_brain3D(GRUP, groupInfo);
                end
            end

            %% plot: stacked freq. bands (stacked FB per clz, SNR, *lags, nCh, nSubj)
            if ismember('stacked_FB', groupInfo.plotStyle)
                groupInfo.outDir = [params.storage.dir_results filesep params.storage.outName filesep 'ch2roi' groupInfo.dirSuffix filesep anatAtlas_name filesep 'stackedFreqBands'];
                groupInfo.figName = [anatAtlas_name '_' sel_SOA '_' sel_GW groupInfo.figSuffix];
                ch2roi_plot_stackFreq(GRUP, groupInfo);
            end

            %% plot: spectra
            if ismember('spectra_SNR', groupInfo.plotStyle)
                groupInfo.outDir = [params.storage.dir_results filesep params.storage.outName filesep 'ch2roi' groupInfo.dirSuffix filesep anatAtlas_name filesep 'spectra_SNR'];
                groupInfo.figName = [anatAtlas_name '_' sel_SOA '_' sel_GW groupInfo.figSuffix];
                if groupInfo.nSpectra > 20
                    ch2roi_plot_spectra(GRUP, groupInfo);
                else
                    ch2roi_plot_spectra_oneFig(GRUP, groupInfo);
                end
            end
            
        end     % of GM/WM
    end     % of soa (signif.)
end     % of atlas

        










%% old stuff (recycle? delete?) (now = Jul22)
% ----->>> select anatomical areas grouping  <<<-------
% for implementation, see: anatomicalAreas_getList.m
% excel tables in: ..\neuroHammer\plot_3Dbrain
% AA_type = 'neurologists_anatomy';                 % -> BrainAtlas_labels.xlsx
% AA_type = 'neurologists_laterality';              % -> BrainAtlas_labels.xlsx
% AA_type = 'neurologists_grouping';                  % -> BrainAtlas_clusters_fromAdam.xlsx
% AA_type = 'neurologists_grouping_laterality';     % -> BrainAtlas_clusters_fromAdam.xlsx

% BP: F. Kanout
% Default.nii = default mode network (DMN)
% Frontoparietal.nii = central executive network (CEN)
