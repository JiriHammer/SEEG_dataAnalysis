function chGroups = ch2roi_getSelChnls(params, subjTag)
% returns selected channels for subjTag
% based in criteria in params.connectivity.selectedChnls, where: 
%  - params.connectivity.selectedChnls{nn,1} = brainAtlas (e.g., 'Yeo7', 'Yeo17', 'Mars', ...)
%  - params.connectivity.selectedChnls{nn,2} = aaName (e.g., 'Default', 'Dorsal Attention', ...)
%  - params.connectivity.selectedChnls{nn,3} = significance of activation (e.g., 'was_sgnf', 'any_sgnf', 'non_sgnf')
%  - params.connectivity.selectedChnls{nn,4} = freq.bands of activation (e.g., 'beta', 'loGamma', 'hiGamma')
% algorithm:
%  - sets: groupInfo struct
%  - calls: ch2roi_wrapper -> GRUP struct
% output:
%  - selCh_all = indices to trialsData of activated channels
%  - chGroups = struct, where:

% (c) Jiri, Oct23

%% default settings
if ~isfield(params, 'P_level')
    params.P_level = 0.05;      % significance leel
end

if ~isfield(params, 't_siglen')
    params.t_siglen = 0.10;     % in [s], min. time of cont. significance
end

%% init
chGroups = struct;

%% get selected channels for each ROI (nn = neural network)
for nn = 1:size(params.connectivity.selectedChnls,1)
    
    %% vars to pass
    brainAtlas = params.connectivity.selectedChnls{nn,1};
    sel_AA = params.connectivity.selectedChnls{nn,2};
    sel_sgnf = params.connectivity.selectedChnls{nn,3};
    sel_FB = params.connectivity.selectedChnls{nn,4};
        
    %% design groupInfo & call ch2roi_wrapper -> GRUP
    params.storage.subjList = [];
    params.storage.subjList{1} = subjTag;    
    groupInfo = struct(...
        'brainAtlas', {{brainAtlas}}, ... % brain atlases: 'Yeo7','Yeo7_addLat','Yeo17','Yeo17_addLat','Mars','Mars_addLat'
        'freqBands', {sel_FB}, ...  % freq. bands (FB)
        'signifOfActiv', {{sel_sgnf}}, ...   % significance of activation
        'significanceLevel', params.P_level, ...      % significance level for ROI (e.g. DMN)
        'significanceTimeThr', params.t_siglen, ...    % in [s], significance time threshold (ch is signif. if its significance > threshold) 
        'GW_selection', {{'any'}}, ...      % grey / white matter: choices: 'gm','wm','csf','bone', 'any'
        'GW_threshold', 60, ...    % probability, in [%]
        'plotStyle', {{'allFB_allROI'}}, ... 'oneFB_allROI','allFB_allROI','spectra_SNR'
        'maxDistToROI', 10, ...   % in [mm], maximum distance (tolerance) to region of interest
        'minSubjInROI', 1, ...  % minimum of subjects (patients) in a ROI (if a ROI has less patients, it is not considered in plots)
        'subtractAvgResponse', false, ...    % subtracts mean response over classes (conditions) from each clz
        'dirSuffix', '_DTF_singleSubj', ... % adds suffix to 'ch2roi' directory name -> ch2roi_outDirSuffix
        'figSuffix', ['_' subjTag], ... % adds suffix to figure name
        'time2plot', [-2.3,2.8] ... % in [s], w.r.t. cutting point (= 0)
        );
    clear GRUP;
    ch2roi_wrapper;
    
    %% cache file
    cacheFile = [params.storage.dir_results filesep params.storage.outName filesep subjTag filesep 'cacheFile.mat'];
    assert(exist(cacheFile,'file') == 2);

    %% load structure H & selected channels -> H, selCh_H
    clear H selCh_H_resp;
    load(cacheFile, 'H', 'selCh_H_resp');   
    selCh_H = selCh_H_resp; 
    if isfield(params, 'nCh')
        assert(params.nCh == size(selCh_H,2));
    end
    
    %% channel names from H -> H_chNames
    H_chNames = cell(1,size(selCh_H,2));
    for ch = 1:size(selCh_H,2)
        H_chNames{ch} = H.channels(selCh_H(ch)).name;
    end
    
    %% find sel_AA in groupInfo.list_AA -> aa
    [tf, aa] = ismember(sel_AA, groupInfo.list_AA(:,1));
    if ~tf
        error(['Anatomic area = ' sel_AA ' was not found in ' brainAtlas]);
    end
    
    %% union across selected freq. bands (? TO DO: intersect ?)
    ch_name = [];   % channel names
    ch_inds = [];   % indices into trials.data
    ch_inds_H = []; % indices into H
    ch_MNIs = [];
    for fb = 1:size(sel_FB,2)
        [tf, freq] = ismember(sel_FB{fb}, groupInfo.list_FB(:,1));
        assert(tf);
        for ch = 1:size(GRUP.chSubj_inds{aa,freq,1},1)
            i_ch_trials = GRUP.chSubj_inds{aa,freq,1}(ch,1);    % ch index to trials, {aa,freq,1} => significant channel
            sName = GRUP.chSubj_name{aa,freq,1}{ch,2};
            cName = GRUP.chSubj_name{aa,freq,1}{ch,1};
            assert(strcmp(subjTag, sName));
            assert(ismember(cName, H_chNames));
            
            % TBD: older version without checking if a contact from bipolar pair is in the other group
%             if ~ismember(i_ch_trials, ch_inds)
%                 if (nn==1) || (nn == 2 && ~ismember(i_ch_trials, chGroups.ch_inds_trials{1,1})) % in case a channel was assigned to the 1st group already, then do not include!!!, works only for 2 NN !!!
%                     ch_inds = cat(2, ch_inds, i_ch_trials);
%                     ch_name = cat(2, ch_name, {cName});
%                     ch_inds_H = cat(2, ch_inds_H, selCh_H(i_ch_trials));
%                 end
%             end
            if ch2roi_getSelCh_add(i_ch_trials, ch_inds, chGroups, nn, cName)
                ch_inds = cat(2, ch_inds, i_ch_trials);
                ch_name = cat(2, ch_name, {cName});
                ch_inds_H = cat(2, ch_inds_H, selCh_H(i_ch_trials));   
                ch_MNIs = cat(1, ch_MNIs, [H.channels(selCh_H(i_ch_trials)).MNI_x, H.channels(selCh_H(i_ch_trials)).MNI_y, H.channels(selCh_H(i_ch_trials)).MNI_z]);
            end
        end
    end
    chGroups.ch_inds_trials{nn,1} = ch_inds;
    chGroups.ch_names{nn,1} = ch_name;
    chGroups.nn_name{nn,1} = sel_AA;
    chGroups.brainAtlas{nn,1} = brainAtlas;
    chGroups.ch_inds_H{nn,1} = ch_inds_H;
    chGroups.ch_MNIs{nn,1} = ch_MNIs;
end

chGroups.H = H;

