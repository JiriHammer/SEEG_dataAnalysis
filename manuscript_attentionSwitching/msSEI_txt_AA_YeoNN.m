%% ms SEI table - anatomy area (AA) in brainAtlas (Yeo NN: Yeo7 or Yeo17)
% lists how many channels from each Yeo NN fall into an AA
% AA can be:
%   - neurologist label
%   - Mars atlas, ...

% (c) Jiri, Oct23

%% settings: v1 - relaxed conditions (p-level, length of significance, ...)
dirSuffix = '_msFig_AA_v1_relaxCond';
t_sel = [-2.3, 2.3];    % the true sel time will be cut by 1 s on edges (= template length): [-2, 2] -> [-1, 1]
P_level = 0.05;         % significance level
t_siglen = 0.10;        % in [s], significance time threshold (ch is signif. if its significance > threshold) 
brainAtlas = 'Yeo7';

%% settings: v3
% figureSuffix = '_v3_strictCond';
% t_sel = [-1.0, 2.0];    % the true sel time will be cut by 1 s on edges (= template length): [-2, 2] -> [-1, 1]
% P_level = 0.01;         % significance level
% t_siglen = 0.25;        % in [s], significance time threshold (ch is signif. if its significance > threshold) 

%% settings
load_ch2roi = true;             % load again the GRUP struct ?

%% get job-related settings (params): jobExecutor_SEI.m, v11 -> params
% --- specify analysis to be performed -> runner
runner = struct;
runner.newJobs = false; % set to 'false' if the results are already computed
runner.newFDR = false;  % true = updates the FDR correction over all patients (and freq. bands, ...) 

% --- job settings: specific (define spectral method) -> params
spectralMethod = 'v31_stft_session_bip';           % ~ v31 (BIP), similar to v12 (STFT base = sess, no zscore, no rej, FB = loFreq0)
% spectralMethod = 'v30_stft_baseRS_bip';         % ~ v24 (BIP), STFT with RS baseline
% spectralMethod = 'v24_stft_baseRS_car';         % ~ v24 (CAR), STFT with RS baseline
% spectralMethod = 'v24_stft_baseRS_bip';         % ~ v24 (BIP), STFT with RS baseline
% spectralMethod = 'stft_sessions';       % ~ v12
% spectralMethod = 'stft_baseRS';         % ~ v11  (CAR)
% spectralMethod = 'stft_baseRS_bip';     % ~ v11b (BIP)

params = msSEI_getJobParams(runner, spectralMethod, P_level);

%% =============================LOAD DATA==================================
%% ch2roi -> GRUP 
if load_ch2roi
    clear GRUP;
    groupInfo = struct(...
        'brainAtlas', {{brainAtlas}}, ... % brain atlases: 'Yeo7','Yeo7_addLat','Yeo17','Yeo17_addLat','Mars','Mars_addLat'
        'freqBands', {{'loFreq0','loGamma','hiGamma'}}, ...  % freq. bands (FB): 'delta','theta','alpha','beta','loGamma','hiGamma'
        'signifOfActiv', {{'was_sgnf'}}, ...   % significance of activation
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
end

%% assign channels of each AA to NN
fieldName_H = 'neurologyLabel';
% fieldName_H = 'ass_mars_name';
for aa = 1:size(GRUP.chSubj_inds,1)

    %% get NN labels for all channels & subj -> TBL_all
    TBL_all = [];
    c = 1;
    for subj = 1:size(params.storage.subjList,1)
        subjTag = params.storage.subjList{subj,1};

        % --- cache file
        cacheFile = [params.storage.dir_results filesep params.storage.outName filesep subjTag filesep 'cacheFile.mat'];
        assert(exist(cacheFile,'file') == 2);

        % --- load structure H & selected channels
        clear H selCh_H_resp;
        load(cacheFile, 'H', 'selCh_H_resp');   
        selCh_H = selCh_H_resp; 
        
        % --- get all NNs -> TBL_all
        for ch = 1:size(selCh_H,2)
            thisCh = selCh_H(ch);
            if strcmp(list_AA{aa}, H.channels(thisCh).(['ass_' lower(brainAtlas) '_name']))
                TBL_all{c,1} = subj;
                TBL_all{c,2} = ch;
                TBL_all{c,3} = H.channels(thisCh).(fieldName_H);
                c = c+1;
            end
        end

    end % of subj
    assert(size(TBL_all,1) >= size(GRUP.chSubj_inds{aa,1,2},1));
    
    %% unique NNs
    i_tbd = [];
    for nn = 1:size(TBL_all,1)
        if isnan(TBL_all{nn,3})
            i_tbd = cat(1, i_tbd, nn);
        end
    end
    if ~isempty(i_tbd)
        TBL_all(i_tbd,:) = [];
    end
    NN_unq = unique(TBL_all(:,3));
    TBL_nn = cell(size(NN_unq,1), 2, (1+size(GRUP.chSubj_inds,2)));
    
    %% go thru all unique NNs
    for nn = 1:size(NN_unq,1)
        NN_name = NN_unq{nn};
        for c = 1:size(GRUP.chSubj_inds{aa,1,2},1)
            this_chnl = GRUP.chSubj_inds{aa,1,2}(c,1);
            this_subj = GRUP.chSubj_inds{aa,1,2}(c,2);
            i_tbl = (cell2mat(TBL_all(:,1)) == this_subj) & (cell2mat(TBL_all(:,2)) == this_chnl);
            assert(sum(i_tbl,1) <= 1);  % found just 1 channel or no channel
            i_table = find(i_tbl == 1);  % index to TBL_all
            if ~isempty(i_table)
                if strcmp(NN_name, TBL_all{i_table,3})
                    TBL_nn{nn,1,1} = cat(1, TBL_nn{nn,1,1}, this_chnl);
                    TBL_nn{nn,2,1} = cat(1, TBL_nn{nn,2,1}, this_subj);
                end
            end
        end % of c
    end % of nn
    
    %% go throu FBs
    for freq = 1:size(GRUP.chSubj_inds,2)    
        for nn = 1:size(NN_unq,1)
            NN_name = NN_unq{nn};
            for c = 1:size(GRUP.chSubj_inds{aa,freq,1},1)    % go only thru the significant channels
                this_chnl = GRUP.chSubj_inds{aa,freq,1}(c,1);
                this_subj = GRUP.chSubj_inds{aa,freq,1}(c,2);
                i_tbl = (cell2mat(TBL_all(:,1)) == this_subj) & (cell2mat(TBL_all(:,2)) == this_chnl);
                assert(sum(i_tbl,1) <= 1);  % found just 1 channel
                i_table = find(i_tbl == 1);  % index to TBL_all
                if ~isempty(i_table)
                    if strcmp(NN_name, TBL_all{i_table,3})
                        TBL_nn{nn,1,1+freq} = cat(1, TBL_nn{nn,1,1+freq}, this_chnl);
                        TBL_nn{nn,2,1+freq} = cat(1, TBL_nn{nn,2,1+freq}, this_subj);
                    end
                end
            end % of c
        end % of nn        
    end
    
    %% sort based on most channels in NN
    nCh_nn = [];
    for nn = 1:size(NN_unq,1)
        nCh_nn = cat(1, nCh_nn, numel(TBL_nn{nn,1,1}));
    end
    [srt,i_srt] = sort(nCh_nn, 'descend');
    i_srt(srt==0) = [];
    
    %% make a table (cell array)
    TBL_export = cell(size(i_srt,1), 1+size(TBL_nn,3));
    
    %% fill table to export
    row = 1;
    for nn = i_srt'
        col = 1;    % nn name
        TBL_export{row,col} = NN_unq{nn};
        
        col = 2;    % total number of channels
        nC = numel(TBL_nn{nn,1,1});
        nS = numel(unique(TBL_nn{nn,2,1}));
        TBL_export{row,col} = [num2str(nC) ' (' num2str(nS) ')'];  % total number of channels
        
        % FB channels
        for freq = 1:size(GRUP.chSubj_inds,2) 
            col = col+1;
            nC = numel(TBL_nn{nn,1,1+freq});
            nS = numel(unique(TBL_nn{nn,2,1+freq}));
            TBL_export{row,col} = [num2str(nC) ' (' num2str(nS) ')'];  % total number of channels    
        end
        row = row+1;
    end
        
    %% export to excell
    outDir = [params.storage.dir_results filesep params.storage.outName filesep 'ch2roi' groupInfo.dirSuffix filesep anatAtlas_name filesep 'msSEI_txt_AA_YeoNN'];
    if ~exist(outDir, 'dir')
        mkdir(outDir);
    end 
    fileName = ['TBL_' brainAtlas '_into_' fieldName_H '.xlsx'];
    xlswrite([outDir filesep fileName],TBL_export,list_AA{aa});
    disp(['Saved sheet: ' list_AA{aa}]);
end
disp(['Table saved into: ' [outDir filesep fileName]]);
        