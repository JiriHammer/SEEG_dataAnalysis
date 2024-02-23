function GRUP = ch2roi_load_FB(params, dataInfo)
% !!! TBD !!!
% pools results from trigger analysis over:
%   - freq. bands
%   - ROIs (regions of interest): anatomic areas, ...
%   - subjects
% returns:
% D{roi,freq} = 3D: t x ch x clz

% (c) Jiri, Nov21

%% required variables
list_AA = dataInfo.list_AA;     % list of ROIs (anatomic areas)
list_FB = dataInfo.list_FB;     % list of freq. bands (delta, .., beta, .., hiGamma)
list_anatomy_signif = dataInfo.list_anatomy_signif; % cell: 1 x 3, where 3 = list of AA - FB - significance: 

%% default output = all
if ~isfield(dataInfo, 'outVars')
    dataInfo.outVars = {'time','trials','sgnfLags','nChnls_roi','nSubj_roi','chnls_MNI','chnls_hVals', 'SNR_trials', 'spectra'};
end

%% allocation of output
D_trials = cell(size(list_AA,1), size(list_FB,1));       % 2D cell: aa x freq, grouped data, where D{aa,freq} = 3D: t x ch x clz
D_time = cell(size(list_AA,1), size(list_FB,1));         % 2D cell: aa x freq, time axis
TIME_signif = cell(size(list_AA,1), size(list_FB,1));    % 2D cell: aa x freq: count of significant time lags over channels 
CHNL_roi = zeros(size(list_AA,1), size(list_FB,1));      % 2D array: aa x freq: count channels in ROI (no matter if significant)
SUBJ_roi = cell(size(list_AA,1), size(list_FB,1),2);     % 3D cell: aa x freq x 2: catted subjects in ROI, 3rd dim = [signif, all] - use 'unique' to get subj count
GRUP_hVals = cell(size(list_AA,1), size(list_FB,1));     % 2D cell: aa x freq, significance of difference between 2 classes, where D{aa,freq} = 2D: t x ch
GRUP_mni   = cell(size(list_AA,1), size(list_FB,1));     % 2D cell: aa x freq, MNI coors of channels, where D{aa,freq} = 1D: ch x 3
SNR_trials = cell(size(list_AA,1), size(list_FB,1));     % 2D cell: aa x freq, SNR of "tuning curve" (var of means / mean of vars)  over classes, where D{aa,freq} = 2D: t x ch
D_spectra = cell(size(list_AA,1), 1);                    % 1D cell: aa x 1, S{aa} = 3D: freq,time,clz
freqAxis = [];

%% pool over subjects
for subj = 1:size(params.storage.subjList,1)
    subjTag = params.storage.subjList{subj,1};

    %% cache file
    cacheFile = [params.storage.dir_results filesep params.storage.outName filesep subjTag filesep 'cacheFile.mat'];
    assert(exist(cacheFile,'file') == 2);

    %% load structure H & selected channels
    clear H selCh_H_resp;
    load(cacheFile, 'H', 'selCh_H_resp');    
    selCh_H = selCh_H_resp;

    %% needed??? update headers H for yeo7 & yeo17 atlas info
%     mni_center = 'singleCh';
%     updateHeaders_atlasInfo(params, selCh_H_resp, H, cacheFile, mni_center);

    %% anatomical assignements from isarg_atlas
    if size(selCh_H,1) == 2
        mni_center = 'chPair'; 
    else
        mni_center = 'singleCh';
    end
    [ass_isargAtlas, H] = updateHeaders_atlasInfo(params, selCh_H_resp, H, cacheFile, mni_center);

    %% define freq band
    for freq = 1:size(list_FB,1)
        freqBand = list_FB{freq,1};
        assert(ismember(freqBand, list_FB(:,1)));
        
        %% load freq. band activations (e.g. trialsData_beta) -> trialsData
        trialsData = [];
        varName = ['trialsData_' freqBand]; 
        % try to load freq. band activation from cache file
        if ~ismember(varName, who('-file', cacheFile))   % if the variable 'varName' exists
            error(['WARNING: subject = ' subjTag ', variable = ' varName ' was not found in cacheFile: ' cacheFile]);
        else
            clear trialsData;
            load(cacheFile, varName);
            assert(exist(varName,'var') == 1);
            eval(['trialsData = ' varName ';']); 
            clear(varName); 
        end
        assert(size(selCh_H,2) == size(trialsData.hVals,2));
        
        %% load spectra (e.g. spectralData_clz1) -> S_data
        S_data = [];
        for clz = 1:size(trialsData.yVals,3)
            varName = ['spectralData_clz' num2str(clz)];
            if ~ismember(varName, who('-file', cacheFile))   % if the variable 'varName' exists
                disp(['WARNING: subject = ' subjTag ', variable = ' varName ' was not found in cacheFile: ' cacheFile]);
            else
                load(cacheFile, varName);
                assert(exist(varName,'var') == 1);
                eval(['S = ' varName ';']); 
                clear(varName); 
                S_data = cat(4, S_data, S.cVals);     % 4D: freq,time,ch,clz
            end
        end
        S_time = S.xVals;
        S_freq = S.yVals;
        assert(size(selCh_H,2) == size(S_data,3));
        
        %% go through selected ROIs (aa = anatomic areas)
        for aa = 1:size(list_AA,1)
            sel_ROI = list_AA{aa,1};          
            if subj == 1
                TIME_signif{aa,freq} = zeros(size(trialsData.xVals,1),1);
                D_time{aa,freq} = trialsData.xVals;
            else
                assert(size(D_time{aa,freq},1) == size(trialsData.xVals,1));
            end
            
            %% go through selected channels
            for ch = 1:size(selCh_H,2)         
                
                %% --- significant? ---
                useCh_significance = 0;    % if to include the channel based on significance criteria
                if strcmp(list_anatomy_signif{1,3},'was_sgnf') && any(trialsData.hVals(:,ch) == 1)
                    useCh_significance = 1;    % 'was_sgnf';
                elseif strcmp(list_anatomy_signif{1,3},'not_sgnf') && all(trialsData.hVals(:,ch) ~= 1)
                    useCh_significance = 1;    % 'not_sgnf'; 
                elseif strcmp(list_anatomy_signif{1,3},'all')
                    useCh_significance = 1;    % 'all';
                end

                %% --- anatomic area assignment ---  
                if strcmp(list_anatomy_signif{1,1},'Yeo7_networks')   
                    thisAA = anatomicalArea_getName('yeo7_name', list_AA, ass_isargAtlas(ch));
                elseif strcmp(list_anatomy_signif{1,1},'Yeo17_networks')   
                    thisAA = anatomicalArea_getName('yeo17_name', list_AA, ass_isargAtlas(ch));
                elseif strcmp(list_anatomy_signif{1,1},'neurologists_grouping')   
                    thisAA = anatomicalArea_getName('neurologists_grouping', list_AA, H.channels(selCh_H(ch)));
                    thisAA = thisAA{1};     % TO DO, now takes only the first assignment
                end

                %% >>>> !!! group channel !!! <<<<
                if strcmp(thisAA, sel_ROI) && useCh_significance
                    D_trials{aa,freq} = cat(2, D_trials{aa,freq}, trialsData.yVals(:,ch,:));
                    SUBJ_roi{aa,freq,1} = cat(1, SUBJ_roi{aa,freq,1}, subj);
                    i_t = find(trialsData.hVals(:,ch) == 1);
                    if ~isempty(i_t)
                        TIME_signif{aa,freq}(i_t) = TIME_signif{aa,freq}(i_t) + 1;    % signif. times
                    end
                    D_spectra{aa} = cat(4, D_spectra{aa}, squeeze(S_data(:,:,ch,:)));
                end

                %% selected ROI (no matter if significant)
                if strcmp(thisAA, sel_ROI)
                    CHNL_roi(aa,freq) = CHNL_roi(aa,freq) + 1;                  % count channels
                    SUBJ_roi{aa,freq,2} = cat(1, SUBJ_roi{aa,freq,2}, subj);    % cat subj index
                    GRUP_mni{aa,freq} = cat(1, GRUP_mni{aa,freq}, [H.channels(selCh_H(ch)).MNI_x, H.channels(selCh_H(ch)).MNI_y, H.channels(selCh_H(ch)).MNI_z]);                   % MNI coors
                    GRUP_hVals{aa,freq} = cat(2, GRUP_hVals{aa,freq}, trialsData.hVals(:,ch)); 
                    SNR_trials{aa,freq} = cat(2, SNR_trials{aa,freq}, var(nanmean(trialsData.yVals,2),0,3)./mean(nanvar(trialsData.yVals,0,2),3));    % 1D: t x 1
                end

            end     % end of ch
        end     % end of aa
    end     % end of subj
end     % end of freq

%% output other variables -> GRUP_info
GRUP = struct;

% time axis of grouped trials data
if ismember('time', dataInfo.outVars)
    GRUP.time = D_time;     % 2D cell: aa x freq, where D{aa,freq} = 1D: t x 1
end

% grouped trials data
if ismember('trials', dataInfo.outVars)
    GRUP.trials = D_trials;     % 2D cell: aa x freq, where D{aa,freq} = 3D: t x ch x clz
end

% count of significant time lags over grouped channels 
if ismember('sgnfLags', dataInfo.outVars)
    GRUP.sgnfLags = TIME_signif;    % 2D cell: aa x freq
end

% count of channels in ROI (no matter if significant) 
if ismember('nChnls_roi', dataInfo.outVars)
    GRUP.nChnls_roi = CHNL_roi;    % 2D cell: aa x freq
end

% count of suibjects in ROI (signif. & all) 
if ismember('nSubj_roi', dataInfo.outVars)
    GRUP.nSubj_roi = SUBJ_roi;    % 3D cell: aa x freq x 2: 3rd dim = [signif, all]
end

% MNI coors of channels in ROI (no matter if significant) 
if ismember('chnls_MNI', dataInfo.outVars)
    GRUP.chnls_MNI = GRUP_mni;    % 2D cell: aa x freq, where D{aa,freq} = 1D: ch x 3
end

% significance of difference between 2 classes (works only for two classes!) 
if ismember('chnls_hVals', dataInfo.outVars)
    GRUP.chnls_hVals = GRUP_hVals;    % 2D cell: aa x freq, where D{aa,freq} = 2D: t x ch
end

% SNR of "tuning curve" (var of means / mean of vars)  over classes
if ismember('SNR_trials', dataInfo.outVars)
    GRUP.SNR_trials = SNR_trials;    % 2D cell: aa x freq, where D{aa,freq} = 2D: t x ch
end

% spectra (mean over channels)
if ismember('spectra', dataInfo.outVars)
    GRUP.spectra_data = cell(size(list_AA,1), 1);
    for aa = 1:size(list_AA,1)
        if ~isempty(D_spectra{aa})
            GRUP.spectra_data{aa} = nanmean(D_spectra{aa},4);    % 1D cell: aa x 1, where D{aa} = 3D: freq x t x clz
        end
    end
    GRUP.spectra_freq = S_freq;
    GRUP.spectra_time = S_time;
    GRUP.spectra_snr = cell(size(list_AA,1), 1);
    for aa = 1:size(list_AA,1)
        GRUP.spectra_snr{aa} = var(nanmean(D_spectra{aa},4),0,3)./mean(nanvar(D_spectra{aa},0,4),3);    % SNR over channels, 2D: f x t
    end
end

% class names
GRUP.info.clzNames = trialsData.info.clzNames;     % assumes same for all FBs
GRUP.info.clzColor = trialsData.info.colors;     % assumes same for all FBs
