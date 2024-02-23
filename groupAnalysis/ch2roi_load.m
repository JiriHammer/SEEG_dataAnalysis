function GRUP = ch2roi_load(params, groupInfo)
% pools results from trigger analysis over:
%   - freq. bands
%   - ROIs (regions of interest): anatomic areas, ...
%   - subjects
% spectra are evaluated based on significance criteria of the last FB specified (typically hiGamma)
% returns: GRUP struct

% (c) Jiri, Nov21

disp('Grouping channels activations ...');

%% required variables
list_AA = groupInfo.list_AA;     % list of ROIs (anatomic areas)
list_FB = groupInfo.list_FB;     % list of freq. bands (delta, .., beta, .., hiGamma)
list_anatomy_signif = groupInfo.list_anatomy_signif; % cell: 1 x 3, where 3 = list of AA - FB - significance: 
brainAtlas_name = list_anatomy_signif{1,1};

%% default variables
if ~isfield(groupInfo, 'significanceLevel')
    significanceLevel = 0.001;      % significance on group level
else
    significanceLevel = groupInfo.significanceLevel;
end

%% allocation of output
D_trials = cell(size(list_AA,1), size(list_FB,1));          % 2D cell: aa x freq, grouped data, where D{aa,freq} = 3D: t x ch x clz
D_time = cell(size(list_AA,1), size(list_FB,1));            % 2D cell: aa x freq, time axis
TIME_signif = cell(size(list_AA,1), size(list_FB,1));       % 2D cell: aa x freq: count of significant time lags over channels 
CHNL_roi = zeros(size(list_AA,1), size(list_FB,1));         % 2D array: aa x freq: count channels in ROI (no matter if significant)
SUBJ_roi = cell(size(list_AA,1), size(list_FB,1),2);        % 3D cell: aa x freq x 2: catted subjects in ROI, 3rd dim = [signif, all] - use 'unique' to get subj count
GRUP_hVals = cell(size(list_AA,1), size(list_FB,1));        % 2D cell: aa x freq, significance of difference between 2 classes, where D{aa,freq} = 2D: t x ch
GRUP_sel_hVals = cell(size(list_AA,1), size(list_FB,1));    % 2D cell: aa x freq, significance of difference between 2 classes for selected channels, where D{aa,freq} = 2D: t x ch
GRUP_mni_all   = cell(size(list_AA,1), size(list_FB,1));    % 2D cell: aa x freq, MNI coors of channels, where D{aa,freq} = 1D: ch x 3, all channels in ROI (no matter of significant activation)
GRUP_mni_sel   = cell(size(list_AA,1), size(list_FB,1));    % 2D cell: aa x freq, MNI coors of channels, where D{aa,freq} = 1D: ch x 3, selected channels in ROI (based on significe of activation)
GRUP_chSubj_inds= cell(size(list_AA,1), size(list_FB,1),2); % 3D cell: aa x freq x 2: 3rd dim = [signif, all], where GRUP.ind_chSubj{aa,freq,1} = [ch, subj], where ch = index to GRUP.trials{aa,freq}(t,ch,clz) or selCh_H_resp
GRUP_chSubj_name= cell(size(list_AA,1), size(list_FB,1),2); % 3D cell: aa x freq x 2: 3rd dim = [signif, all], where GRUP.ind_chSubj{aa,freq,1} = {chName, subjTag}
SNR_trials = cell(size(list_AA,1), size(list_FB,1));     % 2D cell: aa x freq, SNR of "tuning curve" (var of means / mean of vars)  over classes, where D{aa,freq} = 2D: t x ch
D_spectra = cell(size(list_AA,1), 1);                    % 1D cell: aa x 1, S{aa} = 3D: freq,time,clz
paraTimes = cell(size(list_AA,1), size(list_FB,1));      % 2D cell: aa x freq, paradigm times (RT, ...) for each channel (!), 3D: clz, responseType (RTs), ch
H_trials = cell(size(list_AA,1), size(list_FB,1));       % 2D cell: aa x freq, significance between 2 classes, where D{aa,freq} = 3D: t x 1 
P_trials = cell(size(list_AA,1), size(list_FB,1));       % 2D cell: aa x freq, p-values     between 2 classes, where D{aa,freq} = 3D: t x 1 
BASE_vals = cell(size(list_AA,1), size(list_FB,1));       % 2D cell: aa x freq, p-values     between 2 classes, where D{aa,freq} = 3D: t x 1 

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

    %% anatomical assignements from isarg_atlas
    [ass_isargAtlas, H] = updateHeaders_atlasInfo(params, selCh_H, H, cacheFile);

    %% define freq band
    for freq = 1:size(list_FB,1)
        freqBand = list_FB{freq,1};
        assert(ismember(freqBand, list_FB(:,1)));
        disp([' - loading trialsData, subj = ' subjTag ', freq = ' freqBand ' ...']);
        
        %% load freq. band activations (e.g. trialsData_beta) -> trialsData
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
            if ~isfield(trialsData.info, 'paraTimes')
                trialsData.info.paraTimes = nan;        % required field
            end
        end
        assert(size(selCh_H,2) == size(trialsData.hVals,2));
        
        %% load baseline values (from RS) 
        BASE_freq = [];
        varName = 'baseData';
        if ismember(varName, who('-file', cacheFile))   % if the variable 'varName' exists
            clear(varName);
            load(cacheFile, varName);       
            assert(exist(varName,'var') == 1);
            eval(['B = ' varName ';']);     
            clear(varName);
            [tf, i_fb] = ismember(freqBand, params.triggering.freqBands(:,1));
            assert(tf);
            i_freq = closestval(B.xVals, params.triggering.freqBands{i_fb,2}(1)):closestval(B.xVals, params.triggering.freqBands{i_fb,2}(2));
            BASE_freq = [mean(B.yVals(i_freq,:),1); mean(B.yErrs(i_freq,:),1)];    % 2D: ch x 2 (2 = mean & SD)
        end
        
        %% load spectra (e.g. spectralData_clz1) -> S_data
        if freq == 1                                % load only once
            S_data = [];
            for clz = 1:size(trialsData.yVals,3)
                varName = ['spectralData_clz' num2str(clz)];
                if ~ismember(varName, who('-file', cacheFile))   % if the variable 'varName' exists
                    disp(['WARNING: subject = ' subjTag ', variable = ' varName ' was not found in cacheFile: ' cacheFile]);
                else
                    disp([' - loading spectralData, subj = ' subjTag ', clz = ' num2str(clz) ' ...']);
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
        end
        
        %% go through selected ROIs (aa = anatomic areas)
        disp([' - assignment to ROIs, subj = ' subjTag ', freq = ' freqBand ' ...']);
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
                chName = H.channels(selCh_H(ch)).name;
                
                %% --- significant? ---
                useCh_significance = ch2roi_getSignificance(groupInfo, list_anatomy_signif{1,3}, trialsData.hVals(:,ch), trialsData.xVals);
                
                %% --- anatomic area assignment ---  
                thisAA = ch2roi_assignAA(brainAtlas_name, ass_isargAtlas(ch), groupInfo, list_AA);
                
                %% --- grey/white matter assignment ---  
                useCh_GW = ch2roi_assignGW(list_anatomy_signif{1,4}, ass_isargAtlas(ch), groupInfo.GW_threshold);
                
                %% >>>> !!! group channel !!! <<<<
                if strcmp(thisAA, sel_ROI) && useCh_significance && useCh_GW
                    D_trials{aa,freq} = cat(2, D_trials{aa,freq}, trialsData.yVals(:,ch,:));
                    SUBJ_roi{aa,freq,1} = cat(1, SUBJ_roi{aa,freq,1}, subj);  % index of subj & ch
                    SNR_trials{aa,freq} = cat(2, SNR_trials{aa,freq}, var(squeeze(trialsData.yVals(:,ch,:)),0,2)./mean(repmat(trialsData.info.nClz, [size(trialsData.yErrs,1),1,1]) .* squeeze(trialsData.yErrs(:,ch,:)).^2,2));    % 1D: t x 1, SNR = var(clz means)/mean(clz vars)
                    i_t = find(trialsData.hVals(:,ch) == 1);
                    if ~isempty(i_t)
                        TIME_signif{aa,freq}(i_t) = TIME_signif{aa,freq}(i_t) + 1;    % signif. times
                    end
                    if freq == size(list_FB,1)  % !!! only for the last FB !!!
                        D_spectra{aa} = cat(4, D_spectra{aa}, squeeze(S_data(:,:,ch,:)));
                    end
                    if ~isempty(BASE_freq)
                        BASE_vals{aa,freq} = cat(2, BASE_vals{aa,freq}, BASE_freq(:,ch));
                    end
                    GRUP_mni_sel{aa,freq} = cat(1, GRUP_mni_sel{aa,freq}, [H.channels(selCh_H(ch)).MNI_x, H.channels(selCh_H(ch)).MNI_y, H.channels(selCh_H(ch)).MNI_z]);                   % MNI coors
                    GRUP_chSubj_inds{aa,freq,1} = cat(1, GRUP_chSubj_inds{aa,freq,1}, [ch, subj]);              % channel, subj indices of selected channels
                    GRUP_chSubj_name{aa,freq,1} = cat(1, GRUP_chSubj_name{aa,freq,1}, {chName, subjTag});       % channel, subj names   of selected channels
                    GRUP_sel_hVals{aa,freq} = cat(2, GRUP_sel_hVals{aa,freq}, trialsData.hVals(:,ch)); 
                    paraTimes{aa,freq} = cat(3, paraTimes{aa,freq}, trialsData.info.paraTimes);                 % paradigm times (RT, ...), paraTimes{aa,freq} = 3D: clz, responseType (RTs), ch
                end

                %% selected ROI (no matter if significant)
                if strcmp(thisAA, sel_ROI)
                    CHNL_roi(aa,freq) = CHNL_roi(aa,freq) + 1;                  % count channels
                    SUBJ_roi{aa,freq,2} = cat(1, SUBJ_roi{aa,freq,2}, [subj, ch]);    % cat subj index
                    GRUP_mni_all{aa,freq} = cat(1, GRUP_mni_all{aa,freq}, [H.channels(selCh_H(ch)).MNI_x, H.channels(selCh_H(ch)).MNI_y, H.channels(selCh_H(ch)).MNI_z]);                   % MNI coors
                    GRUP_chSubj_inds{aa,freq,2} = cat(1, GRUP_chSubj_inds{aa,freq,2}, [ch, subj]);              % channel, subj indices of all channels in ROI
                    GRUP_chSubj_name{aa,freq,2} = cat(1, GRUP_chSubj_name{aa,freq,2}, {chName, subjTag});       % channel, subj names   of all channels in ROI
                    GRUP_hVals{aa,freq} = cat(2, GRUP_hVals{aa,freq}, trialsData.hVals(:,ch)); 
%                     SNR_trials{aa,freq} = cat(2, SNR_trials{aa,freq}, var(nanmean(trialsData.yVals,2),0,3)./mean(nanvar(trialsData.yVals,0,2),3));    % 1D: t x 1
                end

            end     % end of ch
        end     % end of aa
    end     % end of freq   
end     % end of subj

%% output -> GRUP
GRUP = struct;

% time axis of grouped trials data
GRUP.time = D_time;     % 2D cell: aa x freq, where D{aa,freq} = 1D: t x 1

% grouped trials data
GRUP.trials = D_trials;     % 2D cell: aa x freq, where D{aa,freq} = 3D: t x ch x clz

% significances between classes (works only for 2 classes)
significanceLevel = 0.001;  % quick hack for msSEI (msSEI_fig_AA.m) (now=Jan24)
for aa = 1:size(H_trials,1)
    for freq = 1:size(H_trials,2)
        if size(D_trials{aa,freq},3) == 2
            vals = reshape(D_trials{aa,freq}, [size(D_trials{aa,freq},1), 1, size(D_trials{aa,freq},2)*size(D_trials{aa,freq},3)]); % 3D: t x 1 x chnls-labels, where chnls-labels = [chnls-clz1, chnls-clz2]
            labels = cat(1, ones(size(D_trials{aa,freq},2),1), 2*ones(size(D_trials{aa,freq},2),1));    % hardcoded for 2 classes, labels = 1 & 2
            [P_trials{aa,freq}, H_trials{aa,freq}] = getSignificance(vals, labels, significanceLevel);
        end
    end
end
GRUP.H_trials = H_trials;
GRUP.P_trials = P_trials;

% count of significant time lags over grouped channels 
GRUP.sgnfLags = TIME_signif;    % 2D cell: aa x freq

% TBD???: trials: subj mean -> GRUP.trials_subjMean
D_trials_subjMean = cell(size(list_AA,1), size(list_FB,1));          % 2D cell: aa x freq, grouped data, where D{aa,freq} = 3D: t x ch x clz
for aa = 1:size(D_trials,1)
    for freq = 1:size(D_trials,2)
        for subj = 1:size(params.storage.subjList,1)
            i_subj = [];
            if ~isempty(GRUP_chSubj_inds{aa,freq,1})
                i_subj = find(GRUP_chSubj_inds{aa,freq,1}(:,2) == subj);
            end
            if ~isempty(i_subj)
                D_trials_subjMean{aa,freq} = cat(2, D_trials_subjMean{aa,freq}, mean(D_trials{aa,freq}(:,i_subj,:),2));
            end
        end
    end
end
GRUP.trials_subjMean = D_trials_subjMean;

% TBD???: significances between classes (works only for 2 classes) of trials_subjMean
H_trials_subjMean = cell(size(list_AA,1), size(list_FB,1));       % 2D cell: aa x freq, significance between 2 classes, where D{aa,freq} = 3D: t x 1 
P_trials_subjMean = cell(size(list_AA,1), size(list_FB,1));       % 2D cell: aa x freq, p-values     between 2 classes, where D{aa,freq} = 3D: t x 1 
for aa = 1:size(D_trials_subjMean,1)
    for freq = 1:size(D_trials_subjMean,2)
        if size(D_trials_subjMean{aa,freq},3) == 2  % works only for 2 classes
            vals = reshape(D_trials_subjMean{aa,freq}, [size(D_trials_subjMean{aa,freq},1), 1, size(D_trials_subjMean{aa,freq},2)*size(D_trials_subjMean{aa,freq},3)]); % 3D: t x 1 x subj-labels, where subj-labels = [subj-clz1, subj-clz2]
            labels = cat(1, ones(size(D_trials_subjMean{aa,freq},2),1), 2*ones(size(D_trials_subjMean{aa,freq},2),1));    % hardcoded for 2 classes
            [P_trials_subjMean{aa,freq}, H_trials_subjMean{aa,freq}] = getSignificance(vals, labels, significanceLevel);
        end
    end
end
GRUP.H_trials_subjMean = H_trials_subjMean;
GRUP.P_trials_subjMean = P_trials_subjMean;

% count of channels in ROI (no matter if significant) 
GRUP.nChnls_roi = CHNL_roi;    % 2D cell: aa x freq

% count of subjects in ROI (signif. & all) 
GRUP.nSubj_roi = SUBJ_roi;    % 3D cell: aa x freq x 2: 3rd dim = [signif, all]

% MNI coors of channels in ROI (no matter if significant) 
GRUP.chnls_MNI_all = GRUP_mni_all;    % 2D cell: aa x freq, where D{aa,freq} = 1D: ch x 3

% MNI coors of channels in ROI (using significance criteria) 
GRUP.chnls_MNI_sel = GRUP_mni_sel;    % 2D cell: aa x freq, where D{aa,freq} = 1D: ch x 3

% all channels: significance of difference between 2 classes (works only for two classes!) 
GRUP.chnls_hVals = GRUP_hVals;    % 2D cell: aa x freq, where D{aa,freq} = 2D: t x ch

% selected channels: significance of difference between 2 classes (works only for two classes!) 
GRUP.sel_hVals = GRUP_sel_hVals;    % 2D cell: aa x freq, where D{aa,freq} = 2D: t x ch

% selected channels: channel & subj indices 
GRUP.chSubj_inds = GRUP_chSubj_inds;    % 3D cell: aa x freq x 2: 3rd dim = [signif, all], where GRUP.ind_chSubj{aa,freq,1} = [subj, ch], where ch = index to GRUP.trials{aa,freq}(t,ch,clz) or selCh_H_resp

% selected channels: channel & subj indices 
GRUP.chSubj_name = GRUP_chSubj_name;    % 3D cell: aa x freq x 2: 3rd dim = [signif, all], where GRUP.ind_chSubj{aa,freq,1} = [subj, ch], where ch = index to GRUP.trials{aa,freq}(t,ch,clz) or selCh_H_resp


% SNR of "tuning curve" (var of means / mean of vars)  over classes
GRUP.SNR_trials = SNR_trials;    % 2D cell: aa x freq, where D{aa,freq} = 2D: t x ch

% spectra (mean over channels)
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

% info: class names
GRUP.info.clzNames = trialsData.info.clzNames;     % assumes same for all FBs
GRUP.info.clzColor = trialsData.info.colors;     % assumes same for all FBs
GRUP.info.txtOnFig = [...
    'brain atlas = ' brainAtlas_name ...
    ', significance = ' list_anatomy_signif{1,3} ...
    ', GM/WM = ' list_anatomy_signif{1,4} ...
    ', classes = ' ...
    ];
for clz = 1:size(trialsData.yVals,3)
    clr = trialsData.info.colors(clz,:);
    GRUP.info.txtOnFig = [GRUP.info.txtOnFig ...
        '\color[rgb]{' num2str(clr(1)) ' ' num2str(clr(2)) ' ' num2str(clr(3)) '}' trialsData.info.clzNames{clz,1} ' '];        
end

% paradigm times (RT, ...)
GRUP.paraTimes = paraTimes;

% GM / WM selection
GRUP.info.sel_gw = list_anatomy_signif{1,4};

% baseline values (from RS)
for aa = 1:size(BASE_vals,1)
    for freq = 1:size(BASE_vals,2)
        if ~isempty(BASE_vals{aa,freq})
            BASE_vals{aa,freq} = [mean(BASE_vals{aa,freq}(1,:),2);  ... % mean over ch-means
                                  sqrt(mean(BASE_vals{aa,freq}(2,:).^2,2))];  % mean over ch-SDs
        end
    end
end   
GRUP.BASE_vals = BASE_vals;
