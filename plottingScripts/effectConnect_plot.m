function [FB_y, S, freqAxis, timeAxis] = effectConnect_plot(subjTag)
% standalone script for plotting effective connectivity results (DTF)
% called by: effectConnect_plot_avgOverSubj.m
% WARNING: quick-fix exploratory script (to get mean DTF over subj)

% (c) Jiri, Sep22

%% batch job (uncomment function and use as a script)
% % set subjects
% params_default;
% params.paradigm.usedParadigm = 'switchEI';
% params.paradigm.specificType = 'switchEI';
% params.storage.dir_shareData = 'G';
% [params.storage.pathBeg, params.storage.subjList] = get_subjectList(params);
% subjList = params.storage.subjList;
% % subjList = { 
% %     '20_PR3'; ...  
% % };
% 
% % plot
% for subj = 1:size(subjList,1)
%     subjTag = subjList{subj,1};
%     effectConnect_plot;
% end

%% load DTF
% if ~exist('subjTag','var'), subjTag = '20_PR4'; end
pathBeg = 'G:\dox\ms_switch_EI\data\v17_dtf_test\switchin_EI_IE_car';
cacheFile = [pathBeg filesep subjTag filesep 'cacheFile.mat'];
assert(exist(cacheFile,'file') == 2);
clear EC_info DTF trialsData_erp;
load(cacheFile, 'EC_info', 'DTF', 'trialsData_erp', 'params', 'H');

%% required variables
clzLabels = EC_info.clzLabels;  
S = [];
FB_y = [];
freqAxis = [];
timeAxis = [];

%% what to plot?
plotEC_ch2ch_spectra = false;
plotEC_ch2ch_freqBands = false;
plotEC_selGroups_spectra = true;
plotEC_selGroups_freqBands = true;

%% channel selection to ROI -> selCh_groups{roi}
params.connectivity.selectedChnls = {...    % each row = tuple: atlasName, areaName (must match names in anatomicalAreas_getList.m) 
    'Yeo7','Default'; ...
    'Yeo7','Dorsal Attention'; ...  
    };  
[COH_selCh, selCh_groups] = ch2roi_selCh_atlasRoi(params, size(trialsData_erp.yVals,2));

% COH_selCh = [];
% selCh_groups = cell(size(params.connectivity.selectedChnls,1),1);
% % go for each ROI (= 1 row in params.connectivity.selectedChnls)
% for roi = 1:size(params.connectivity.selectedChnls,1)
%     brainAtlas_name = params.connectivity.selectedChnls{roi,1};
%     sel_ROI = params.connectivity.selectedChnls{roi,2};
%     if strcmp(brainAtlas_name, 'all')       % special case: all channels
%         selCh_groups{roi} = 1:size(trials.data,2);
%         COH_selCh = 1:size(trials.data,2);
%     else    
%         % load H & selCh_H
%         load(params.storage.cacheFile, 'H', 'selCh_H_resp');
%         selCh_H = selCh_H_resp;
%         assert(size(trialsData_erp.yVals,2) == size(selCh_H,2));
%         
%         % anatomical assignements from isarg_atlas
%         if size(selCh_H,1) == 2
%             mni_center = 'chPair'; 
%         else
%             mni_center = 'singleCh';
%         end
%         [ass_isargAtlas, H] = updateHeaders_atlasInfo(params, selCh_H_resp, H, params.storage.cacheFile, mni_center);
% 
%         % assign each channel
%         for ch = 1:size(selCh_H,2)
%             % --- anatomic area assignment ---  
%             thisAA = ch2roi_assignAA(brainAtlas_name, ass_isargAtlas(ch), params.connectivity);                      
%             if strcmp(thisAA, sel_ROI)                  % is it the selected ROI?
%                 selCh_groups{roi} = cat(2, selCh_groups{roi}, ch);
%                 COH_selCh = cat(2, COH_selCh, ch);
%             end
%         end
% 
%     end
% end
% COH_selCh = sort(unique(COH_selCh));
% if isempty(COH_selCh)
%     disp('No channels found for connectivity computation. Check selection in: params.connectivity.selectedChnls');
%     return;
% end

%% plot connectivity: selected groups of channels (each group to each group)
if plotEC_selGroups_spectra
    COH_selCh = EC_info.selCh;
    
    %% number of non-empty channel groups
    N_groups = 0;
    i_selCh_groups = cell(size(selCh_groups,1),1);      % indices of channels in COH_selCh
    for g = 1:size(selCh_groups,1)
        if ~isempty(selCh_groups{g})
            N_groups = N_groups+1;
            for ch = 1:size(selCh_groups{g},2)
                i_ch = find(COH_selCh == selCh_groups{g}(ch));
                i_selCh_groups{g} = cat(2, i_selCh_groups{g}, i_ch);
            end
        end
    end
    L = ones(N_groups,N_groups);
%     L = tril(ones(N_groups,N_groups));
%     L = L(:);

    if N_groups < 2
        disp('Not enough groups found for grouped connectivity plotting!');
        return;
    end

    %% init struct spectra (for connectivity plotting)
    freqAxis = EC_info.freq;     % assumes row vector, 1 x N
    timeAxis = EC_info.time;
    cLims = [-0.1, 0.1];
    clzLabels = EC_info.clzLabels;  
    S = cell(N_groups*N_groups, size(clzLabels,2));

    %% === (2a) imagesc connectivity a-la spectra: plot for each group =======
    vars2plot = {'DTF'};     % variables to plot
    for v = 1:size(vars2plot,2)
        disp(['connectivity analysis: plotting: ' vars2plot{v} ' ...']);
        M = eval(vars2plot{v});
        assert(size(M,5) == size(clzLabels,2));

        %% figure
        f = fig_make;
        nCols = numel(find(L == 1));  % each group to each group
        nRows = size(clzLabels,2);  % classes below each other
        nPlot = 1;

        % group-by-group connectivity -> C
        for clz = 1:size(clzLabels,2)
            g = 1;
            for ch_i = 1:size(selCh_groups,1)
                selCh_i = i_selCh_groups{ch_i};         % indices of channels in COH matrix
                groupName_i = params.connectivity.selectedChnls{ch_i,2};
                if ~isempty(selCh_i)
                    for ch_j = 1:size(selCh_groups,1)
                        selCh_j = i_selCh_groups{ch_j}; % indices of channels in COH matrix
                        groupName_j = params.connectivity.selectedChnls{ch_j,2};
                        if ~isempty(selCh_j)
                            % selected channels
                            nCh = size(M,1);
                            i_chMatrix = false(nCh,nCh);                % init
                            i_chMatrix(selCh_i,selCh_j) = true;         % selected channels (symmetric along diagonal, incl. diagonal)
                            i_subdiagonal = tril(ones(nCh,nCh),-1);     % below main diagonal
                            i_chMatrix = i_subdiagonal & i_chMatrix;    % deselect the diagonal
                            inds2use = repmat(i_chMatrix, [1,1,size(M,3),size(M,4)]);

                            % data
                            C = M(:,:,:,:,clz);                         % 4D: ch x ch x F x T 
                            C(~inds2use) = nan;
                            C = squeeze(nanmean(nanmean(C,2),1));       % 2D: F x T
                            S{g,clz} = C;
                            g = g+1;

                            % plot
                            subplot(nRows, nCols, nPlot);
                            hold on;
                            h = imagesc(timeAxis, freqAxis, C, cLims);
                            colormap(gca,brewermap(256,'*RdBu'));   
                            axis tight;
                            title({[vars2plot{v} ': ' groupName_i ', ' groupName_j ]; ['clz = ' trialsData_erp.info.clzNames{clz,1}]});
                            xlabel('time (s)');
                            ylabel('freq (Hz)');
                            xticks = [-2:2];
                            for t = xticks
                                plot([t t], ylim, '--k');   % time ticks
                            end
                            yticks = [10, 30, 50, 100, 120];
                            for t = yticks
                                plot(xlim, [t t], '--k');   % freq ticks
                            end
                            nPlot = nPlot+1;
                            
                        end
                    end
                end
            end
        end

        %% save figure
%         figname = [vars2plot{v} '_chGroups_spectra'];
%         outDir = [params.storage.outputDir filesep 'coherency_channelGroups'];
        figname = [subjTag '_' vars2plot{v} '_chGroups_spectra'];
        outDir = [params.storage.dir_results filesep params.storage.outName filesep 'coherency_channelGroups'];                
        fig_save(f, figname, outDir);
        close(f); 
    end
end

FB_y = cell(size(params.triggering.freqBands,1), N_groups*N_groups);
if plotEC_selGroups_freqBands

    %% === (2b) plotband connectivity in freq. bands: plot for each group ====
    vars2plot = {'DTF'};     % variables to plot
    for v = 1:size(vars2plot,2)
        disp(['connectivity analysis: plotting: ' vars2plot{v} ' ...']);
        M = eval(vars2plot{v});
        assert(size(M,5) == size(clzLabels,2));

        %% figure
        f = fig_make;
        nRows = size(params.triggering.freqBands,1);    % freq. bands
        nCols = numel(find(L == 1));                      % each group to each group
        nPlot = 1;

        % group-by-group connectivity -> C
        for freq = 1:size(params.triggering.freqBands,1)
            selFreq = params.triggering.freqBands{freq,2};
            i_fr = closestval(freqAxis,selFreq(1)):closestval(freqAxis,selFreq(2));
            g = 1;
            for ch_i = 1:size(selCh_groups,1)
                selCh_i = i_selCh_groups{ch_i};         % indices of channels in COH matrix
                groupName_i = params.connectivity.selectedChnls{ch_i,2};
                if ~isempty(selCh_i)
                    for ch_j = 1:size(selCh_groups,1)
                        selCh_j = i_selCh_groups{ch_j}; % indices of channels in COH matrix
                        groupName_j = params.connectivity.selectedChnls{ch_j,2};
                        if ~isempty(selCh_j)
                            % selected channels
                            nCh = size(M,1);
                            i_chMatrix = false(nCh,nCh);                % init
                            i_chMatrix(selCh_i,selCh_j) = true;         % selected channels (symmetric along diagonal, incl. diagonal)
                            i_subdiagonal = tril(ones(nCh,nCh),-1);     % below main diagonal
                            i_chMatrix = i_subdiagonal & i_chMatrix;    % deselect the diagonal
                            inds2use = repmat(i_chMatrix, [1,1,size(M,3),size(M,4),size(M,5)]);

                            % data
                            C = M;                                        % 5D: ch x ch x F x T x clzMean
                            C(~inds2use) = nan;                    % use chnl mask to select channels 

                            % data
                            y_avg = squeeze(nanmean(nanmean(nanmean(C(:,:,i_fr,:,:),1),2),3));     % 2D: T x clz
                            y_sem = squeeze(sem(nanmean(nanmean(C(:,:,i_fr,:,:),1),2),3));      % 2D: T x clz (SEM over freq. only!)
                            FB_y{freq,g} = y_avg;       % 2D: T x clz
                            g = g+1;

                            % data: old impl.
%                             y_avg = squeeze(nanmean(nanmean(nanmean(M(selCh_i,selCh_j,i_fr,:,:),1),2),3));     % 2D: T x clz
%                             y_sem = squeeze(sem(mean(mean(M(selCh_i,selCh_j,i_fr,:,:),1),2),3));      % 2D: T x clz (SEM over freq. only!)

                            % plot
                            subplot(nRows, nCols, nPlot);
                            hold on;
                            for clz = 1:size(clzLabels,2)
                                h_ax = plotband(timeAxis, y_avg(:,clz), y_sem(:,clz), trialsData_erp.info.colors(clz,:));
                            end
                            axis tight;
                            title({[vars2plot{v} ': ' groupName_i ', ' groupName_j]; ['freq = ' params.triggering.freqBands{freq,1}]});
                            xlabel('time (s)');
                            ylabel(vars2plot{v});
                            xticks = [-2:2];
                            for t = xticks
                                plot([t t], ylim, '--k');   % time ticks
                            end
                            nPlot = nPlot+1;
                        end
                    end
                end
            end
        end

        %% save figure
%         figname = [vars2plot{v} '_chGroups_freqBands'];
%         outDir = [params.storage.outputDir filesep 'coherency_channelGroups'];
        figname = [subjTag '_' vars2plot{v} '_chGroups_freqBands'];
        outDir = [params.storage.dir_results filesep params.storage.outName filesep 'coherency_channelGroups'];        
        fig_save(f, figname, outDir);
        close(f);    
    end

end    



%% plot connectivity: each channel to all channels (takes quite a lot of time: 100 chnls = 100 figures)
if plotEC_ch2ch_spectra   

    %% init struct spectra (for connectivity plotting)
    spectra = struct;
    spectra.freq = EC_info.freq;     % assumes row vector, 1 x N
    spectra.time = EC_info.time;
    spectra.labels = EC_info.clzLabels;
    spectra.clzNames = trialsData_erp.info.clzNames;
    spectra.clzNames{1,2} = trialsData_erp.info.clzNumber(1,1);
    spectra.clzNames{2,2} = trialsData_erp.info.clzNumber(2,1);
    spectra.info = trialsData_erp.info;
    COH_selCh = EC_info.selCh;
    
    %% (1a) connectivity a-la spectra: plot for each channel
    vars2plot = {'DTF'};     % variables to plot
    for v = 1:size(vars2plot,2)
        disp(['connectivity analysis: plotting: ' vars2plot{v} ' ...']);
        M = eval(vars2plot{v});
        for ch = 1:size(M,1)
           
            assert(size(M,5) == size(clzLabels,2));
            for clz = 1:length(clzLabels)           % plot for each clzMean
        
                % spectra -> struct spectralData
                D = squeeze(M(ch,:,:,:,:));        % 4D:  Ch x F x T x clzMean
                D = permute(D, [2, 3, 1, 4]);           % 4D:  F x T x Ch x clzMean
                spectra.data = D;
                params.selCh_H = COH_selCh;
                params.H = H;
                spectralData = prepare2plotSpectra(params, spectra, clz);
        
                % plot spectralData
                spectralData.info.text = [vars2plot{v} ', ' spectralData.info.text ', selCh = ' H.channels(COH_selCh(ch)).name];
                spectralData.info.figName = [vars2plot{v} '_ch' num2str(ch) '_clz' num2str(clz)];
                spectralData.info.outDir = [params.storage.outputDir filesep 'connectivity_' vars2plot{v}];
                plotSpectra(params, spectralData);    
            end
        end
    end
end

if plotEC_ch2ch_freqBands
    COH_selCh = EC_info.selCh;
    
    %% (1b) connectivity in freq. bands (FBs): plot for each channel
    vars2plot = {'DTF'};     % variables to plot
    for v = 1:size(vars2plot,2)
        disp(['connectivity analysis: plotting: ' vars2plot{v} ' ...']);
        M = eval(vars2plot{v});
        for ch = 1:size(M,1)
            for fb = 1:size(params.triggering.freqBands,1)
            
                % select frequencies
                freqBand = params.triggering.freqBands{fb,2};
                freqBandName = params.triggering.freqBands{fb,1};
                i_fr = closestval(EC_info.freq,freqBand(1)):closestval(EC_info.freq,freqBand(2));
            
                % re-define the trials struct
                D = squeeze(M(ch,:,:,:,:));             % 4D: Ch x F x T x clzMean
                D = squeeze(mean(D(:,i_fr,:,:),2));     % 3D: Ch x T x clzMean
                D = permute(D, [2, 1, 3]);              % 3D: T x Ch x clzMean
                trials_fb = struct;
                trials_fb.data = D;
                trials_fb.time = EC_info.time;
                assert(size(trials_fb.data,1) == size(trials_fb.time,1));
    %             trials.rejected = squeeze(spectra.rejected(1,:,:,:));   % assumes same rejection for all freq. bins of spectrograms
                trials_fb.labels = EC_info.clzLabels;
                trials_fb.clzNames = trialsData_erp.info.clzNames;
                trials_fb.clzNames{1,2} = trialsData_erp.info.clzNumber(1,1);
                trials_fb.clzNames{2,2} = trialsData_erp.info.clzNumber(2,1);                
                trials_fb.info = trialsData_erp.info;
    
                % arranges data from 'trials' to 'trialsData' structure
                trials_fb.selFreq = freqBand;
                trials_fb.freqBandName = freqBandName;
                trials_fb.figname = [vars2plot{v} '_ch' num2str(ch) '_' num2str(fb) '_' freqBandName '_' params.storage.subjTag];
                trials_fb.get_pVals = false;
                trials_fb.get_pValsBase = false;            
                params.plot_triggering.paraTimes = [];  % hack!
                trialsData = prepare2plotTrials(params, trials_fb);
            
                % plot: AVG +/- SEM
                trialsData.info.text = [vars2plot{v} ', selCh = ' H.channels(COH_selCh(ch)).name ', ' trialsData.info.text];
                trialsData.info.outDir = [params.storage.outputDir filesep vars2plot{v} '_' num2str(fb) '_' freqBandName];
                plotTrials(params, trialsData);  
            end     % of fb
        end     % of ch
    end     % of var2plot

end
