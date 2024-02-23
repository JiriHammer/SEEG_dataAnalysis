%function ch2roi_stackedFreqBands(params, groupInfo)
% groups channels into defined ROIs (anatomic areas), averages & plots
% FIG: 
%   - subplots = freqBands x ROIs
%   - 1 subplot = 1 class imagesc of stacked FB (y-axis) x time (x-axis) x activations (c-vals)
% based on the following criteria:
%   - selected anatomical areas (groups)
%   - selected freq. band
%   - significance (difference of conditions, from baseline (TO DO), or all)
% settings (subjects, directory names, ...) are stored in 'params'

% (c) Jiri, Sep21

disp('Grouping FB-channels activations ...');

%% settings: activations
time2plot = [-2,2];     % TO DO ...

% list_FB = cat(1, {'erp', [nan, nan]}, params.triggering.freqBands);   % add ERP to freq. bands
list_FB = params.triggering.freqBands;   % add ERP to freq. bands
freqAxis = 1:size(list_FB,1);

%% ============================IMPLEMENTATION==============================
disp('Grouping channels to ROIs (anatomical areas) ...');

%% load & plot frequency bands + ERP (over subjects and anatomy areas)
for grup = 1:size(groupInfo.anatomy_signif_list,1)
    list_anatomy_signif = groupInfo.anatomy_signif_list(grup,:);
    
    %% define ROIs (anatom. areas)
    anatAtlas = list_anatomy_signif{1,1};
    list_AA = anatomicalAreas_getList(params, anatAtlas);

    %% >>> get group data -> G <<<
    dataInfo.list_AA = list_AA;
    dataInfo.list_FB = list_FB;
    dataInfo.list_anatomy_signif = list_anatomy_signif;
%     [G, G_info] = ch2roi_load(params, dataInfo);
    GRUP = ch2roi_load(params, dataInfo);

    %% figure
    f = fig_make;
    nCols = size(list_AA,1);
    nRows = 2+4;
    nPlot = 1;

    % ---------FB activations for each class (E-I, I-E)---------------
    cLims = [-1, 1];
    for clz = 1:size(GRUP.trials{1,1},3)
        for aa = 1:size(list_AA,1)
            cVals = zeros(size(freqAxis,2), size(GRUP.trials{1,1},1));    % init, 2D: freq x time
            for freq = 1:size(GRUP.trials,2)
                if ~isempty(GRUP.trials{aa,freq})
                    S(1,:) = nanmean(GRUP.trials{aa,freq}(:,:,clz),2);    % 1D: 1 x t
                    cVals(freq,:) = S;                         % 2D: freq x time
                    clear S;
                end
            end

            % subplot
            subplot(nRows, nCols, nPlot);
            hold on;

            % >>> plot spectra <<<
            h = imagesc(GRUP.time{aa,freq}, freqAxis, cVals, cLims);
            colormap(gca,brewermap(256,'*RdBu'));     

            % axis setting
            axis tight;
            yLims = get(gca, 'ylim');
            xLims = get(gca, 'xlim');        
            plot([0 0], yLims, '--k');

            % title
            title([list_AA{aa} ', ' GRUP.info.clzNames{clz,1}]);
            %ylabel('freq [Hz]');
            if mod(nPlot,nCols) == 1
                set(gca,'YTick', freqAxis);
                set(gca, 'YTickLabel', list_FB(:,1));
            else
                set(gca,'YTick', []);
            end
            if nPlot >= (nRows-1)*nCols
                xlabel('time [s]');
            end

            nPlot = nPlot+1;
        end
    end    
    
    % ------------------ SNR ------------------
    cLims = [-1, 1];
    for aa = 1:size(list_AA,1)
        snrVals = zeros(size(freqAxis,2), size(GRUP.trials{1,1},1));    % init, 2D: freq x time
        for freq = 1:size(GRUP.trials,2)
            if ~isempty(GRUP.trials{aa,freq})
%                 S(1,:) = mean(GRUP.SNR_trials{aa,freq},2);    % mean over channels, 1D: 1 x t
                S(1,:) = var(nanmean(GRUP.trials{aa,freq},2),0,3)./mean(nanvar(GRUP.trials{aa,freq},0,2),3);    % SNR over channels, 1D: 1 x t
                snrVals(freq,:) = S;                         % 2D: freq x time
                clear S;
            end
        end
        
        % subplot
        subplot(nRows, nCols, nPlot);
        hold on;

        % >>> plot spectra <<<
        h = imagesc(GRUP.time{aa,freq}, freqAxis, snrVals, cLims);
        colormap(gca,brewermap(256,'PRGn'));     

        % axis setting
        axis tight;
        yLims = get(gca, 'ylim');
        xLims = get(gca, 'xlim');        
        plot([0 0], yLims, '--k');

        % title
        title(['SNR: ' list_AA{aa,1}]);
        if mod(nPlot,nCols) == 1
            set(gca,'YTick', freqAxis);
            set(gca, 'YTickLabel', list_FB(:,1));
        else
            set(gca,'YTick', []);
        end
        if nPlot >= (nRows-1)*nCols
            xlabel('time [s]');
        end
            
        nPlot = nPlot+1;
    end
    
    % ------------------ time: signif. differences ------------------
    nSigLags = [];
    for aa = 1:size(list_AA,1)
        for freq = 1:size(list_FB,1)
            nSigLags = cat(1, nSigLags, GRUP.sgnfLags{aa,freq}(:));
        end
    end
    cLims = [0, max([1e-6, prctile(nSigLags/100,95)])];     % 500 is a normalization factor = nSubj x nLags (TO DO!)
    for aa = 1:size(list_AA,1)
        
        % arrange valkues to time x freq "spectra"
        X = zeros(size(GRUP.trials,2), size(GRUP.trials{1,1},1));    % init, 2D: freq x time
        for freq = 1:size(GRUP.trials,2)
            if ~isempty(GRUP.trials{aa,freq})
                X(freq,:) = GRUP.sgnfLags{aa,freq}./size(GRUP.trials{aa,freq},2);    % 2D: freq x time, normalized by number of channels
            end
        end        
        
        % subplot
        subplot(nRows, nCols, nPlot);
        hold on;

        % >>> plot spectra <<<
        h = imagesc(GRUP.time{aa,freq}, freqAxis, X, cLims);
        colormap(gca,brewermap(256,'Purples'));     

        % axis setting
        axis tight;
        yLims = get(gca, 'ylim');
        xLims = get(gca, 'xlim');        
        plot([0 0], yLims, '--k');

        % title
        title('N signif. time lags');
        if mod(nPlot,nCols) == 1
            set(gca,'YTick', freqAxis);
            set(gca, 'YTickLabel', list_FB(:,1));
        else
            set(gca,'YTick', []);
        end
        if nPlot >= (nRows-1)*nCols
            xlabel('time [s]');
        end
            
        nPlot = nPlot+1;
    end
    
    % ------------------ nChnls in AAs ------------------
    for aa = 1:size(list_AA,1)
        
        % subplot
        subplot(nRows, nCols, nPlot);
        hold on;
        
        nChnls = zeros(1,size(GRUP.trials,2));
        for freq = 1:size(GRUP.trials,2)
%             nChnls = cat(2, nChnls, size(G{aa,freq},2)/AA_chnls(aa,freq)*100);
            if GRUP.nChnls_roi(aa,freq) ~= 0
                nChnls(freq) = size(GRUP.trials{aa,freq},2)/GRUP.nChnls_roi(aa,freq)*100;
            end
        end
        
%         histogram(nChnls, size(freqAxis,2), 'Orientation','horizontal');
        h = barh(freqAxis, nChnls);

        % axis setting
        axis tight;
        set(gca, 'xlim', [0, 100]);

        % title
        title(['signif. chnls (tot = ' num2str(GRUP.nChnls_roi(aa,freq)) ')']);
        if mod(nPlot,nCols) == 1
            set(gca,'YTick', freqAxis);
            set(gca, 'YTickLabel', list_FB(:,1));
        else
            set(gca,'YTick', []);
        end
%         if nPlot >= (nRows-1)*nCols
%             xlabel('N chnls (%)');
%         end
        xlabel('N chnls (%)');
        
        nPlot = nPlot+1;
    end
    
    % ------------------ nSubj ------------------
    for aa = 1:size(list_AA,1)
        
        % subplot
        subplot(nRows, nCols, nPlot);
        hold on;
        
        nSubj = [];
        for freq = 1:size(GRUP.trials,2)
            nSubj = cat(2, nSubj, length(unique(GRUP.nSubj_roi{aa,freq,1}))/length(unique(GRUP.nSubj_roi{aa,freq,2}))*100);
        end
        
%         histogram(nChnls, size(freqAxis,2), 'Orientation','horizontal');
        h = barh(freqAxis, nSubj);

        % axis setting
        axis tight;
        set(gca, 'xlim', [0, 100]);

        % title
        title(['signif. subj (tot = ' num2str(length(unique(GRUP.nSubj_roi{aa,freq,2}))) ')']);
        if mod(nPlot,nCols) == 1
            set(gca,'YTick', freqAxis);
            set(gca, 'YTickLabel', list_FB(:,1));
        else
            set(gca,'YTick', []);
        end
        if nPlot >= (nRows-1)*nCols
            xlabel('N subj. (%)');
        end
        
        nPlot = nPlot+1;
    end

    %% save fig
    outDir = [params.storage.dir_results filesep params.storage.outName filesep 'ch2roi' filesep 'stackedFreqBands'];
    figname = ['stackedFreq_' list_anatomy_signif{1} '_' list_anatomy_signif{3}];
    fig_save(f, figname, outDir);
    %close(f);    
        
end     % end of grup
