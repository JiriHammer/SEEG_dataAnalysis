function ch2roi_plot_stackFreq(GRUP, groupInfo)
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

%% required variables
list_FB = groupInfo.list_FB;
list_AA = groupInfo.list_AA;
freqAxis = 1:size(list_FB,1);

%% used ROIs -> nClz & nSpectra (!!! taken from the significance criteria of last FB !!!)
sel_AA = [];    
for aa = 1:size(list_AA,1) 
    if ~isempty(GRUP.trials{aa,end})
        % must contain at least minSubj
        if size(unique(GRUP.nSubj_roi{aa,end,1}),1) >= groupInfo.minSubj  
            sel_AA = cat(2, sel_AA, aa); 
        end
    end
end
if isempty(sel_AA)
    disp('Not enough channels in ROIs found.');
    return;
end
nClz = size(GRUP.spectra_data{sel_AA(1)},3);

%% selected time
% time2plot = [-2,2];
time2plot = groupInfo.time2plot;

%% figure
f = fig_make;
nCols = size(sel_AA,2);  % size(list_AA,1);
nRows = 2+4;
nPlot = 1;
marg_h = [0.1 0.1];
marg_w = [0.04 0.04];
gap = [0.05, 0.02];

%% ---------FB activations for each class (E-I, I-E)---------------
% cLims = [-1, 1];
cLims = [-0.2, 0.2];
for clz = 1:nClz
    for aa = sel_AA % 1:size(list_AA,1)
        cVals = zeros(size(freqAxis,2), size(GRUP.time{aa,1},1));    % init, 2D: freq x time
        for freq = 1:size(GRUP.trials,2)
            if ~isempty(GRUP.trials{aa,freq})
                if size(unique(GRUP.nSubj_roi{aa,freq,1}),1) >= groupInfo.minSubj
                    S(1,:) = nanmean(GRUP.trials{aa,freq}(:,:,clz),2);    % 1D: 1 x t
                    cVals(freq,:) = S;                         % 2D: freq x time
                    clear S;
                end
            end
        end

        % subplot
        if nCols > 20
            subplot(nRows, nCols, nPlot);
        else
            subtightplot(nRows, nCols, nPlot, gap, marg_h, marg_w);
        end
        hold on;
        set(gca, 'ButtonDownFcn','call_copy');   % click on plot to see it bigger   
        
        % >>> plot spectra <<<
        i_t = closestval(GRUP.time{aa,freq},time2plot(1)):closestval(GRUP.time{aa,freq},time2plot(2));
        h = imagesc(GRUP.time{aa,freq}(i_t), freqAxis, cVals(:,i_t), cLims);
        colormap(gca,brewermap(256,'*RdBu'));     

        % axis setting
        axis tight;
        yLims = get(gca, 'ylim');
        xLims = get(gca, 'xlim');        
        plot([0 0], yLims, '--k');

        % title
        if clz == 1
            title({[list_AA{aa}], [GRUP.info.clzNames{clz,1}]});
        else
            title([GRUP.info.clzNames{clz,1}]);
        end
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

%% ------------------ SNR ------------------
cLims = [-2, 2];
for aa = sel_AA % 1:size(list_AA,1)
    snrVals = zeros(size(freqAxis,2), size(GRUP.time{aa,1},1));    % init, 2D: freq x time
    for freq = 1:size(GRUP.trials,2)
        if ~isempty(GRUP.trials{aa,freq})
            if size(unique(GRUP.nSubj_roi{aa,freq,1}),1) >= groupInfo.minSubj
%                 S(1,:) = mean(GRUP.SNR_trials{aa,freq},2);    % mean over channels, 1D: 1 x t
                S(1,:) = var(nanmean(GRUP.trials{aa,freq},2),0,3)./mean(nanvar(GRUP.trials{aa,freq},0,2),3);    % SNR over channels, 1D: 1 x t
                snrVals(freq,:) = S;                         % 2D: freq x time
                clear S;
            end
        end
    end

    % subplot
    if nCols > 20
        subplot(nRows, nCols, nPlot);
    else
        subtightplot(nRows, nCols, nPlot, gap, marg_h, marg_w);
    end
    hold on;
    set(gca, 'ButtonDownFcn','call_copy');   % click on plot to see it bigger   
    
    % >>> plot spectra <<<
    i_t = closestval(GRUP.time{aa,freq},time2plot(1)):closestval(GRUP.time{aa,freq},time2plot(2));
    h = imagesc(GRUP.time{aa,freq}(i_t), freqAxis, snrVals(:,i_t), cLims);
    colormap(gca,brewermap(256,'PRGn'));     

    % axis setting
    axis tight;
    yLims = get(gca, 'ylim');
    xLims = get(gca, 'xlim');        
    plot([0 0], yLims, '--k');

    % title
%     title(['SNR: ' list_AA{aa,1}]);
    title('SNR');
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

%% ------------------ time: signif. differences ------------------
nSigLags = [];
for aa = sel_AA % 1:size(list_AA,1)
    for freq = 1:size(list_FB,1)
        nSigLags = cat(1, nSigLags, GRUP.sgnfLags{aa,freq}(:));
    end
end
cLims = [0, max([1e-6, prctile(nSigLags/100,95)])];     % 500 is a normalization factor = nSubj x nLags (TO DO!)
for aa = sel_AA % 1:size(list_AA,1)

    % arrange values to time x freq "spectra"
    X = zeros(size(GRUP.trials,2), size(GRUP.time{aa,1},1));    % init, 2D: freq x time
    for freq = 1:size(GRUP.trials,2)
        if ~isempty(GRUP.trials{aa,freq})
            if size(unique(GRUP.nSubj_roi{aa,freq,1}),1) >= groupInfo.minSubj
                X(freq,:) = GRUP.sgnfLags{aa,freq}./size(GRUP.trials{aa,freq},2);    % 2D: freq x time, normalized by number of channels
            end
        end
    end        

    % subplot
    if nCols > 20
        subplot(nRows, nCols, nPlot);
    else
        subtightplot(nRows, nCols, nPlot, gap, marg_h, marg_w);
    end
    hold on;
    set(gca, 'ButtonDownFcn','call_copy');   % click on plot to see it bigger   
    
    % >>> plot spectra <<<
    i_t = closestval(GRUP.time{aa,freq},time2plot(1)):closestval(GRUP.time{aa,freq},time2plot(2));
    h = imagesc(GRUP.time{aa,freq}(i_t), freqAxis, X(:,i_t), cLims);
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

%% ------------------ nChnls in AAs ------------------
for aa = sel_AA % 1:size(list_AA,1)

    % subplot
    if nCols > 20
        subplot(nRows, nCols, nPlot);
    else
        subtightplot(nRows, nCols, nPlot, gap, marg_h, marg_w);
    end
    hold on;
    set(gca, 'ButtonDownFcn','call_copy');   % click on plot to see it bigger   
    
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

%% ------------------ nSubj ------------------
for aa = sel_AA % 1:size(list_AA,1)

    % subplot
    if nCols > 20
        subplot(nRows, nCols, nPlot);
    else
        subtightplot(nRows, nCols, nPlot, gap, marg_h, marg_w);
    end
    hold on;
    set(gca, 'ButtonDownFcn','call_copy');   % click on plot to see it bigger   
    
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

%% --------------- text on upper part of the figure ----------------
if isfield(GRUP.info, 'txtOnFig')
    tx = axes('visible','off', 'position',[0 0 1 1]);
    mytitle = GRUP.info.txtOnFig;
    mytitle = strrep(mytitle, '_','\_');
    text(0.016, 0.98, mytitle, 'fontsize', 14, 'fontw', 'bold');
end

%% save fig
fig_save(f, groupInfo.figName, groupInfo.outDir);
close(f);    
