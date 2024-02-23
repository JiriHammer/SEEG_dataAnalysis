function ch2roi_plot_AUC_allFB_allROI(GRUP, groupInfo)
% plots AUC (area under the curve), or sum in certain time window, for
%    - clz 1 vs. baseline
%    - clz 2 vs. baseline
%    - clz 1 vs. clz 2

% (c) Jiri, Jun23

%% settings: AUC
if ~isfield(groupInfo, 'AUC_settings')
    groupInfo.AUC_settings.baselineWindow = [-0.5 0];
    groupInfo.AUC_settings.timeWindow = [0 0.5];
end
t_bsln = groupInfo.AUC_settings.baselineWindow;
t_wind = groupInfo.AUC_settings.timeWindow;

%% required variables
list_FB = groupInfo.list_FB;
list_AA = groupInfo.list_AA;

%% figure
f = fig_make;
nRows = size(list_FB,1);    % rows = FBs
nCols = size(groupInfo.clzNames,1)+1;   % cols = clz vs. baseline + clz1 vs. clz2
% nCols = size(list_AA,1);
% nCols = max(groupInfo.nROI);
nPlot = 1;
clrs = GRUP.info.clzColor;
marg_h = [0.1 0.1];
marg_w = [0.04 0.04];
gap = [0.05, 0.02];
clrs = colorPalette(size(list_AA,1));

% plot
for freq = 1:size(list_FB,1)
    freqBand = list_FB{freq,1};
    
    % --- plot: clz vs. baseline ---
    for clz = 1:size(groupInfo.clzNames,1)
        subtightplot(nRows, nCols, nPlot, gap, marg_h, marg_w);
        hold on;

        % pool sum of activations over AAs
        for aa = 1:size(list_AA,1)
            % if nSubj_roi >= minSubj
            if ~isempty(GRUP.trials{aa,freq})
                if size(unique(GRUP.nSubj_roi{aa,freq,1}),1) >= groupInfo.minSubj
                    % data: baseline
                    i_t = closestval(GRUP.time{aa,freq},t_bsln(1)):closestval(GRUP.time{aa,freq},t_bsln(2));                
                    s2 = sum(GRUP.trials{aa,freq}(i_t,:,clz),1)./diff(t_bsln);
                    x = aa-0.2;
                    h = boxplot(s2, 'Positions',x, 'Symbol','w', 'Widths',0.1, 'Colors',clrs(aa,:)); 

                    % data: time window of clz
                    i_t = closestval(GRUP.time{aa,freq},t_wind(1)):closestval(GRUP.time{aa,freq},t_wind(2));                
                    s2 = sum(GRUP.trials{aa,freq}(i_t,:,clz),1)./diff(t_wind);
                    x = aa+0.2;
                    h = boxplot(s2, 'Positions',x, 'Symbol','w', 'Widths',0.1, 'Colors',clrs(aa,:));                 
                end
            end
        end
        axis tight;
        xlim([0, size(list_AA,1)+1]);
        set(gca, 'XTick', [1:size(list_AA,1)]);
        set(gca, 'XTickLabel', list_AA(:,1));
        set(gca,'XTickLabelRotation',45);
        grid on;
        ylabel(['AUC: ' freqBand]);
        title(['baseline vs. ' groupInfo.clzNames{clz,1}]);
        nPlot = nPlot+1;
    end
    
    % --- plot: clz1 vs. clz2 vs. clz3 ... ---
    subtightplot(nRows, nCols, nPlot, gap, marg_h, marg_w);
    hold on;
    
    % pool sum of activations over AAs
    for aa = 1:size(list_AA,1)
        % if nSubj_roi >= minSubj
        if ~isempty(GRUP.trials{aa,freq})
            if size(unique(GRUP.nSubj_roi{aa,freq,1}),1) >= groupInfo.minSubj
                nClz = size(groupInfo.clzNames,1);
                for clz = 1:size(groupInfo.clzNames,1)
                    % data: time window of clz
                    i_t = closestval(GRUP.time{aa,freq},t_wind(1)):closestval(GRUP.time{aa,freq},t_wind(2));                
                    s2 = sum(GRUP.trials{aa,freq}(i_t,:,clz),1)./diff(t_wind);
                    x = aa+0.2*(clz-nClz/2);
                    h = boxplot(s2, 'Positions',x, 'Symbol','w', 'Widths',0.1/(nClz/2), 'Colors',clrs(aa,:));  
                end
            end
        end
    end
    axis tight;
    xlim([0, size(list_AA,1)+1]);
    set(gca, 'XTick', [1:size(list_AA,1)]);
    set(gca, 'XTickLabel', list_AA(:,1));
    set(gca,'XTickLabelRotation',45);
    grid on;
    ylabel(['AUC: ' freqBand]);
    title_str = [groupInfo.clzNames{1,1}];
    for clz = 2:size(groupInfo.clzNames,1)
        title_str = [title_str ' vs. ' groupInfo.clzNames{clz,1}];
    end
    title(title_str);
    nPlot = nPlot+1;    
    
end % of for freq = ...

% text on upper part of the figure
if isfield(GRUP.info, 'txtOnFig')
    tx = axes('visible','off', 'position',[0 0 1 1]);
    mytitle = GRUP.info.txtOnFig;
    mytitle = strrep(mytitle, '_','\_');
    text(0.016, 0.98, mytitle, 'fontsize', 14, 'fontw', 'bold');
end

%% save fig
% figname = [list_anatomy_signif{1} '_' list_anatomy_signif{2} '_' list_anatomy_signif{3}];
fig_save(f, groupInfo.figName, groupInfo.outDir);
close(f);    
