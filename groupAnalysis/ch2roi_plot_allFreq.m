function ch2roi_plot_allFreq(GRUP, groupInfo)
% groups channels into defined ROIs (anatomic areas), averages & plots
% FIG: 
%   - subplots = freqBands x ROIs
%   - 1 subplot: time x clz
% based on the following criteria:
%   - selected anatomical areas (groups)
%   - selected freq. band
%   - significance (difference of conditions, from baseline (TO DO), or all)
% settings (subjects, directory names, ...) are stored in 'params'

% (c) Jiri, Aug21, Jul22

%% settings: default
if ~isfield(groupInfo, 'subtractAvgResponse')
    groupInfo.subtractAvgResponse = false;
end

%% settings: activations
time2plot = groupInfo.time2plot;

%% required variables
list_FB = groupInfo.list_FB;
list_AA = groupInfo.list_AA;

%% figure
f = fig_make;
nRows = size(list_FB,1);
nCols = size(list_AA,1);
% nCols = max(groupInfo.nROI);
nPlot = 1;
clrs = GRUP.info.clzColor;
marg_h = [0.1 0.1];
marg_w = [0.04 0.04];
gap = [0.05, 0.02];

% plot
for freq = 1:size(list_FB,1)
    freqBand = list_FB{freq,1};
    for aa = 1:size(list_AA,1)
        
        % axes
        if nCols > 20
            subplot(nRows, nCols, nPlot);
        else
            subtightplot(nRows, nCols, nPlot, gap, marg_h, marg_w);
        end
        hold on;      
        box on;
        set(gca, 'ButtonDownFcn','call_copy');   % click on plot to see it bigger 
        
        % title, labels
        str_chInfo = ['C=' num2str(round(size(GRUP.trials{aa,freq},2)/GRUP.nChnls_roi(aa,freq)*100)) '%' , ... 
                    ', P=' num2str(size(unique(GRUP.nSubj_roi{aa,freq}),1))];
        if nPlot <= nCols
            if nCols <= 15  % long AA string
                title({list_AA{aa,1}; ['C=' num2str(GRUP.nChnls_roi(aa,freq))]; str_chInfo});
            else            % short AA string (to avoid overlap)
                title({shortenChannelName(list_AA{aa,1},7); ['P=' num2str(GRUP.nChnls_roi(aa,freq))]; str_chInfo});
            end
        else
            title(str_chInfo);
        end
        if nPlot > nCols*(nRows-1)
            xlabel('time [s]');
        end
        if mod(nPlot,nCols) == 1
            ylabel(freqBand);
        end
                
        % >>> plot (if nSubj_roi >= minSubj) <<<
        if ~isempty(GRUP.trials{aa,freq})
            if size(unique(GRUP.nSubj_roi{aa,freq,1}),1) >= groupInfo.minSubj
               
                % data to plot
                i_t = closestval(GRUP.time{aa,freq},time2plot(1)):closestval(GRUP.time{aa,freq},time2plot(2));
                x = GRUP.time{aa,freq}(i_t);
                y_avg = mean(GRUP.trials{aa,freq}(i_t,:,:),2);
                y_sem = sem(GRUP.trials{aa,freq}(i_t,:,:),2);
                if groupInfo.subtractAvgResponse
                    y_avg = y_avg - mean(y_avg,3);  % subtract AVG (non-specific) response (mean over clz)
                end
                
                % >>> plot grouped channels activity in GRUP.trials{aa,freq} <<<
                for clz = 1:size(GRUP.trials{aa,freq},3)
                    plotband(x, y_avg(:,:,clz), y_sem(:,:,clz), clrs(clz,:));
                end
                axis tight;
                
                % plot the AVG (non-specific) response (mean over clz)
                if ~groupInfo.subtractAvgResponse
                    y_nonSpec = mean(y_avg,3);  % mean over clz
                    plot(x, y_nonSpec, ':k', 'LineWidth',2);
                end
            
                % plot baseline values
                if isfield(GRUP, 'BASE_vals')
                    if ~isempty(GRUP.BASE_vals{aa,freq})
                        B_avg = GRUP.BASE_vals{aa,freq}(1);
                        B_sem = GRUP.BASE_vals{aa,freq}(2);
    %                     plotband(xlim, [B_avg,B_avg], [B_sem,B_sem], 'g');
                        plot(xlim, [B_avg,B_avg], 'g', 'LineWidth',2);
                    end
                end
            
                % plot significance of grouped activations
                if isfield(GRUP, 'H_trials')
                    plot2axes_signif_filled(GRUP.time{aa,freq}(i_t), GRUP.H_trials{aa,freq}(i_t), GRUP.P_trials{aa,freq}(i_t));
                end

                % plot grouped paradigm times (RT, ...)
                if isfield(GRUP, 'paraTimes')
                    plot2axes_paraTimes(GRUP.paraTimes{aa,freq}, clrs);
                end

                box on;
                grid on;
                grid minor;
            end % of if nSubj > minSubj
        end % of if ~ismepty
        nPlot = nPlot+1;
    end % for aa = ...
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
