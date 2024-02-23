function ch2roi_plot_selFreq(GRUP, groupInfo)
% groups channels into defined ROIs (anatomic areas), averages & plots
% based on the following criteria:
%   - selected anatomical areas (groups)
%   - selected freq. band
%   - significance (difference of conditions, from baseline (TO DO), or all)
% settings (subjects, directory names, ...) are stored in 'params'

% (c) Jiri, Aug21, Jul22
    
%% required variables
freq = groupInfo.i_freq;
freqBand = groupInfo.freqBand;
list_AA = groupInfo.list_AA;
time2plot = groupInfo.time2plot;

%% figure
f = fig_make;
[nRows, nCols] = getSubplotLayout(groupInfo.nROI(freq));
nPlot = 1;
clrs = GRUP.info.clzColor;
marg_h = [0.1 0.1];
marg_w = [0.04 0.04];
gap = [0.08, 0.03];
fontSize = 12;
% gap = [0.10, 0.04];
% fontSize = 16;

% plot ROIs
for aa = 1:size(list_AA,1)
    if ~isempty(GRUP.trials{aa,freq})
        if size(unique(GRUP.nSubj_roi{aa,freq,1}),1) >= groupInfo.minSubj
            
            % axes
%             subplot(nRows, nCols, nPlot);
            subtightplot(nRows, nCols, nPlot, gap, marg_h, marg_w);
            hold on;
            set(gca, 'ButtonDownFcn','call_copy');   % click on plot to see it bigger   
            set(gca, 'FontSize',fontSize);
            
            % data to plot
            i_t = closestval(GRUP.time{aa,freq},time2plot(1)):closestval(GRUP.time{aa,freq},time2plot(2));
            x = GRUP.time{aa,freq}(i_t);
            y_avg = mean(GRUP.trials{aa,freq}(i_t,:,:),2);  % mean over channels
            y_sem = sem(GRUP.trials{aa,freq}(i_t,:,:),2);   % sem  over channels
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
            
            % plot RT
            if isfield(GRUP, 'paraTimes')
                plot2axes_paraTimes(GRUP.paraTimes{aa,freq}, clrs);
            end
            
            % title & labels
            title({strrep([list_AA{aa,1}], '_', '\_'), ['N_c_h = ' num2str(size(GRUP.trials{aa,freq},2)), ', N_P = ' num2str(size(unique(GRUP.nSubj_roi{aa,freq,1}),1))]}, 'FontSize',fontSize);
            if nPlot > nCols*(nRows-1)
                xlabel('time [s]', 'FontSize',fontSize);
            end
            if mod(nPlot,nCols) == 1
                ylabel(freqBand, 'FontSize',fontSize);
            end
            grid on;
            grid minor;
            box on;
            
            nPlot = nPlot+1;
        end
    end
end

% text on upper part of the figure
if isfield(GRUP.info, 'txtOnFig')
    tx = axes('visible','off', 'position',[0 0 1 1]);
    mytitle = ['FB = ' freqBand ', ' GRUP.info.txtOnFig];
    mytitle = strrep(mytitle, '_','\_');
    text(0.016, 0.98, mytitle, 'fontsize', 14, 'fontw', 'bold');
end

%% save
fig_save(f, groupInfo.figName, groupInfo.outDir);
close(f);    

