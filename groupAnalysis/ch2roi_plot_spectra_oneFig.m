function ch2roi_plot_spectra_oneFig(GRUP, groupInfo)
% plots spectra + SNR of grouped channels into defined ROIs (anatomic areas)
% based on the following criteria:
%   - selected anatomical areas (groups)
%   - selected freq. band
%   - significance (difference of conditions, from baseline (TO DO), or all)
% settings (subjects, directory names, ...) are stored in 'params'

% (c) Jiri, Aug21, Aug22

%% required variables
list_AA = groupInfo.list_AA;

%% selected time & freq
sel_freq = [0, 120];
sel_time = groupInfo.time2plot;
    
%% selected indices: time & freq
i_fr = [closestval(GRUP.spectra_freq, sel_freq(1)):closestval(GRUP.spectra_freq, sel_freq(2))];
i_t = [closestval(GRUP.spectra_time, sel_time(1)):closestval(GRUP.spectra_time, sel_time(2))];
    
%% used ROIs -> nClz & nSpectra
sel_AA = [];    
for aa = 1:size(list_AA,1) 
    if ~isempty(GRUP.spectra_data{aa})
        % must contain at least minSubj, !!! taken from the significance criteria of last FB !!!
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

%% plot figure
% aa_found = [];    
% for aa = 1:size(list_AA,1) 
%     if ~isempty(GRUP.spectra_data{aa}), aa_found = cat(2, aa_found, aa); end
% end
% nClz = size(GRUP.spectra_data{aa_found(1)},3);

f = fig_make;
nRows = nClz + 1;   % nClz + SNR
nCols = size(sel_AA,2);    % AAs
nPlot = 1;
fontSize = 10;
marg_h = [0.1 0.13];
marg_w = [0.04 0.04];
gap = [0.07, 0.02];

% --------------- spectra of ROIs ------------------
%     cLims = [-1, 1];
cLims = [-0.1, 0.1];    % z-score
for aa = sel_AA   % 1:size(list_AA,1)
    if ~isempty(GRUP.spectra_data{aa})
        for clz = 1:nClz

            % axes
%             subplot(nRows, nCols, nPlot+(clz-1)*nCols);
            subtightplot(nRows, nCols, nPlot+(clz-1)*nCols, gap, marg_h, marg_w);
            hold on;

            % plot spectra
            imagesc(GRUP.spectra_time(i_t), GRUP.spectra_freq(i_fr), GRUP.spectra_data{aa}(i_fr,i_t,clz), cLims);
            colormap(gca,brewermap(256,'*RdBu'));
            axis tight;

            % plot grouped paradigm times (RT, ...) - based on last FB
            if isfield(GRUP, 'paraTimes')
                if size(GRUP.paraTimes{aa,end},1) > 1
                    plot2axes_paraTimes(GRUP.paraTimes{aa,end}(clz,:,:), GRUP.info.clzColor(clz,:), 0);
                end
            end            
            
            % ticks
            plot([0 0], ylim, '--k');
            freq_ticks = [8, 12, 30, 50, 100];
            for t = freq_ticks
                plot(xlim, [t, t], '--k');
            end

            % labels
            if clz == 1
                title({list_AA{aa,1}; ...
                       ['N_c_h = ' num2str(size(GRUP.trials{aa},2))]; ['N_P = ' num2str(size(unique(GRUP.nSubj_roi{aa}),1))]; ...
                       GRUP.info.clzNames{clz}}, 'FontSize',fontSize);
            else
                title(GRUP.info.clzNames{clz}, 'FontSize',fontSize);
            end
            xlabel('time (s)', 'FontSize',fontSize);
            if mod(nPlot,nCols) == 1
                ylabel('freq (Hz)', 'FontSize',fontSize);
            else
                set(gca, 'YTick', []);
            end
        end
        nPlot = nPlot+1;
    end
end

% --------------- SNR over chnls in ROIs ------------------
nPlot = nCols * nClz + 1;
cLims = [-1, 1];
%     cLims = [-0.1, 0.1];    % z-score
for aa = 1:size(list_AA,1)
    if ~isempty(GRUP.spectra_snr{aa})
        % axes
%         subplot(nRows, nCols, nPlot);
        subtightplot(nRows, nCols, nPlot, gap, marg_h, marg_w);
        hold on;

        % plot spectra
        imagesc(GRUP.spectra_time(i_t), GRUP.spectra_freq(i_fr), GRUP.spectra_snr{aa}(i_fr,i_t), cLims);
        colormap(gca,brewermap(256,'PRGn'));
        axis tight;
            
        % ticks
        plot([0 0], ylim, '--k');
        freq_ticks = [8, 12, 30, 50, 100];
        for t = freq_ticks
            plot(xlim, [t, t], '--k');
        end

        % labels
        title('SNR', 'FontSize',fontSize);
        xlabel('time [s]', 'FontSize',fontSize);
        if mod(nPlot,nCols) == 1
            ylabel('freq (Hz)', 'FontSize',fontSize);
        else
            set(gca, 'YTick', []);
        end
        nPlot = nPlot+1;
    end
end

% --------------- text on upper part of the figure ----------------
if isfield(GRUP.info, 'txtOnFig')
    tx = axes('visible','off', 'position',[0 0 1 1]);
    mytitle = GRUP.info.txtOnFig;
    mytitle = strrep(mytitle, '_','\_');
    text(0.016, 0.98, mytitle, 'fontsize', 14, 'fontw', 'bold');
end

%% save
figName = groupInfo.figName;
fig_save(f, figName, groupInfo.outDir, 'res',600);
close(f);    

