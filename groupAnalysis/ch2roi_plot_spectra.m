function ch2roi_plot_spectra(GRUP, groupInfo)
% groups channels into defined ROIs (anatomic areas), averages & plots
% based on the following criteria:
%   - selected anatomical areas (groups)
%   - selected freq. band
%   - significance (difference of conditions, from baseline (TO DO), or all)
% settings (subjects, directory names, ...) are stored in 'params'

% (c) Jiri, Aug21

%% required variables
list_AA = groupInfo.list_AA;

%% selected time & freq
sel_freq = [0, 120];
sel_time = groupInfo.time2plot;
i_fr = [closestval(GRUP.spectra_freq, sel_freq(1)):closestval(GRUP.spectra_freq, sel_freq(2))];
i_t = [closestval(GRUP.spectra_time, sel_time(1)):closestval(GRUP.spectra_time, sel_time(2))];

%% used ROIs -> nClz & nSpectra
% aa_found = [];    
% for aa = 1:size(list_AA,1) 
%     if ~isempty(GRUP.spectra_data{aa}), aa_found = cat(2, aa_found, aa); end
% end
% nClz = size(GRUP.spectra_data{aa_found(1)},3);
% % nSpectra = size(aa_found,2);  % TO DO
% nSpectra = aa_found(end);   % plots with gaps

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

%% get color limits
cVals = [];
cVals_snr = [];
for aa = sel_AA
    if ~isempty(GRUP.spectra_data{aa})
        if size(unique(GRUP.nSubj_roi{aa,end,1}),1) >= groupInfo.minSubj  
            % spectral vals
            S = GRUP.spectra_data{aa}(i_fr,i_t,:);
            cVals = cat(1, cVals, S(:));
            
            % SNR vals
            S_snr = GRUP.spectra_snr{aa}(i_fr,i_t);
            cVals_snr = cat(1, cVals_snr, S_snr(:));            
        end
    end
end
prc_S = [prctile(cVals,5), prctile(cVals,95)];
maxAbs_S = max(abs(prc_S));
cLims_S = [-maxAbs_S, maxAbs_S];

prc_S = [prctile(cVals_snr,5), prctile(cVals_snr,95)];
maxAbs_S = max(abs(prc_S));
cLims_snr = [-maxAbs_S, maxAbs_S];

%% plot figure: max 20 cols per page
for aa = sel_AA
    if ~isempty(GRUP.spectra_data{aa})
        if size(unique(GRUP.nSubj_roi{aa,end,1}),1) >= groupInfo.minSubj  
            
            %% figure
            f = fig_make;
            nRows = 1;
            nCols = nClz + 1;   % nClz + SNR
            nPlot = 1;
            fontSize = 10;

            %% --------------- spectra of ROI for each clz ------------------
            % cLims = [-1, 1];
            cLims = 2*cLims_S;    % z-score            
            for clz = 1:nClz

                % axes
                subplot(nRows, nCols, nPlot);
                hold on;

                % plot spectra
                imagesc(GRUP.spectra_time(i_t), GRUP.spectra_freq(i_fr), GRUP.spectra_data{aa}(i_fr,i_t,clz), cLims);
                colormap(gca,brewermap(256,'*RdBu'));
                axis tight;

                % plot grouped paradigm times (RT, ...) - based on last FB
                if isfield(GRUP, 'paraTimes')
                    plot2axes_paraTimes(GRUP.paraTimes{aa,end}(clz,:,:), GRUP.info.clzColor(clz,:), 0);
                end  
            
                % ticks
                plot([0 0], ylim, '--k');
                freq_ticks = [8, 12, 30, 50, 100];
                for t = freq_ticks
                    plot(xlim, [t, t], '--k');
                end

                % labels
                if clz == 1
                    title([strrep(list_AA{aa,1}, '_', '\_') ...
                           ', N_c_h = ' num2str(size(GRUP.trials{aa,end},2)) ...
                           ', N_P = ' num2str(size(unique(GRUP.nSubj_roi{aa,end,1}),1)) ', '...
                           GRUP.info.clzNames{clz}], 'FontSize',fontSize);
                else
                    title(GRUP.info.clzNames{clz}, 'FontSize',fontSize);
                end
                xlabel('time (s)', 'FontSize',fontSize);
                if mod(nPlot,nCols) == 1
                    ylabel('freq (Hz)', 'FontSize',fontSize);
                else
                    set(gca, 'YTick', []);
                end
                nPlot = nPlot+1;
            end
            
            %% --------------- SNR over chnls in the ROI ------------------
%             cLims = [-1, 1];
            cLims = 2*cLims_snr;
            % axes
            subplot(nRows, nCols, nPlot);
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

            %% save
            figName = [groupInfo.figName '_' list_AA{aa,1}];
            fig_save(f, figName, groupInfo.outDir, 'res',600);
            close(f);    
            
        end % of if nSubj > minSubj
    end % of if ~isempty
end % for aa


