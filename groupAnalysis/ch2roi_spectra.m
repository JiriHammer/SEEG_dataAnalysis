%function ch2roi_selFreq(params, groupInfo)
% groups channels into defined ROIs (anatomic areas), averages & plots
% based on the following criteria:
%   - selected anatomical areas (groups)
%   - selected freq. band
%   - significance (difference of conditions, from baseline (TO DO), or all)
% settings (subjects, directory names, ...) are stored in 'params'

% (c) Jiri, Aug21

disp('Grouping channels spectra ...');

%% selected time & freq
sel_freq = [0, 120];
sel_time = [-3, 3;];

%% load & plot frequency bands + ERP (over subjects and anatomy areas)
for grup = 1:size(groupInfo.anatomy_signif_list,1)
    list_anatomy_signif = groupInfo.anatomy_signif_list(grup,:);
    
    %% define freq band
    freqBand = list_anatomy_signif{1,2};
    list_FB{1} = freqBand;

    %% define ROIs (anatom. areas)
    anatAtlas = list_anatomy_signif{1,1};
    list_AA = anatomicalAreas_getList(params, anatAtlas);
    
    %% >>> get group data -> G <<<
    dataInfo.list_AA = list_AA;
    dataInfo.list_FB = list_FB;
    dataInfo.list_anatomy_signif = list_anatomy_signif;
    GRUP = ch2roi_load(params, dataInfo);
    
    %% selected indices: time & freq
    i_fr = [closestval(GRUP.spectra_freq, sel_freq(1)):closestval(GRUP.spectra_freq, sel_freq(2))];
    i_t = [closestval(GRUP.spectra_time, sel_time(1)):closestval(GRUP.spectra_time, sel_time(2))];
    
    %% plot figure
    aa_found = [];    
    for aa = 1:size(list_AA,1) 
        if ~isempty(GRUP.spectra_data{aa}), aa_found = cat(2, aa_found, aa); end
    end
    nClz = size(GRUP.spectra_data{aa_found(1)},3);
    
    f = fig_make;
    nRows = nClz + 1;   % nClz + SNR
    nCols = size(aa_found,2);    % AAs
    nPlot = 1;
    fontSize = 10;
    
    % --------------- spectra of ROIs ------------------
%     cLims = [-1, 1];
    cLims = [-0.1, 0.1];    % z-score
    for aa = 1:size(list_AA,1)
        if ~isempty(GRUP.spectra_data{aa})
            for clz = 1:nClz
                
                % axes
                subplot(nRows, nCols, nPlot+(clz-1)*nCols);
                hold on;
                
                % plot spectra
                imagesc(GRUP.spectra_time(i_t), GRUP.spectra_freq(i_fr), GRUP.spectra_data{aa}(i_fr,i_t,clz), cLims);
                colormap(gca,brewermap(256,'*RdBu'));
                axis tight;
                
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
        end
    end
    
    %% save
    outDir = [params.storage.dir_results filesep params.storage.outName filesep 'ch2roi' filesep 'spectra'];
    figname = [list_anatomy_signif{1} '_' list_anatomy_signif{3}];
    fig_save(f, figname, outDir, 'res',600);
    %close(f);    
    
end     % end of grup
