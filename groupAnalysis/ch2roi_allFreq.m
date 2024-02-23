%function ch2roi_allFreq(params, groupInfo)
% TBD
% groups channels into defined ROIs (anatomic areas), averages & plots
% FIG: 
%   - subplots = freqBands x ROIs
%   - 1 subplot: time x clz
% based on the following criteria:
%   - selected anatomical areas (groups)
%   - selected freq. band
%   - significance (difference of conditions, from baseline (TO DO), or all)
% settings (subjects, directory names, ...) are stored in 'params'

% (c) Jiri, Aug21

disp('Grouping channels activations ...');

%% settings: activations
time2plot = [-2,2];
% allFreqBands = cat(1, {'erp', [nan, nan]}, params.triggering.freqBands);   % add ERP to freq. bands
list_FB = cat(1, {'erp', [nan, nan]}, params.triggering.freqBands);   % add ERP to freq. bands

%% load & plot frequency bands + ERP (over subjects and anatomy areas)
for grup = 1:size(groupInfo.anatomy_signif_list,1)
    list_anatomy_signif = groupInfo.anatomy_signif_list(grup,:);
    
    %% define ROIs (anatom. areas)
    anatAtlas = list_anatomy_signif{1,1};
    list_AA = anatomicalAreas_getList(params, anatAtlas);
%     AA_selected = AA_list(:,1);    

    %% >>> get group data -> G <<<
    dataInfo.list_AA = list_AA;
    dataInfo.list_FB = list_FB;
    dataInfo.list_anatomy_signif = list_anatomy_signif;
    GRUP = ch2roi_load(params, dataInfo);
    
    %% figure
    f = fig_make;
    nRows = size(list_FB,1);
    nCols = size(list_AA,1);
    nPlot = 1;
    clrs = {'b','r','m','c'};
    
    %% plot figure
    for freq = 1:size(list_FB,1)
        freqBand = list_FB{freq,1};
        for aa = 1:size(list_AA,1)
            if ~isempty(GRUP.trials{aa})
                subplot(nRows, nCols, nPlot);
                hold on;
                i_t = closestval(GRUP.time{aa,freq},time2plot(1)):closestval(GRUP.time{aa,freq},time2plot(2));
                for clz = 1:size(GRUP.trials{aa},3)
                    plotband(GRUP.time{aa,freq}(i_t), mean(GRUP.trials{aa,freq}(i_t,:,clz),2), sem(GRUP.trials{aa,freq}(i_t,:,clz),2), clrs{clz});
                end
                
                % title
                str_chInfo = ['N_c_h = ' num2str(round(size(GRUP.trials{aa,freq},2)/GRUP.nChnls_roi(aa,freq)*100)) '%' , ... 
                            ', N_P = ' num2str(size(unique(GRUP.nSubj_roi{aa,freq}),1))];
                if nPlot <= nCols
                    title({[list_AA{aa,1} ', N_c_h = ' num2str(GRUP.nChnls_roi(aa,freq))]; str_chInfo});
                else
                    title(str_chInfo);
                end
                if nPlot >= nCols*(nRows-1)
                    xlabel('time [s]');
                end
                if mod(nPlot,nCols) == 1
                    ylabel(freqBand);
                end
%                 plot([0 0], ylim, '--k');
                box on;
                axis tight;
                grid on;
                yLims = get(gca, 'ylim');
                if yLims(1) > -0.5 && yLims(2) < 0.5
                    set(gca, 'ylim', [-0.5, 0.5]);
                end
            end
            nPlot = nPlot+1;
        end
    end

    % save fig
    outDir = [params.storage.dir_results filesep params.storage.outName filesep 'ch2roi' filesep 'allFreq'];
    figname = [list_anatomy_signif{1} '_' list_anatomy_signif{2} '_' list_anatomy_signif{3}];
    fig_save(f, figname, outDir);
    %close(f);    
        
end     % end of grup



%% old... (TO DO: revive?) selected anatomical areas
% AA_list = anatomicalAreas_getList(params, 'neurologists_grouping');
% AA_selected = {'Inferior frontal gyrus', 'Middle frontal gyrus', ...
%     'Middle temporal gyrus', 'Superior temporal gyrus', 'Fusiform gyrus', 'Medial occipital cortex',...
%     'Precentral gyrus','Postcentral gyrus'};
