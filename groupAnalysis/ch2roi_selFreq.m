%function ch2roi_selFreq(params, groupInfo)
% groups channels into defined ROIs (anatomic areas), averages & plots
% based on the following criteria:
%   - selected anatomical areas (groups)
%   - selected freq. band
%   - significance (difference of conditions, from baseline (TO DO), or all)
% settings (subjects, directory names, ...) are stored in 'params'

% (c) Jiri, Aug21

disp('Grouping channels activations ...');

%% old: selected anatomical areas
% AA_list = anatomicalAreas_getList(params, 'neurologists_grouping');
% AA_selected = {'Inferior frontal gyrus', 'Middle frontal gyrus', ...
%     'Middle temporal gyrus', 'Superior temporal gyrus', 'Fusiform gyrus', 'Medial occipital cortex',...
%     'Precentral gyrus','Postcentral gyrus'};

%% load & plot frequency bands + ERP (over subjects and anatomy areas)
for grup = 1:size(groupInfo.anatomy_signif_list,1)
    list_anatomy_signif = groupInfo.anatomy_signif_list(grup,:);
    
    %% define freq band
    freqBand = list_anatomy_signif{1,2};
    list_FB{1} = freqBand;

    %% define ROIs (anatom. areas)
    anatAtlas = list_anatomy_signif{1,1};
    list_AA = anatomicalAreas_getList(params, anatAtlas);
    
    %% >>> get group data -> GRUP <<<
    dataInfo.list_AA = list_AA;
    dataInfo.list_FB = list_FB;
    dataInfo.list_anatomy_signif = list_anatomy_signif;
    GRUP = ch2roi_load(params, dataInfo);
    
    %% plot figure
    f = fig_make;
    [nRows, nCols] = getSubplotLayout(size(list_AA,1));
%     nRows = 2;
%     nCols = 4;
    nPlot = 1;
    clrs = GRUP.info.clzColor;
    for aa = 1:size(list_AA,1)
        if ~isempty(GRUP.trials{aa})
            subplot(nRows, nCols, nPlot);
            hold on;
            for clz = 1:size(GRUP.trials{aa},3)
                plotband(GRUP.time{aa}, mean(GRUP.trials{aa}(:,:,clz),2), sem(GRUP.trials{aa}(:,:,clz),2), clrs(clz,:));
            end
            title([list_AA{aa,1}, ', N_c_h = ' num2str(size(GRUP.trials{aa},2)), ', N_P = ' num2str(size(unique(GRUP.nSubj_roi{aa}),1))]);
            xlabel('time [s]');
            ylabel(freqBand);
            axis tight;
            grid on;
            box on;
        end
        nPlot = nPlot+1;
    end
            
    %% save
    outDir = [params.storage.dir_results filesep params.storage.outName filesep 'ch2roi' filesep 'selFreq'];
    figname = [list_anatomy_signif{1} '_' list_anatomy_signif{2} '_' list_anatomy_signif{3}];
    fig_save(f, figname, outDir);
    %close(f);    
    
end     % end of grup
