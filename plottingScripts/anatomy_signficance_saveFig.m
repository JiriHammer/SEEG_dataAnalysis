function anatomy_signficance_saveFig(params, plotInfo, groupInfo, chAss_isargAtlas)
% saves figure in plotInfo.fig to different folders based on significance and anatomical assignments
% output dir = slicesActivity_selected
% for example: ...\slicesActivity_selected\Yeo7_networks\Attention_hiGamma_not_sgnf ---

% (c) Jiri, Jun20

if ~isfield(groupInfo, 'anatomy_signif_list')
    return;
end

% go through all groups (rows) of the list
anatomy_signif_list = groupInfo.anatomy_signif_list;
for grp = 1:size(anatomy_signif_list,1)
    
    %% --- select frequency -> freqTag
    [tf, i_freq] = ismember(anatomy_signif_list{grp,2}, plotInfo.freqBands(:,1));
    assert(tf);
    freqTag = anatomy_signif_list{grp,2};
    
    %% --- significant? -> sgnfTag
    sgnfTag = [];
    if strcmp(anatomy_signif_list{grp,3},'was_sgnf') && (plotInfo.was_significant(i_freq) == 1)
        sgnfTag = 'was_sgnf';
    elseif strcmp(anatomy_signif_list{grp,3},'not_sgnf') && (plotInfo.was_significant(i_freq) == 0)
        sgnfTag = 'not_sgnf'; 
    elseif strcmp(anatomy_signif_list{grp,3},'any_sgnf')
        sgnfTag = 'any_sgnf';
    end
    
    %% --- anatomic area assignment -> anatTag  
    brainAtlas_name = anatomy_signif_list{grp,1};   
    thisAA = ch2roi_assignAA(brainAtlas_name, chAss_isargAtlas, groupInfo);
    
    if strcmp(thisAA, 'n.a.')
        anatTag = 'notAssigned';
    else
        anatTag = thisAA;
    end
    
    %% --- save figure -> slicesActivity_selected, e.g.: \Yeo7_networks\Attention_hiGamma_not_sgnf ---
    if ~isempty(sgnfTag)
        outDir =[params.storage.dir_results filesep params.storage.outName filesep ...
            'slicesActivity_selected' filesep ...
            anatomy_signif_list{grp,1} filesep anatTag '_' freqTag '_' sgnfTag];
        print_format = params.plot_brainSlices.printFormats;
        print_res = params.plot_brainSlices.printResolution;
        fig_save(plotInfo.fig, plotInfo.figname, outDir, 'res',print_res, 'format',print_format);
    end
end
