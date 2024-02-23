% function poolSubjects_brain3D_activity(params, groupInfo)
% plots 3D brain & activations (ERP + freq.band + spectra)
% pools together subject processed by a given job
% settings (subjects, directory names, ...) are stored in 'params'

% (c) Jiri, Apr18


disp('Plotting channels: 3D brain & activations ...');

%% camera views: AXIAL, SAG, COR
sel_viewAngles =  [
    0, 90;          % axial = 0,90
    90, 0;          % sagittal = 100,10
    -90, 0;          % sagittal = 100,10
...    180,0; ...      % coronal front = 180,0
...    0,0; ...        % coronal back = 0,0
    ];
% params.plot_brain3D.backgroundColor = 'k';
params.plot_brain3D.backgroundColor = 'w';

%% load & plot frequency bands + ERP (over subjects and anatomy areas)
for grup = 1:size(groupInfo.anatomy_signif_list,1)
    list_anatomy_signif = groupInfo.anatomy_signif_list(grup,:);
    
    %% define freq band
    freqBand = list_anatomy_signif{1,2};
    list_FB{1} = freqBand;

    %% selected ROIs (anatom. areas)
%     list_AA = list_anatomy_signif{1,4};
%     for aa = 1:size(list_AA,1)
%         list_AA{aa,2} = {list_AA{aa,1}};        % quick fix ...
%     end
    list_AA = anatomicalAreas_getList(params, list_anatomy_signif{1,4}); 
    
    %% >>> get group data -> G <<<
    dataInfo.list_AA = list_AA;
    dataInfo.list_FB = list_FB;
    dataInfo.list_anatomy_signif = list_anatomy_signif;
    GRUP = ch2roi_load(params, dataInfo);
    
    %% get brain volumes
    [brainVols, clrmap, inds_vol] = brainVolumes_fixed(params, params.plot_brain3D);
    
    %% values to plot: trials 
    val2plot = 'trials';   
    assert(isfield(GRUP, val2plot));
    V = GRUP.(val2plot);            % values: cell, V{aa} = 3D: t x ch x clz
    assert(isfield(GRUP, 'chnls_MNI'));
    
    % ch-lims from all data
    all_data = [];
    for aa = 1:size(list_AA,1)
        all_data = cat(1, all_data, V{aa}(:));
    end
    chVals_lims_V = [prctile(all_data,5),prctile(all_data,95)];

    %% values to plot: SNR
    val2plot = 'SNR_trials';   
    assert(isfield(GRUP, val2plot));
    S = GRUP.(val2plot);            % values: cell, V{aa} = 2D: t x ch
    
    % ch-lims from all data
    all_data = [];
    for aa = 1:size(list_AA,1)
        all_data = cat(1, all_data, S{aa}(:));
    end
    chVals_lims_S = [prctile(all_data,5),prctile(all_data,95)];
    
    %% ----------- plot for each time: sel vars (all clz) -------------- 
    time2plot = groupInfo.time2plot;    % selected times
    for t = 1:size(time2plot,2)
        i_t = closestval(GRUP.time{aa}, time2plot(t));
        
        %% === plot figure ===
        f = fig_make;
        nRows = size(V{1},3)*ismember('trials', groupInfo.vals2plot(:,1)) + ...
                size(S{1},3)*ismember('SNR_trials', groupInfo.vals2plot(:,1));            % = nClz   
        nCols = size(sel_viewAngles,1);  % = views: AXIAL, SAG, COR views
        nPlot = 1;
        plotInfo = params.plot_brain3D;
        plotInfo.chVals_asPatches = false;
        plotInfo.chVals_asCircles = true;
        plotInfo.chVals_patchTransp = 0.8;
        
        % ----------- trials -----------
        if ismember('trials', groupInfo.vals2plot(:,1))
            for clz = 1:size(V{1},3)

                % cat all channel values (only for first aa): MNIs + VALues & circle size
                chnls_MNI_VAL = [];
                chnls_SUBJ = [];
                for roi = 1:size(list_AA,1)
                    chnls_MNI_VAL = cat(1, chnls_MNI_VAL, [GRUP.chnls_MNI{roi,1}, V{roi,1}(i_t,:,clz)']);
                    chnls_SUBJ = cat(1, chnls_SUBJ, GRUP.nSubj_roi{roi,1,2});
                    assert(size(chnls_MNI_VAL,1) == size(chnls_SUBJ,1));
                end

                % average over different channels in a larger voxel
                d_mm = 5;       % larger voxel size, in [mm]
                n_subjContribute = 1;   % number of different subjects in each larger voxel
                [chVals_out, chMNIs_out] = voxelAvg_chVals_inCube(chnls_MNI_VAL(:,4), chnls_MNI_VAL(:,1:3), chnls_SUBJ, d_mm, n_subjContribute);
                chnls_MNI_VAL = cat(2, chMNIs_out', chVals_out);

                % define circle size based on chnl values
                plotInfo.def_circle_size = brain3D_getCircleSize(chnls_MNI_VAL, chVals_lims_V, plotInfo.circleSizeLims);
%                 plotInfo.def_circle_size = linTransform(chnls_MNI_VAL(:,4), chVals_lims_V, [-plotInfo.circleSizeLims(2), plotInfo.circleSizeLims(2)]);
%                 plotInfo.def_circle_size = abs(plotInfo.def_circle_size);   % small & high values -> large circles
%                 plotInfo.def_circle_size(plotInfo.def_circle_size<plotInfo.circleSizeLims(1)) = plotInfo.circleSizeLims(1);   % clip close-to zero vals to circle size = 5                      

                % plot different view angles
                for col = 1:size(sel_viewAngles,1)

                    % axes
                    ax = subplot(nRows,nCols,nPlot); 
                    hold on;
                    plotInfo.axHandle = ax;
                    plotInfo.figHandle = f;
                    plotInfo.text_tag = GRUP.info.clzNames{clz,1};

                    % ---plot brain (with chnls values)---
                    vol = 1;            % !!! assumes that first volume is the whole brain (grey matter)
                    if brainVols{vol}.loaded
                        plotInfo.thisVolume = vol;
                        h = brain3D_plot(brainVols{vol}, chnls_MNI_VAL, plotInfo);     % no channel values = []
                    end            

                    % ---plot brain volume (1 at a time)---
    %                 for vol = 2:size(plotInfo.volumes2plot,1)
    %                     if brainVols{vol}.loaded
    %                         plotInfo.thisVolume = vol;
    %                         h = brain3D_plot(brainVols{vol}, [], plotInfo);
    %                     end    
    %                 end

                    % set view angle
                    view(sel_viewAngles(col,:));

                    nPlot = nPlot+1;
                end
            end
        end

        % ----------- SNR ----------- (copy pasted from above trials, TO DO: loop in for cycle over sel. vars)
        if ismember('SNR_trials', groupInfo.vals2plot(:,1))
            % cat all channel values (only for first aa): MNIs + VALues & circle size
            chnls_MNI_VAL = [];
            chnls_SUBJ = [];
            for roi = 1:size(list_AA,1)
                chnls_MNI_VAL = cat(1, chnls_MNI_VAL, [GRUP.chnls_MNI{roi,1}, S{roi,1}(i_t,:)']);
                chnls_SUBJ = cat(1, chnls_SUBJ, GRUP.nSubj_roi{roi,1,2});
                assert(size(chnls_MNI_VAL,1) == size(chnls_SUBJ,1));
            end

            % average over different channels in a larger voxel
            d_mm = 5;       % larger voxel size, in [mm]
            n_subjContribute = 1;   % number of different subjects in each larger voxel
            [chVals_out, chMNIs_out] = voxelAvg_chVals_inCube(chnls_MNI_VAL(:,4), chnls_MNI_VAL(:,1:3), chnls_SUBJ, d_mm, n_subjContribute);
            chnls_MNI_VAL = cat(2, chMNIs_out', chVals_out);

            % define circle size based on chnl values
            plotInfo.def_circle_size = linTransform(chnls_MNI_VAL(:,4), chVals_lims_S, plotInfo.circleSizeLims);
    %         plotInfo.def_circle_size = abs(plotInfo.def_circle_size);   % small & high values -> large circles
            plotInfo.def_circle_size(plotInfo.def_circle_size<plotInfo.circleSizeLims(1)) = plotInfo.circleSizeLims(1);   % clip close-to zero vals to circle size = 5                      
            plotInfo.chVals_colorMap = brewermap(256,'RdPu');

            % plot different view angles
            for col = 1:size(sel_viewAngles,1)

                % axes
                ax = subplot(nRows,nCols,nPlot); 
                hold on;
                plotInfo.axHandle = ax;
                plotInfo.figHandle = f;

                % ---plot brain (with chnls values)---
                vol = 1;            % !!! assumes that first volume is the whole brain (grey matter)
                if brainVols{vol}.loaded
                    plotInfo.thisVolume = vol;
                    h = brain3D_plot(brainVols{vol}, chnls_MNI_VAL, plotInfo);     % no channel values = []
                end            

                % ---plot brain volume (1 at a time)---
    %             for vol = 2:size(plotInfo.volumes2plot,1)
    %                 if brainVols{vol}.loaded
    %                     plotInfo.thisVolume = vol;
    %                     h = brain3D_plot(brainVols{vol}, [], plotInfo);
    %                 end    
    %             end

                % set view angle
                view(sel_viewAngles(col,:));

                nPlot = nPlot+1;
            end
        end
        
        % ----- text -----
        mytitle = [list_anatomy_signif{2} ', T = ' num2str(time2plot(t)) ' s'];
        tx = axes('visible','off', 'position',[0 0 1 1]);
        mytitle = strrep(mytitle, '_','\_');
        text(0.016, 0.97, mytitle, 'fontsize', 14, 'fontw', 'bold');

        %% ---- save fig ----
        outDir = [params.storage.dir_results filesep params.storage.outName filesep 'ch2roi_brain3D' filesep 'trials_SNR'];
        figname = [list_anatomy_signif{2} '_t' num2str(t)];
        fig_save(f, figname, outDir, 'res',600, 'format','png');
        close(f);    

    end     % end of time2plot

    
end



