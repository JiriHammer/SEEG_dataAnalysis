%function poolSubjects_brain3D(params, groupInfo)
% plots 3D brain & activations (ERP + freq.band + spectra)
% pools together subject processed by a given job
% settings (subjects, directory names, ...) are stored in 'params'

% (c) Jiri, Apr18


disp('Plotting channels: 3D brain & activations ...');

%% camera views: AXIAL, SAG, COR
sel_viewAngles =  [
    0, 90;          % axial = 0,90
    90, 0;          % sagittal = 100,10
    180,0; ...      % coronal front = 180,0
    0,0; ...        % coronal back = 0,0
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
    list_AA = list_anatomy_signif{1,4};
    for aa = 1:size(list_AA,1)
        list_AA{aa,2} = {list_AA{aa,1}};        % quick fix ...
    end
%     list_AA = anatomicalAreas_getList(params, anatAtlas);   % TO DO: register DMN, DAN, ...
    
    %% >>> get group data -> G <<<
    dataInfo.list_AA = list_AA;
    dataInfo.list_FB = list_FB;
    dataInfo.list_anatomy_signif = list_anatomy_signif;
    GRUP = ch2roi_load(params, dataInfo);
    
    %% get brain volumes
    [brainVols, clrmap, inds_vol] = brainVolumes_fixed(params, params.plot_brain3D);
    
    %% ==== T-RESOLVED VALUES: values to plot & MNI coors ====
    for v = 1:size(groupInfo.vals2plot,1)
        % values to plot  
        val2plot = groupInfo.vals2plot{v,1};
        assert(isfield(GRUP, val2plot));
        V = GRUP.(val2plot);            % values: cell, V{aa} = 2D: t x ch
        assert(isfield(GRUP, 'chnls_MNI'));

        % selected times
        time2plot = groupInfo.time2plot;
        
        % ch-lims
        all_data = [];
        for aa = 1:size(list_AA,1)
            all_data = cat(1, all_data, V{aa}(:));  % cat all data
        end
        if strcmp(groupInfo.vals2plot{v,2}, 'prctile5')
            chVals_lims = [prctile(all_data,5),prctile(all_data,95)];
        elseif strcmp(groupInfo.vals2plot{v,2}, 'minmax')
            chVals_lims = [min(all_data),max(all_data)];
        else    % user defined
            assert(isnumeric(groupInfo.vals2plot{v,2}));
            assert(groupInfo.vals2plot{v,2}(1) < groupInfo.vals2plot{v,2}(2))
            chVals_lims = groupInfo.vals2plot{v,2};
        end

        % plot: clz & sel times
        for clz = 1:size(V{1},3)
            for t = 1:size(time2plot,2)

                %% plot figure
                f = fig_make;
                nRows = size(list_AA,1);   
                nCols = size(sel_viewAngles,1);  % AXIAL, SAG, COR views
                nPlot = 1;
                for aa = 1:size(list_AA,1)
                    i_t = closestval(GRUP.time{aa}, time2plot(t));
                    for col = 1:size(sel_viewAngles,1)
                        plotInfo = params.plot_brain3D;
                        
                        % axes
                        ax = subplot(nRows,nCols,nPlot); 
                        hold on;

                        % MNIs + VALues
                        chnls_MNI_VAL = cat(2, GRUP.chnls_MNI{aa,1}, V{aa,1}(i_t,:,clz)');
    %                     chVals_lims = [prctile(chnls_MNI_VAL(:,4),5),prctile(chnls_MNI_VAL(:,4),95)];

                        % re-define colormap
                        plotInfo.chVals_colorMap = groupInfo.vals2plot{v,3};
                        
                        % defined circle sizes
                        if strcmp(groupInfo.vals2plot{v,4}, 'circleSize_absScaled')
                            plotInfo.def_circle_size = linTransform(chnls_MNI_VAL(:,4), chVals_lims, [-plotInfo.circleSizeLims(2), plotInfo.circleSizeLims(2)]);
                            plotInfo.def_circle_size = abs(plotInfo.def_circle_size);   % small & high values -> large circles
                            plotInfo.def_circle_size(plotInfo.def_circle_size<plotInfo.circleSizeLims(1)) = plotInfo.circleSizeLims(1);   % clip close-to zero vals to circle size = 5
                        elseif strcmp(groupInfo.vals2plot{v,4}, 'circleSize_linScaled')
                            plotInfo.def_circle_size = linTransform(chnls_MNI_VAL(:,4), chVals_lims, [plotInfo.circleSizeLims(1), plotInfo.circleSizeLims(2)]);
                            plotInfo.def_circle_size(plotInfo.def_circle_size<plotInfo.circleSizeLims(1)) = plotInfo.circleSizeLims(1);   % clip close-to zero vals to circle size = 5
                        end                           
                            
                        % ---plot brain (no chnls values)---
                        plotInfo.axHandle = ax;
                        plotInfo.figHandle = f;
                        vol = 1;            % !!! assumes that first volume is the whole brain (grey matter)
                        if brainVols{vol}.loaded
                            plotInfo.thisVolume = vol;
                            h = brain3D_plot(brainVols{vol}, [], plotInfo);     % no channel values = []
                        end            

                        % ---plot brain volume (1 at a time)---
                        vol = aa+1;
                        plotInfo.chVals_lims = chVals_lims;     % computed automatically if left out
                        if brainVols{vol}.loaded
                            plotInfo.thisVolume = vol;
                            h = brain3D_plot(brainVols{vol}, chnls_MNI_VAL, plotInfo);
                        end    

                        % set view angle
                        view(sel_viewAngles(col,:));

                        nPlot = nPlot+1;
                    end
                end

                % text
                mytitle = [list_anatomy_signif{2} ', ' GRUP.info.clzNames{clz,1} ': T = ' num2str(time2plot(t)) ' s'];
                tx = axes('visible','off', 'position',[0 0 1 1]);
                mytitle = strrep(mytitle, '_','\_');
                text(0.016, 0.97, mytitle, 'fontsize', 14, 'fontw', 'bold');

                %% ---- save fig ----
                outDir = [params.storage.dir_results filesep params.storage.outName filesep 'ch2roi_brain3D' filesep val2plot];
                figname = ['tRes_signif_chnls_' list_anatomy_signif{2} '_clz' num2str(clz) '_t' num2str(t)];
                fig_save(f, figname, outDir, 'res',600, 'format','png');
                close(f);    

            end     % end of time2plot
        end     % end of clz
    end     % end of var2plot
    
    %% ==== ALL SIGNIF CHNLS: values to plot & MNI coors ===
    val2plot = 'chnls_hVals';   
    assert(isfield(GRUP, val2plot));
    V = GRUP.(val2plot);            % values
    assert(isfield(GRUP, 'chnls_MNI'));
    
    % plot figure
    f = fig_make;
    nRows = size(list_AA,1);   
    nCols = size(sel_viewAngles,1);  % AXIAL, SAG, COR views
    nPlot = 1;
    for aa = 1:size(list_AA,1)
        for col = 1:size(sel_viewAngles,1)
        
            % axes
            ax = subplot(nRows,nCols,nPlot); 
            hold on;

            % MNIs + VALues
            chnls_MNI_VAL = cat(2, GRUP.chnls_MNI{aa,1}, any(V{aa,1},1)');
    %         chVals_lims = [prctile(chnls_MNI_VAL(:,4),5),prctile(chnls_MNI_VAL(:,4),95)];
            chVals_lims = [min(chnls_MNI_VAL(:,4)),max(chnls_MNI_VAL(:,4))];

            % ---plot brain (no chnls values)---
            plotInfo = params.plot_brain3D;
            plotInfo.axHandle = ax;
            plotInfo.figHandle = f;
            vol = 1;            % !!! assumes that first volume is the whole brain (grey matter)
            if brainVols{vol}.loaded
                plotInfo.thisVolume = vol;
                h = brain3D_plot(brainVols{vol}, [], plotInfo);
            end            

            % ---plot brain volume (1 at a time)---
            vol = aa+1;
            plotInfo.chVals_lims = chVals_lims;     % computed automatically if left out
            if brainVols{vol}.loaded
                plotInfo.thisVolume = vol;
                h = brain3D_plot(brainVols{vol}, chnls_MNI_VAL, plotInfo);
            end    
            
            % set view angle
            view(sel_viewAngles(col,:));

            nPlot = nPlot+1;
        end
    end
    
    % ---- save fig ----
    outDir = [params.storage.dir_results filesep params.storage.outName filesep 'ch2roi_brain3D'];
    figname = ['signif_chnls_' list_anatomy_signif{2}];
    fig_save(f, figname, outDir, 'res',600, 'format','png');
    %close(f);    
    
    %% ==== T-RESOLVED SIGNIF CHNLS: values to plot & MNI coors ====
    % values to plot
    val2plot = 'chnls_hVals';   
    assert(isfield(GRUP, val2plot));
    V = GRUP.(val2plot);            % values: cell, V{aa} = 2D: t x ch
    assert(isfield(GRUP, 'chnls_MNI'));
    
    % selected times
    tStep = 0.1;
    selTime = [-3, 3];
    time2plot = selTime(1):tStep:selTime(2);
    
    for t = 1:size(time2plot,2)

        %% plot figure
        f = fig_make;
        nRows = size(list_AA,1);   
        nCols = size(sel_viewAngles,1);  % AXIAL, SAG, COR views
        nPlot = 1;
        for aa = 1:size(list_AA,1)
            i_t = closestval(GRUP.time{aa}, time2plot(t));
            for col = 1:size(sel_viewAngles,1)

                % axes
                ax = subplot(nRows,nCols,nPlot); 
                hold on;

                % MNIs + VALues
                chnls_MNI_VAL = cat(2, GRUP.chnls_MNI{aa,1}, V{aa,1}(i_t,:)');
        %         chVals_lims = [prctile(chnls_MNI_VAL(:,4),5),prctile(chnls_MNI_VAL(:,4),95)];
%                 chVals_lims = [min(chnls_MNI_VAL(:,4)),max(chnls_MNI_VAL(:,4))];
                chVals_lims = [0 1];        % hardcoded for hVals !!!

                % ---plot brain (no chnls values)---
                plotInfo = params.plot_brain3D;
                plotInfo.axHandle = ax;
                plotInfo.figHandle = f;
                vol = 1;            % !!! assumes that first volume is the whole brain (grey matter)
                if brainVols{vol}.loaded
                    plotInfo.thisVolume = vol;
                    h = brain3D_plot(brainVols{vol}, [], plotInfo);
                end            

                % ---plot brain volume (1 at a time)---
                vol = aa+1;
                plotInfo.chVals_lims = chVals_lims;     % computed automatically if left out
                if brainVols{vol}.loaded
                    plotInfo.thisVolume = vol;
                    h = brain3D_plot(brainVols{vol}, chnls_MNI_VAL, plotInfo);
                end    

                % set view angle
                view(sel_viewAngles(col,:));

                nPlot = nPlot+1;
            end
        end

        % text
        mytitle = ['T = ' num2str(time2plot(t)) ' s'];
        tx = axes('visible','off', 'position',[0 0 1 1]);
        mytitle = strrep(mytitle, '_','\_');
        text(0.016, 0.97, mytitle, 'fontsize', 14, 'fontw', 'bold');
        
        % ---- save fig ----
        outDir = [params.storage.dir_results filesep params.storage.outName filesep 'ch2roi_brain3D'];
        figname = ['tRes_signif_chnls_' list_anatomy_signif{2} '_t' num2str(t)];
        fig_save(f, figname, outDir, 'res',600, 'format','png');
        %close(f);    
        
    end     % end of time2plot
end



