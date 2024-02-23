function ch2roi_plot_allFB_brain3D(GRUP, groupInfo)
% plots one ROI (brain3D + channels) with all activations
% FIG: 
%   - 3D brain plot with highlighted channels in MNI
%   - subplots = freqBands
%   - 1 subplot: time x clz
% based on the following criteria:
%   - selected anatomical areas (groups)
%   - selected freq. band
%   - significance (difference of conditions, from baseline (TO DO), or all)
% settings (subjects, directory names, ...) are stored in 'groupInfo'

% (c) Jiri, Aug21, Jul22

%% settings: activations
time2plot = groupInfo.time2plot;

%% required variables
list_FB = groupInfo.list_FB;
list_AA = groupInfo.list_AA;
aa = groupInfo.i_roi;       % selected anatomic area (ROI)
if ~isempty(GRUP.chnls_MNI_sel{aa,1})
    disp(['Grouping analysis, plotting: ' list_AA{aa,1} ' ...']);
else
    disp(['Grouping analysis, no channels found in: ' list_AA{aa,1} ' !']);
    return;
end

%% settings: brain 3D -> params
params_default;
params.plot_brain3D = struct(...        % user interface structure: holds most (but not all!) of the user settings, see also inside brain3D_plot.m
    'plot_brain3D', true, ...           % 3D brain model
    'volumes2plot', {{...               % volumes to plot (coregistered .nii files)
        'icbm152', [1 1 1], 0.07;  ...     % 1x3 cell: name, color, transparency (other choices for norm brain = wc1T1_colin27)
        }}, ...          
    'chVals_asCircles', true, ...       % plots channel values as color circles
    'chVals_asPatches', false, ...       % each channel value is mapped onto surface (but may be overwritten by latter points!), see also distanceThreshold for focality inside brain3D_plot.m
    'chVals_colorMap', [0 0 0; 0 0 0], ...    % colormap of channel values
    'chVals_patchTransp', 0.2, ...      % transparency of channel values
    'circleSizeLims', [5 5], ...       % in points, maps chnl values to circle size
    'size_interpolate', 1.0, ...        % in [mm], voxel size to which the brain is interpolated
    'backgroundColor', 'w', ...         % figure / axes background color, good choices: 'w', 'k'
    'text_LR', false, ...                % marks L / R side
    'text_tag', '', ...    % string
    'viewAngle', [0, 90], ...           % axial = 0,90; sagittal = 100,10; coronal front = 0,180; coronal back = 0,0
    'visible_axis', false ...           % if axis & their labels are visible
); 

%% add volume2plot: ROI
vol = 2;    % where vol = 1 -> whole brain 
nROIs = size(list_AA{aa,2},2);
if nROIs == 1
    clrs = [0 1 0];     % green
else
    clrs = colorPalette(3+nROIs);
    clrs(1:3,:) = [];
end
for roi = 1:nROIs
    roi2load = list_AA{aa,2}{1,roi};
%     if strcmp(roi2load, 'Visual central')
%         why;
%     end
    fileName = [params.path2others.normBrains filesep anatomicalAreas_atlas4AA(roi2load) filesep roi2load '.nii'];
    if exist(fileName,'file') == 2
        params.plot_brain3D.volumes2plot{vol,1} = roi2load;     % name  old: ['w' roi2load];   
        params.plot_brain3D.volumes2plot{vol,2} = clrs(roi,:);  % color
        params.plot_brain3D.volumes2plot{vol,3} = 0.3;          % transparency
        vol = vol+1;
    else
        disp([' - warning: ROI = ' roi2load ' not found in normalized files: ' params.path2others.normBrains]);
    end
end

%% camera views: AXIAL, SAG, COR
% sel_viewAngles =  [
%     0, 90;          % axial top = 0,90
%     0, -90;         % axial bottom ?
%     90, 0;          % sagittal left
%     -90, 0;          % sagittal right
%     180,0; ...      % coronal front = 180,0
%     0,0; ...        % coronal back = 0,0
%     ];
sel_viewAngles =  {
    'axial top', [0, 90];          % axial top = 0,90
    'sagittal R', [90, 5];          % sagittal left
    'sagittal L', [-90, 5];          % sagittal right
    };

%% get brain volumes
[brainVols, clrmap, inds_vol] = brainVolumes_fixed(params, params.plot_brain3D);
    
%% MNIs + VALues
chnls_MNI_VAL = cat(2, GRUP.chnls_MNI_all{aa,1}, ones(size(GRUP.chnls_MNI_all{aa,1},1),1));

%% figure
f = fig_make;
% nRows = size(list_FB,1);
% % nCols = size(list_AA,1);
% nCols = max(groupInfo.nROI);
nRows = max(6,size(list_FB,1));  % min = 6 for plotting
nCols = 2*size(sel_viewAngles,1)+1;                       % view angles (W = 2 cols) + FB (W = 1 col)
nPlot = 1;
% clrs = {'b','r','m','c','y','g','k'};
clrs = GRUP.info.clzColor;
marg_h = [0.1 0.1];
marg_w = [0.06 0.06];
gap = [0.05, 0.02];

%% plot figure: brain 3D
v_subplots = 1:nCols*nRows;
m_subplots = reshape(v_subplots, nCols, nRows)';
for v = 1:size(sel_viewAngles,1)
    % subplot
    i_subplots = [m_subplots(:,nPlot), m_subplots(:,nPlot+1)];
    ax = subtightplot(nRows, nCols, i_subplots(:), gap, marg_h, marg_w);
    hold on;
    
    % ---plot brain (no chnls values)---
    plotInfo = params.plot_brain3D;
    plotInfo.axHandle = ax;
    plotInfo.figHandle = f;
    plotInfo.chVals_lims = [0, 1];
    for vol = 1:size(brainVols,2)
        if brainVols{vol}.loaded
            plotInfo.thisVolume = vol;

            % plot all hnl values at once
            if strcmp(sel_viewAngles{v,1}, 'sagittal R')
                plotInfo.text_LR = false;
                i_ch = chnls_MNI_VAL(:,1) > 0;
                plot_chnls_MNI_VAL = chnls_MNI_VAL(i_ch,:);
                h = brain3D_plot(brainVols{vol}, plot_chnls_MNI_VAL, plotInfo); 
            elseif strcmp(sel_viewAngles{v,1}, 'sagittal L')
                plotInfo.text_LR = false;
                i_ch = chnls_MNI_VAL(:,1) <= 0;
                plot_chnls_MNI_VAL = chnls_MNI_VAL(i_ch,:);
                h = brain3D_plot(brainVols{vol}, plot_chnls_MNI_VAL, plotInfo);
            else
                plotInfo.text_LR = true;
                h = brain3D_plot(brainVols{vol}, chnls_MNI_VAL, plotInfo); 
            end
        end
    end    

    % set view angle
    view(sel_viewAngles{v,2});
    title(sel_viewAngles{v,1});
    nPlot = nPlot+2;        % jump by two columns
end
    
%% plot figure: plotband FB
i_subplots = [nCols:nCols:nCols*nRows];
for freq = 1:size(list_FB,1)
    freqBand = list_FB{freq,1};
    if ~isempty(GRUP.trials{aa,freq})
        if size(unique(GRUP.nSubj_roi{aa,freq,1}),1) >= groupInfo.minSubj

            % subplot
            subtightplot(nRows, nCols, i_subplots(freq), gap, marg_h, marg_w);
            hold on;

            % >>> plot FB activations <<<
            i_t = closestval(GRUP.time{aa,freq},time2plot(1)):closestval(GRUP.time{aa,freq},time2plot(2));
            for clz = 1:size(GRUP.trials{aa},3)
                plotband(GRUP.time{aa,freq}(i_t), mean(GRUP.trials{aa,freq}(i_t,:,clz),2), sem(GRUP.trials{aa,freq}(i_t,:,clz),2), clrs(clz,:));
            end

            % title
            str_chInfo = [freqBand ': N_c_h=' num2str(round(size(GRUP.trials{aa,freq},2)/GRUP.nChnls_roi(aa,freq)*100)) '%' , ... 
                        ', N_P=' num2str(size(unique(GRUP.nSubj_roi{aa,freq}),1))];
            if nPlot <= nCols
                title({list_AA{aa,1}; ['N_c_h=' num2str(GRUP.nChnls_roi(aa,freq))]; str_chInfo});
            else
                title(str_chInfo);
            end
            if nPlot > nCols*(nRows-1)
                xlabel('time [s]');
            end
            if mod(nPlot,nCols) == 1
                ylabel(freqBand);
            end
%                 plot([0 0], ylim, '--k');
            box on;
            axis tight;
            grid on;

            % plot grouped paradigm times (RT, ...)
            if isfield(GRUP, 'paraTimes')
                plot2axes_paraTimes(GRUP.paraTimes{aa,freq}, clrs);
            end
            
        end % of if nSubj > minSubj
    end % of if ~ismepty
    nPlot = nPlot+1;

end % of for freq = ...

% text on upper part of the figure
if isfield(GRUP.info, 'txtOnFig')
    tx = axes('visible','off', 'position',[0 0 1 1]);
    mytitle = ['ROI = ' list_AA{aa,1} ', ' GRUP.info.txtOnFig];
    mytitle = strrep(mytitle, '_','\_');
    text(0.016, 0.98, mytitle, 'fontsize', 14, 'fontw', 'bold');
end

%% save fig
% figname = [list_anatomy_signif{1} '_' list_anatomy_signif{2} '_' list_anatomy_signif{3}];
fig_save(f, groupInfo.figName, groupInfo.outDir, 'format','png');
close(f);    
