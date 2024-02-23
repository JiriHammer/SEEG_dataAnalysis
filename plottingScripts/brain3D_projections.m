function brain3D_projections(params, vals, plotInfo)
% plots iEEG channel activations on 3D transparent brain model
% prints axial, sagittal and coronal views

% (c) Jiri, Apr17

%% color maps
clrmap.brain = gray(128);                           % colormap for brain (T1, T2, CT, ...)
if isfield(params.plot_brainTopo, 'colorMap')
    clrmap.chnls = params.plot_brainTopo.colorMap;
else
    % clrmap.chnls = jet(128);                            % default colormap for values: jet
    % clrmap.chnls = getColorMap('bwr', 128);             % colormap for values: blue - white - red
    clrmap.chnls = getColorMap('bcwwmr', 128);          % colormap for values: blue - cyan - white - magenta - red
end
clrmap.fig = cat(1, clrmap.brain, clrmap.chnls);    % colormap of the figure
alphaVal = 0.1;                                     % transparency of the colored values (0 = opaque)

%% bkg/text colors
if ~isfield(plotInfo, 'fig_bkgClr')
    plotInfo.fig_bkgClr = 'k';
end
if strcmp(plotInfo.fig_bkgClr, 'k')
    plotInfo.txtClr = 'w';
else
    plotInfo.txtClr = 'k';
end
if ~isfield(plotInfo, 'visible_axis')
    plotInfo.visible_axis = true;
end

%% output directory
if isfield(plotInfo, 'outDir')
    outDir = [plotInfo.outDir  filesep 'projections_3D'];
else
    outDir = [params.storage.outputDir filesep 'projections_3D'];
end
if ~exist(outDir, 'dir')
    mkdir(outDir);
end  

%% figure name
if isfield(plotInfo, 'figName')
    figname = plotInfo.figName;
else
    figname = 'notNamed';
end

%% pass info from loaded brain (see getBrainData.m)
VI = plotInfo.brain.VI;      % interpolated volume
xi = plotInfo.brain.xi;      % interpolated x-axis, in [mm] of MNI coors
yi = plotInfo.brain.yi;      % interpolated y-axis, in [mm] of MNI coors
zi = plotInfo.brain.zi;      % interpolated z-axis, in [mm] of MNI coors
assert(size(plotInfo.chnlsMNI,2) == size(vals,1));
V = linTransform(VI, [min(VI(:)), max(VI(:))], [0, 1]);

%% isosurface
if ~isfield(plotInfo.brain, 'fv')
    separationThreshold = 0.5;                  % (value 0.5 separates gray matter from the dark background). Other values 0 - 1 may work also fine
    plotInfo.brain.fv = isosurface(V, separationThreshold);    % surface, vertex 
end

%% figure
% --- figure
f = figure('units','normalized','outerposition',[0 0 1 1]);
% f = figure('visible','on', 'Color',plotInfo.fig_bkgClr);
% set(f, 'Position', [100 100 round(1.0*1120) round(1.0*840)]);
% set(f, 'Position', [1 -479 2880 1463]);
set(f, 'Color', 'k');
set(f, 'InvertHardcopy','off');                 % preserves black background
set(f, 'Colormap', clrmap.fig);
set(gcf, 'Renderer','OpenGL');
opengl('hardware');

nRows = 4;
nCols = 7;

%% channel color limits
if ~isfield(plotInfo, 'chnl_clims')
    clims = [min(vals(:)),max(vals(:))];
else
    clims = plotInfo.chnl_clims;
end
%circle_size = 40;

%% colorbar
plotInfo.clims = clims;
plotInfo.inds = size(clrmap.brain,1)+[1,size(clrmap.chnls,1)];
clrBar_axes(plotInfo);

%% 3D model data -> plotInfo
plotInfo.vals = vals;
plotInfo.clrmap = clrmap;
plotInfo.alphaVal = alphaVal;
plotInfo.circle_size = 32;
angle_camLight = 10;

%%  >>>>>>>>>>>>>>>> axial (top) view: 3D brain model <<<<<<<<<<<<<<<<<<<<<
subplot(nRows, nCols, [1,2,3,8,9,10,15,16,17,22,23,24]);
hold on;
plotInfo.circle_size = 60;
%brain3D_model(plotInfo);
plot_brain3D_patch(plotInfo);
view(0,90);
camlight(0+angle_camLight,90+angle_camLight);
title('axial (top) view', 'Color',plotInfo.txtClr, 'fontsize',14, 'fontw','bold');
drawnow;

%%  >>>>>>>>>>>>>>>> coronal (front) view: 3D brain model <<<<<<<<<<<<<<<<<<<<<
subplot(nRows, nCols, [4,5,11,12]);
hold on;
plotInfo.circle_size = 40;
%brain3D_model(plotInfo);
plot_brain3D_patch(plotInfo);
view(180,0);
camlight(180+angle_camLight,0+angle_camLight);
title('coronal (front) view', 'Color',plotInfo.txtClr, 'fontsize',14, 'fontw','bold');
drawnow;

%%  >>>>>>>>>>>>>>>> sagittal (R-side) view: 3D brain model <<<<<<<<<<<<<<<<<<<<<
subplot(nRows, nCols, [6,7,13,14]);
hold on;
%brain3D_model(plotInfo);
plot_brain3D_patch(plotInfo);
view(100,10);
camlight(100+angle_camLight,10+angle_camLight);
title('sagittal (R-side) view', 'Color',plotInfo.txtClr, 'fontsize',14, 'fontw','bold');
drawnow;

%%  >>>>>>>>>>>>>>>> coronal (back) view: 3D brain model <<<<<<<<<<<<<<<<<<<<<
subplot(nRows, nCols, [18,19,25,26]);
hold on;
%brain3D_model(plotInfo);
plot_brain3D_patch(plotInfo);
view(0,0);
camlight(0+angle_camLight,0+angle_camLight);
title('coronal (back) view', 'Color',plotInfo.txtClr, 'fontsize',14, 'fontw','bold');
drawnow;

%%  >>>>>>>>>>>>>>>> sagittal (L-side) view: 3D brain model <<<<<<<<<<<<<<<<<<<<<
subplot(nRows, nCols, [20,21,27,28]);
hold on;
%brain3D_model(plotInfo);
plot_brain3D_patch(plotInfo);
view(80+180,10);
camlight(80+180+angle_camLight,10+angle_camLight);
title('sagittal (L-side) view', 'Color',plotInfo.txtClr, 'fontsize',14, 'fontw','bold');
drawnow;

%% text
if isfield(plotInfo, 'text')
    tx = axes('visible','off', 'position',[0 0 1 1]);
    mytitle = plotInfo.text;
    mytitle = strrep(mytitle, '_','\_');
    text(0.016, 0.97, mytitle, 'fontsize', 14, 'fontw', 'bold', 'Color',plotInfo.txtClr);
end

%% save
disp('Saving figure ...');
set(f, 'InvertHardcopy','off');                 % preserves black background
set(f, 'PaperPositionMode','auto');
%saveas(f, [outDir filesep figname '.fig']);
if params.plot_brainTopo.printResolution == 0
    print(f, '-dpng','-r0', [outDir filesep figname '.png']);
else
    print(f, '-dpng','-r600', [outDir filesep figname '.png']);
end
close(f);    
disp(['Figure: ' figname ' stored in: ' outDir]);


