function plot_brain3D_wrapper(params, vals, plotInfo)
% plots iEEG channel activations on 3D transparent brain model

% (c) Jiri, Jan17

%% color maps
clrmap.brain = gray(128);                           % colormap for brain (T1, T2, CT, ...)
if isfield(params.plot_brainTopo, 'colorMap')
    clrmap.chnls = params.plot_brainTopo.colorMap;
else
    clrmap.chnls = jet(128);                            % default colormap for values: jet
    % clrmap.chnls = getColorMap('bwr', 128);             % colormap for values: blue - white - red
    % clrmap.chnls = getColorMap('bcwwmr', 128);          % colormap for values: blue - cyan - white - magenta - red
end
clrmap.fig = cat(1, clrmap.brain, clrmap.chnls);    % colormap of the figure
alphaVal = 0.1;                                     % transparency of the colored values (0 = opaque)

%% bkg/text colors
if ~isfield(plotInfo, 'fig_bkgClr')
    plotInfo.fig_bkgClr = 'w';
end
if strcmp(plotInfo.fig_bkgClr, 'w')
    plotInfo.txtClr = 'k';
else
    plotInfo.txtClr = 'w';
end
if ~isfield(plotInfo, 'visible_axis')
    plotInfo.visible_axis = true;
end

%% output directory
if isfield(plotInfo, 'outDir')
    outDir = [plotInfo.outDir  filesep 'model_3D'];
else
    outDir = [params.storage.outputDir filesep 'model_3D'];
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
if ~isfield(plotInfo, 'figHandle')
    f = figure('units','normalized','outerposition',[0 0 1 1]);
    %set(f, 'Position', [100 100 round(1.0*1120) round(1.0*840)]);
    %set(f, 'Position', [1 -479 2880 1463]);
    set(f, 'visible','on', 'Color',plotInfo.fig_bkgClr);
    set(f, 'InvertHardcopy','off');                 % preserves black background
%     set(f, 'Colormap', clrmap.fig);
    set(gcf, 'Renderer','OpenGL');    
else
    f = plotInfo.figHandle;
end
%set(f, 'Colormap', clrmap.fig);
f.Colormap = clrmap.fig;
opengl('hardware');

%% text
if isfield(plotInfo, 'text')
    tx = axes('visible','off', 'position',[0 0 1 1]);
    mytitle = plotInfo.text;
    mytitle = strrep(mytitle, '_','\_');
    text(0.016, 0.97, mytitle, 'fontsize', 14, 'fontw', 'bold', 'Color',plotInfo.txtClr);
end

%% channel color limits
if ~isfield(plotInfo, 'chnl_clims')
    clims = [min(vals(:)),max(vals(:))];
else
    clims = plotInfo.chnl_clims;
end

%% colorbar
plotInfo.clims = clims;
plotInfo.inds = size(clrmap.brain,1)+[1,size(clrmap.chnls,1)];
%clrBar_axes(plotInfo);

%% axes
if ~isfield(plotInfo, 'axHandle')
    model_axes = axes('visible','on', 'position',[0.05 0.05 0.9 0.9]);  % position
else
    model_axes = plotInfo.axHandle;
end
axes(model_axes);
hold on;

%% >>>>>>>>>>>>> plot: 3D brain model + color channel values <<<<<<<<<<<<<<
plotInfo.vals = vals;
plotInfo.clrmap = clrmap;
plotInfo.alphaVal = alphaVal;
if ~isfield(plotInfo, 'circle_size'), plotInfo.circle_size = 40; end
if ~isfield(plotInfo, 'alpha_val'), plotInfo.alpha_val = 0.9; end
plot_brain3D_patch(plotInfo);

%% specify view
if isfield(plotInfo, 'viewAngle')
    view(plotInfo.viewAngle);
end
    
%% snapshots from different views
for v = 1:size(params.plot_brainTopo.model_views,1)
    thisView = params.plot_brainTopo.model_views(v,:);
    view(thisView(1),thisView(2));
    drawnow;
    
    % save snapshot
    set(f, 'InvertHardcopy','off');                 % preserves black background
    thisFigName = [figname '_view' num2str(v)];
    if plotInfo.printResolution == 0
        print(f, '-dpng','-r0', [outDir filesep thisFigName '.png']);
    else
        print(f, '-dpng','-r600', [outDir filesep thisFigName '.png']);
    end    
end

%% animation loop (animated gif)
if params.plot_brainTopo.plot_animation
    outDir_gif = [outDir filesep 'gif_animations'];
    if ~exist(outDir_gif, 'dir')
        mkdir(outDir_gif);
    end  
    filename = [outDir_gif filesep figname '.gif'];
    az = [-180:20:180];
    for n = 1:size(az,2)
        view(az(n),0);
        drawnow;
        frame = getframe(f);
        %frame = getframe(f, get(f, 'Position'));
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        if n == 1
            imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
        else
            imwrite(imind,cm,filename,'gif','WriteMode','append');
        end
    end
end

%% save
if ~isfield(plotInfo, 'saveFigLater')
    set(f, 'InvertHardcopy','off');                 % preserves black background
    set(f, 'PaperPositionMode','auto');
    %saveas(f, [outDir filesep figname '.fig']);
    if params.plot_brainTopo.printResolution == 0
        print(f, '-dpng','-r0', [outDir filesep figname '.png']);
    else
        print(f, '-dpng','-r600', [outDir filesep figname '.png']);
    end
    %close(f);    
    disp(['Figure: ' figname ' stored in: ' outDir]);
end



%% 3D brain model (Radek)
% p=patch(fv); % vytvoøení 3D modelu
% isonormals(x,y,z,V,p); % vyhlazení okrajù
