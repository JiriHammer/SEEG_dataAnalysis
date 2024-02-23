function brainSlices_plot(brainVolumes, mni_vox, plotInfo)
% plots a single slice (axial, sagital, coronal) into specified axes or a new figure is created
% MRI (or CT) volumes are mapped to colormap/transparency based on their values (intensities)
% 
% input
%   - brainVolumes{vol} = 
%   struct with fields:
%             hdr: [1×1 struct]
%             vol: [177×209×176 double]
%             xyz: [3×6510768 double]
%              VI: [177×209×176 double]
%              xi: [1×177 double]
%              yi: [1×209 double]
%              zi: [1×176 double]
%     voxSize_new: 1
%          loaded: 1
%           cInds: [177×209×176 double]
%           aVals: [177×209×176 double]
%          vox_ix: 92
%          vox_iy: 120
%          vox_iz: 136
%   - mni_vox = MNI coordinates (x,y,z) for the cross-hair position
%   - plotInfo =  struct (must contain fields: colormap & plot_brainSlices)

% (c) Jiri, Nov21
% see also: brainSlices_vols_fixed.m & brainSlices_vols_chSpec.m
% based on previous: plot_ieegChnlsInBrain.m

%% defaults
if ~isfield(plotInfo, 'labelSlices'), plotInfo.labelSlices = true; end
if ~isfield(plotInfo, 'fontSize'), plotInfo.fontSize = 14; end
if ~isfield(plotInfo, 'txt_pos'), plotInfo.txt_pos = [0.016, 0.98]; end
if ~isfield(plotInfo, 'onlySlices'), plotInfo.onlySlices = false; end   % quick fix (TO DO: better axes allocation)

%% figure
if ~isfield(plotInfo, 'fig')
    f = fig_make;
    
    nRows = 1;
    nCols = size(plotInfo.plot_brainSlices.slicePlanes,2);
    nPlot = 1:nCols;
    marg_h = [0.05 0.05];
    marg_w = [0.04 0.04];
    gap = [0.03, 0.03];    
    fontSize = 14;
else
    f = plotInfo.fig;   % existing figure
    nRows = plotInfo.nRows;
    nCols = plotInfo.nCols;    
    nPlot = plotInfo.p_subPlot;
    marg_h = plotInfo.marg_h;
    marg_w = plotInfo.marg_w;
    gap = plotInfo.gap;    
    fontSize = plotInfo.fontSize;
end

%% settings
nVolumes = size(brainVolumes,2);    % number of brain volumes to plot
clrmap = plotInfo.clrmap;           % color map

%% MNI axis
mni_x = brainVolumes{1}.xi;      % interpolated x-axis, in [mm] of MNI coors
mni_y = brainVolumes{1}.yi;      % interpolated y-axis, in [mm] of MNI coors
mni_z = brainVolumes{1}.zi;      % interpolated z-axis, in [mm] of MNI coors

n = 1;

%% AXIAL (xy) slice
if ismember('axial', plotInfo.plot_brainSlices.slicePlanes)
    
    % axes
    if ~plotInfo.onlySlices
        subtightplot(nRows, nCols, nPlot(n,:), gap, marg_h, marg_w);
    else
        subtightplot(nRows, nCols, nPlot(:,n), gap, marg_h, marg_w);
    end
    set(gca,'Position', get(gca,'Position')+[0,-0.10,0.05,0.05]);
    hold on;    
    colormap(gca, clrmap);

    % >>> plot slice of brain volumes <<<
    for vol = 1:nVolumes
        if ~isempty(brainVolumes{vol})
            iz = brainVolumes{vol}.vox_iz;    % z-index of the slice
            V = brainVolumes{vol}.VI;         % volume slice (not really used for plotting)
            x = brainVolumes{vol}.xi;
            y = brainVolumes{vol}.yi;
            Z = ones(size(V,1),size(V,2));  % = 1

            cInds = brainVolumes{vol}.cInds(:,:,iz);          % coloring of indices
            aVals = brainVolumes{vol}.aVals(:,:,iz);          % transparency of indices
            h_vol = image(x, y, Z);                     % !!! plot slice !!!
            set(h_vol, 'CData', cInds', 'CDataMapping','direct', 'AlphaData',aVals');
        end
    end
    axis image;    

    % MNI coor (= channel)
    plot([-mni_vox(1), -mni_vox(1)],get(gca,'ylim'), 'y');
    plot(get(gca,'xlim'),[mni_vox(2), mni_vox(2)], 'y');

    % left / right orientation
    txt_L = mni_x(1) + (mni_x(end)-mni_x(1))/100*5;        % ~ 5% offset from left  side
    txt_R = mni_x(end) - (mni_x(end)-mni_x(1))/100*5;      % ~ 5% offset from right side
    txt_U = mni_y(end) - (mni_y(end)-mni_y(1))/100*5;      % ~ 5% offset from upper side
    text(txt_L, txt_U, 'R', 'fontsize',ceil(1.5*fontSize), 'fontw','bold', 'HorizontalAlignment','center', 'Color','w');
    text(txt_R, txt_U, 'L', 'fontsize',ceil(1.5*fontSize), 'fontw','bold', 'HorizontalAlignment','center', 'Color','w');

    % slices labels
    if plotInfo.labelSlices
        title(['AXIAL: MNI(z) = ' num2str(mni_vox(3))], 'FontWeight','bold', 'FontSize', ceil(1.2*fontSize));
        xlabel('x-MNI', 'FontSize', fontSize);
        ylabel('y-MNI', 'FontSize', fontSize);
        set(gca, 'FontSize', fontSize);
    else
        set(gca, 'XTick',[], 'YTick',[]);
    end
    n = n+1;
end

%% SAGITTAL (yz) slice
if ismember('sagittal', plotInfo.plot_brainSlices.slicePlanes)
    
    % axes
    if ~plotInfo.onlySlices
        subtightplot(nRows, nCols, nPlot(n,:), gap, marg_h, marg_w);
    else
        subtightplot(nRows, nCols, nPlot(:,n), gap, marg_h, marg_w);
    end
    set(gca,'Position', get(gca,'Position')+[0,-0.10,0.05,0.05]);
    hold on;    
    colormap(gca, clrmap);

    % >>> plot slice of brain volumes <<<
    for vol = 1:nVolumes
        if ~isempty(brainVolumes{vol})
            ix = brainVolumes{vol}.vox_ix;    % z-index of the slice
            V = brainVolumes{vol}.VI;         % volume slice (not really used for plotting)
            y = brainVolumes{vol}.yi;
            z = brainVolumes{vol}.zi;
            Z = ones(size(V,2),size(V,3));  % = 1

            cInds = squeeze(brainVolumes{vol}.cInds(ix,:,:));          % coloring of indices
            aVals = squeeze(brainVolumes{vol}.aVals(ix,:,:));          % transparency of indices
            h_vol = image(y, z, Z);                     % !!! plot slice !!!
            set(h_vol, 'CData', cInds', 'CDataMapping','direct', 'AlphaData',aVals');
        end
    end
    axis image;    

    % MNI coor (= channel)
    plot([mni_vox(2), mni_vox(2)],get(gca,'ylim'), 'y');
    plot(get(gca,'xlim'),[mni_vox(3), mni_vox(3)], 'y');

    % front / back orientation
    txt_L = mni_y(1) + (mni_y(end)-mni_y(1))/100*5;        % ~ 5% offset from left  side
    txt_R = mni_y(end) - (mni_y(end)-mni_y(1))/100*5;      % ~ 5% offset from right side
    txt_U = mni_z(end) - (mni_z(end)-mni_z(1))/100*5;      % ~ 5% offset from upper side
    text(txt_L, txt_U, 'B', 'fontsize',ceil(1.5*fontSize), 'fontw','bold', 'HorizontalAlignment','center', 'Color','w');
    text(txt_R, txt_U, 'F', 'fontsize',ceil(1.5*fontSize), 'fontw','bold', 'HorizontalAlignment','center', 'Color','w');

    % slices labels
    if plotInfo.labelSlices
        title(['SAGITTAL: MNI(x) = ' num2str(mni_vox(1))], 'FontWeight','bold', 'FontSize', ceil(1.2*fontSize));
        xlabel('y-MNI', 'FontSize', fontSize);
        ylabel('z-MNI', 'FontSize', fontSize);
        set(gca, 'FontSize', fontSize);
    else
        set(gca, 'XTick',[], 'YTick',[]);
    end
    n = n+1;
end

%% CORONAL (xz) slice
if ismember('coronal', plotInfo.plot_brainSlices.slicePlanes)
    
    % axes
    if ~plotInfo.onlySlices
        subtightplot(nRows, nCols, nPlot(n,:), gap, marg_h, marg_w);
    else
        subtightplot(nRows, nCols, nPlot(:,n), gap, marg_h, marg_w);
    end
    set(gca,'Position', get(gca,'Position')+[0,-0.10,0.05,0.05]);
    hold on;    
    colormap(gca, clrmap);

    % >>> plot slice of brain volumes <<<
    for vol = 1:nVolumes
        if ~isempty(brainVolumes{vol})
            iy = brainVolumes{vol}.vox_iy;    % z-index of the slice
            V = brainVolumes{vol}.VI;         % volume slice (not really used for plotting)
            x = brainVolumes{vol}.xi;
            z = brainVolumes{vol}.zi;
            Z = ones(size(V,1),size(V,3));  % = 1

            cInds = squeeze(brainVolumes{vol}.cInds(:,iy,:));          % coloring of indices
            aVals = squeeze(brainVolumes{vol}.aVals(:,iy,:));          % transparency of indices
            h_vol = image(-x, z, Z);                     % !!! plot slice !!!
            set(h_vol, 'CData', cInds', 'CDataMapping','direct', 'AlphaData',aVals');
        end
    end
    axis image;    

    % MNI coor (= channel)
    plot([mni_vox(1), mni_vox(1)],get(gca,'ylim'), 'y');
    plot(get(gca,'xlim'),[mni_vox(3), mni_vox(3)], 'y');

    % left / right orientation
    txt_L = mni_x(1) + (mni_x(end)-mni_x(1))/100*5;        % ~ 5% offset from left  side
    txt_R = mni_x(end) - (mni_x(end)-mni_x(1))/100*5;      % ~ 5% offset from right side
    txt_U = mni_z(end) - (mni_z(end)-mni_z(1))/100*5;      % ~ 5% offset from upper side
    text(txt_L, txt_U, 'L', 'fontsize',ceil(1.5*fontSize), 'fontw','bold', 'HorizontalAlignment','center', 'Color','w');
    text(txt_R, txt_U, 'R', 'fontsize',ceil(1.5*fontSize), 'fontw','bold', 'HorizontalAlignment','center', 'Color','w');

    % slices labels
    if plotInfo.labelSlices
        title(['CORONAL: MNI(y) = ' num2str(mni_vox(2))], 'FontWeight','bold', 'FontSize', ceil(1.2*fontSize));
        xlabel('x-MNI', 'FontSize', fontSize);
        ylabel('z-MNI', 'FontSize', fontSize);
        set(gca, 'FontSize', fontSize);
    else
        set(gca, 'XTick',[], 'YTick',[]);
    end
    set(gca, 'XDir','reverse');
    n = n+1;
end

%% text
if ~isfield(plotInfo, 'fig')
    if ~isfield(plotInfo.channel, 'neurologyLabel'), plotInfo.channel.neurologyLabel = 'n.a.'; end
    tx = axes('visible','off', 'position',[0 0 1 1]);
    mytitle = ['subject: ' plotInfo.subjTag ', channel = ' plotInfo.channel.name '(' num2str(plotInfo.thisCh) ')'...
        ', MNI = [' num2str(mni_vox,'% 3.0f') '], anatomical atlas = ' plotInfo.channel.ass_brainAtlas ...
        ', cytoarchitec.map = ' plotInfo.channel.ass_cytoarchMap ...
        ', neurologist = ' plotInfo.channel.neurologyLabel];
    mytitle = strrep(mytitle, '_','\_');
    text(plotInfo.txt_pos(1), plotInfo.txt_pos(2), mytitle, 'fontsize', 16, 'fontw', 'bold');
end

%% save fig
if ~isfield(plotInfo, 'fig')
    % output directory
    if isfield(plotInfo, 'outDir')
        outDir = plotInfo.outDir;
    else
        disp('output directory not specified, figure will not be saved.');
        return;
    end

    % figure name
    if isfield(plotInfo, 'figName')
        figname = plotInfo.figName;
    else
        figname = 'notNamed';
    end
    
    fig_save(f, figname, outDir, 'format',plotInfo.plot_brainSlices.printFormats, 'res',plotInfo.plot_brainSlices.printResolution);
    close(f);
end
