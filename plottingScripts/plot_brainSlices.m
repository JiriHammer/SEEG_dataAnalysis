function plot_brainSlices(params, vals, plotInfo)
% plots topologically colorcoded values 'vals' of channels onto brain slices

% (c) Jiri, Jan17

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
alphaVal = 1.0;                                     % transparency of the colored values (1 = opaque)

%% figure name, position
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
voxSize_new = plotInfo.brain.voxSize_new;

%% area for colored channel values: enlarged voxel
enlargedVoxel.size_mm = params.plot_brainTopo.size_coloredCube;                             % size of the enlarged "voxel", in [mm]
enlargedVoxel.size_ind = abs(closestval(xi,enlargedVoxel.size_mm/2) - closestval(xi,0));    % half-size of the enlarged "voxel", in [indices]
if enlargedVoxel.size_ind == 0, enlargedVoxel.size_ind = 1; end
enlargedVoxel.side = -enlargedVoxel.size_ind:enlargedVoxel.size_ind;    % side of the cube, in [indices], w.r.t. its center

%% insert channel values at MNI coors to 3D arrays: color and transparency
chnls_cData = zeros(size(VI));                          % color data values
chnls_aData = zeros(size(VI));                          % fully transparent
mni_coors = [];
for ch = 1:size(plotInfo.chnlsMNI,2)
    mni_coors = cat(1, mni_coors, [-plotInfo.chnlsMNI(1,ch), plotInfo.chnlsMNI(2,ch), plotInfo.chnlsMNI(3,ch)]);
    [ix,iy,iz] = mni2vox(-plotInfo.chnlsMNI(1,ch), plotInfo.chnlsMNI(2,ch), plotInfo.chnlsMNI(3,ch), xi, yi, zi); % index of MNI coor
    
    i_sel_x = ix + enlargedVoxel.side;                  % x-indices of enlarged voxel
    i_sel_x(i_sel_x < 1) = [];                          % indices out of range
    i_sel_x(i_sel_x > length(xi)) = [];
    
    i_sel_y = iy + enlargedVoxel.side;                  % y-indices of enlarged voxel
    i_sel_y(i_sel_y < 1) = [];                          % indices out of range
    i_sel_y(i_sel_y > length(yi)) = [];    
    
    i_sel_z = iz + enlargedVoxel.side;                  % z-indices of enlarged voxel
    i_sel_z(i_sel_z < 1) = [];                          % indices out of range
    i_sel_z(i_sel_z > length(zi)) = [];    
    
    chnls_cData(i_sel_x,i_sel_y,i_sel_z) = vals(ch);    % insert the chnl value in a larger (cubic) voxel
    chnls_aData(i_sel_x,i_sel_y,i_sel_z) = alphaVal;    % transparency of the chnl value in a larger (cubic) voxel
end

%% map values to colormaps
% 3D brain data
%brain_cInds = cVals2cInds(VI, [min(VI(:)),max(VI(:))], [1,size(clrmap.brain,1)]);
brain_cInds = cVals2cInds(VI, [prctile(VI(:),1),prctile(VI(:),99)], [1,size(clrmap.brain,1)]);

% channel values
if ~isfield(plotInfo, 'chnl_clims')
    clims = [min(chnls_cData(:)),max(chnls_cData(:))];
else
    clims = plotInfo.chnl_clims;
end
chnls_cInds = cVals2cInds(chnls_cData, [clims(1),clims(2)], size(clrmap.brain,1)+[1,size(clrmap.chnls,1)]);
Z = ones(size(VI));     % dummy variable for image

% colorbar info
plotInfo.clims = clims;
plotInfo.inds = size(clrmap.brain,1)+[1,size(clrmap.chnls,1)];

%% slicing settings
sliceStep = 2.0;            % in [mm]
n_maxSlices = 120;          % max. number of slices per page

%% ------------------------- AXIAL SLICES ---------------------------------
if ismember('axial', params.plot_brainTopo.slicePlanes)
% --- number of slices
dim = 3;
min_mni = min(mni_coors(:,dim));
max_mni = max(mni_coors(:,dim));
i_sliceStep = closestval(zi, sliceStep) - closestval(zi, 0.0);
i_slices = closestval(zi, min_mni):i_sliceStep:closestval(zi, max_mni);    
n_slices = length(i_slices);
while n_slices > n_maxSlices
    sliceStep = sliceStep + voxSize_new;        % increase slice step
    i_sliceStep = closestval(zi, sliceStep) - closestval(zi, 0.0);
    i_slices = closestval(zi, min_mni):i_sliceStep:closestval(zi, max_mni);    
    n_slices = length(i_slices);
end

% --- figure
% f = figure('visible','on');
% set(f, 'Position', params.plot_brainTopo.figurePosition);
% set(f, 'Position', [1 41 1920 963]);
f = figure('units','normalized','outerposition',[0 0 1 1]);
set(f, 'Colormap', clrmap.fig);

nRows = floor(sqrt(n_slices/2));    % good coverage when: nCols = 2 * nRows
nCols = ceil(2*sqrt(n_slices/2));
if n_slices > nRows * nCols, nRows = nRows+1; end
assert(n_slices <= nRows * nCols);
marg_h = [0.05 0.05];
marg_w = [0.04 0.04];
gap = [0.005, 0.005];

% --- text
if isfield(plotInfo, 'text')
    tx = axes('visible','off', 'position',[0 0 1 1]);
    mytitle = plotInfo.text;
    mytitle = strrep(mytitle, '_','\_');
    text(0.016, 0.97, mytitle, 'fontsize', 16, 'fontw', 'bold');
end

% --- slices
for n = 1:n_slices
    iz = i_slices(end-n+1);
    
    % new subplot   
    %subplot(nRows, nCols, n);
    subtightplot(nRows, nCols, n, gap, marg_h, marg_w);
    hold on;    
    h_brain = image(xi, yi, Z(:,:,iz));
    set(h_brain, 'CData', brain_cInds(:,:,iz)');
    %h_chnls = image(-xi, yi, Z(:,:,iz));
    h_chnls = image(xi, yi, Z(:,:,iz));
    set(h_chnls, 'CData', chnls_cInds(:,:,iz)', 'AlphaData',chnls_aData(:,:,iz)');
    
%     title(['z = ' num2str(zi(iz))]);
    if ceil(n/nCols) == nRows && mod(n,nCols) == 1
        xlabel('x-MNI');
        ylabel('y-MNI');
    else
        set(gca, 'XTick',[], 'YTick',[]);
    end
    axis image;
    
    % click on plot to see it bigger
    set(gca, 'ButtonDownFcn','call_copy');   
    set(h_brain, 'ButtonDownFcn','call_copy');
    set(h_chnls, 'ButtonDownFcn','call_copy');
    
    % left/right orientation
    if n == 1
        txt_L = xi(1) + (xi(end)-xi(1))/100*5;        % ~ 5% offset from left  side
        txt_R = xi(end) - (xi(end)-xi(1))/100*5;      % ~ 5% offset from right side
        txt_U = yi(end) - (yi(end)-yi(1))/100*5;      % ~ 5% offset from upper side
        txt_D = yi(1) + (yi(end)-yi(1))/100*5;        % ~ 5% offset from lower side
        text(txt_L, txt_U, 'R', 'fontsize',12, 'fontw','bold', 'HorizontalAlignment','center', 'Color','w');
        text(txt_R, txt_U, 'L', 'fontsize',12, 'fontw','bold', 'HorizontalAlignment','center', 'Color','w');
    end    
    % MNI coor
    text(txt_L, txt_D,['z = ' num2str(zi(iz))], 'fontsize',12, 'fontw','bold', 'HorizontalAlignment','left', 'Color','w');
end

% --- colorbar
clrBar_axes(plotInfo);

% clrbar_axes = axes('visible','on', 'position',[0.97 0.2 0.01 0.6]);  % position
% clrbar_vals = [clims(1):1/1000*diff(clims):clims(2)]';
% clrbar_xAx = 1:10;
% clrbar_yAx = 1:size(clrbar_vals,1);
% clrbar_inds = cVals2cInds(repmat(clrbar_vals,[1,10]), [clims(1),clims(2)], size(clrmap.brain,1)+[1,size(clrmap.chnls,1)]);
% clrbar_hndl = image(clrbar_xAx, clrbar_yAx, ones(size(clrbar_inds)));
% set(clrbar_hndl, 'CData', clrbar_inds, 'CDataMapping','direct');
% minInd = min(clrbar_yAx); maxInd = max(clrbar_yAx); 
% tickVals = round([minInd: (maxInd-minInd)/4 :maxInd]);                                                              % tick vals/labels
% tickName = cell(1,length(tickVals));
% for tick = 1:length(tickVals)
%     tickName{tick} = num2str(clrbar_vals(tickVals(tick)), '%01.1f');
% end
% % set(clrbar_axes, 'XTick',[], 'YTick',tickVals, 'YTickLabel',tickName);
% set(gca, 'YDir','normal');
% set(gca, 'XTick',[], 'YTick',tickVals, 'YTickLabel',tickName);
% set(clrbar_axes, 'YAxisLocation','right');

% --- output directory
if isfield(plotInfo, 'outDir')
    outDir = [plotInfo.outDir  filesep 'slices_axial'];
else
    outDir = [params.storage.outputDir filesep 'slices_axial'];
end
if ~exist(outDir, 'dir')
    mkdir(outDir);
end  

% --- save
set(f, 'PaperPositionMode','auto');
saveas(f, [outDir filesep figname '.fig']);
if params.plot_brainTopo.printResolution == 0
    print(f, '-dpng','-r0', [outDir filesep figname '.png']);
else
    print(f, '-dpng','-r600', [outDir filesep figname '.png']);
end
close(f);
display(['Figure: ' figname ' stored in: ' outDir]);
end

%% ------------------------- SAGITTAL SLICES ---------------------------------
if ismember('sagittal', params.plot_brainTopo.slicePlanes)
% --- number of slices
dim = 1;
min_mni = min(mni_coors(:,dim));
max_mni = max(mni_coors(:,dim));
i_sliceStep = closestval(xi, sliceStep) - closestval(xi, 0.0);
i_slices = closestval(xi, min_mni):i_sliceStep:closestval(xi, max_mni);    
n_slices = length(i_slices);
while n_slices > n_maxSlices
    sliceStep = sliceStep + voxSize_new;        % increase slice step
    i_sliceStep = closestval(xi, sliceStep) - closestval(xi, 0.0);
    i_slices = closestval(xi, min_mni):i_sliceStep:closestval(xi, max_mni);    
    n_slices = length(i_slices);
end

% --- figure
f = figure('units','normalized','outerposition',[0 0 1 1]);
% f = figure('visible','on');
% set(f, 'Position', params.plot_brainTopo.figurePosition);
% set(f, 'Position', [1 41 1920 963]);
set(f, 'Colormap', clrmap.fig);

nRows = floor(sqrt(n_slices/2));    % good coverage when: nCols = 2 * nRows
nCols = ceil(2*sqrt(n_slices/2));
if n_slices > nRows * nCols, nRows = nRows+1; end
assert(n_slices <= nRows * nCols);
marg_h = [0.05 0.05];
marg_w = [0.04 0.04];
gap = [0.005, 0.005];

% --- text
if isfield(plotInfo, 'text')
    tx = axes('visible','off', 'position',[0 0 1 1]);
    mytitle = plotInfo.text;
    mytitle = strrep(mytitle, '_','\_');
    text(0.016, 0.97, mytitle, 'fontsize', 16, 'fontw', 'bold');
end

% --- slices
for n = 1:n_slices
    ix = i_slices(n);
    
    % new subplot   
    %subplot(nRows, nCols, n);
    subtightplot(nRows, nCols, n, gap, marg_h, marg_w);
    hold on;    
    h_brain = image(yi, zi, squeeze(Z(ix,:,:)));
    set(h_brain, 'CData', squeeze(brain_cInds(ix,:,:))');
    h_chnls = image(yi, zi, squeeze(Z(ix,:,:)));
    set(h_chnls, 'CData', squeeze(chnls_cInds(ix,:,:))', 'AlphaData',squeeze(chnls_aData(ix,:,:))');
    
%     title(['x = ' num2str(-xi(ix))]);
    if ceil(n/nCols) == nRows && mod(n,nCols) == 1
        xlabel('x-MNI');
        ylabel('z-MNI');
    else
        set(gca, 'XTick',[], 'YTick',[]);
    end
    axis image;
    
    % click on plot to see it bigger
    set(gca, 'ButtonDownFcn','call_copy');   
    set(h_brain, 'ButtonDownFcn','call_copy');
    set(h_chnls, 'ButtonDownFcn','call_copy');
    
    % left/right orientation
    if n == 1
        txt_L = yi(1) + (yi(end)-yi(1))/100*5;        % ~ 5% offset from left  side
        txt_M = yi(1) + (yi(end)-yi(1))/100*50;       % ~ 50% offset from left  side
        txt_U = zi(end) - (zi(end)-zi(1))/100*5;      % ~ 5% offset from upper side
        txt_D = zi(1) + (zi(end)-zi(1))/100*5;        % ~ 5% offset from lower side
        text(txt_M, txt_U, 'L', 'fontsize',12, 'fontw','bold', 'HorizontalAlignment','center', 'Color','w');
    end    
    if n == n_slices
        text(txt_M, txt_U, 'R', 'fontsize',12, 'fontw','bold', 'HorizontalAlignment','center', 'Color','w');
    end
    % MNI coor
    text(txt_L, txt_D,['x = ' num2str(xi(ix))], 'fontsize',12, 'fontw','bold', 'HorizontalAlignment','left', 'Color','w');
end

% --- colorbar
clrBar_axes(plotInfo);

% --- output directory
if isfield(plotInfo, 'outDir')
    outDir = [plotInfo.outDir  filesep 'slices_sagittal'];
else
    outDir = [params.storage.outputDir filesep 'slices_sagittal'];
end
if ~exist(outDir, 'dir')
    mkdir(outDir);
end  

% --- save
set(f, 'PaperPositionMode','auto');
saveas(f, [outDir filesep figname '.fig']);
if params.plot_brainTopo.printResolution == 0
    print(f, '-dpng','-r0', [outDir filesep figname '.png']);
else
    print(f, '-dpng','-r600', [outDir filesep figname '.png']);
end
close(f);    
display(['Figure: ' figname ' stored in: ' outDir]);
end

%% ------------------------- CORONAL SLICES ---------------------------------
if ismember('coronal', params.plot_brainTopo.slicePlanes)
% --- number of slices
dim = 2;
min_mni = min(mni_coors(:,dim));
max_mni = max(mni_coors(:,dim));
i_sliceStep = closestval(yi, sliceStep) - closestval(yi, 0.0);
i_slices = closestval(yi, min_mni):i_sliceStep:closestval(yi, max_mni);    
n_slices = length(i_slices);
while n_slices > n_maxSlices
    sliceStep = sliceStep + voxSize_new;        % increase slice step
    i_sliceStep = closestval(yi, sliceStep) - closestval(yi, 0.0);
    i_slices = closestval(yi, min_mni):i_sliceStep:closestval(yi, max_mni);    
    n_slices = length(i_slices);
end

% --- figure
f = figure('units','normalized','outerposition',[0 0 1 1]);
% f = figure('visible','on');
% set(f, 'Position', params.plot_brainTopo.figurePosition);
% set(f, 'Position', [1 41 1920 963]);
set(f, 'Colormap', clrmap.fig);

nRows = floor(sqrt(n_slices/2));    % good coverage when: nCols = 2 * nRows
nCols = ceil(2*sqrt(n_slices/2));
if n_slices > nRows * nCols, nRows = nRows+1; end
assert(n_slices <= nRows * nCols);
marg_h = [0.05 0.05];
marg_w = [0.04 0.04];
gap = [0.005, 0.005];

% --- text
if isfield(plotInfo, 'text')
    tx = axes('visible','off', 'position',[0 0 1 1]);
    mytitle = plotInfo.text;
    mytitle = strrep(mytitle, '_','\_');
    text(0.016, 0.97, mytitle, 'fontsize', 16, 'fontw', 'bold');
end

% --- slices
for n = 1:n_slices
    iy = i_slices(n);
    
    % new subplot   
    %subplot(nRows, nCols, n);
    subtightplot(nRows, nCols, n, gap, marg_h, marg_w);
    hold on;    
    h_brain = image(xi, zi, squeeze(Z(:,iy,:)));
    set(h_brain, 'CData', squeeze(brain_cInds(:,iy,:))');
    h_chnls = image(xi, zi, squeeze(Z(:,iy,:)));
    set(h_chnls, 'CData', squeeze(chnls_cInds(:,iy,:))', 'AlphaData',squeeze(chnls_aData(:,iy,:))');
    
%     title(['y = ' num2str(yi(iy))]);
    if ceil(n/nCols) == nRows && mod(n,nCols) == 1
        xlabel('x-MNI');
        ylabel('z-MNI');
    else
        set(gca, 'XTick',[], 'YTick',[]);
    end
    axis image;
    
    % click on plot to see it bigger
    set(gca, 'ButtonDownFcn','call_copy');   
    set(h_brain, 'ButtonDownFcn','call_copy');
    set(h_chnls, 'ButtonDownFcn','call_copy');
    
    % left/right orientation
    if n == 1
        txt_L = xi(1) + (xi(end)-xi(1))/100*5;        % ~ 5% offset from left  side
        txt_R = xi(end) - (xi(end)-xi(1))/100*5;      % ~ 5% offset from right side
        txt_U = zi(end) - (zi(end)-zi(1))/100*5;      % ~ 5% offset from upper side
        txt_D = zi(1) + (zi(end)-zi(1))/100*5;        % ~ 5% offset from lower side
        text(txt_L, txt_U, 'R', 'fontsize',12, 'fontw','bold', 'HorizontalAlignment','center', 'Color','w');
        text(txt_R, txt_U, 'L', 'fontsize',12, 'fontw','bold', 'HorizontalAlignment','center', 'Color','w');
    end    
    % MNI coor
    text(txt_L, txt_D,['y = ' num2str(yi(iy))], 'fontsize',12, 'fontw','bold', 'HorizontalAlignment','left', 'Color','w');    
end

% --- colorbar
clrBar_axes(plotInfo);

% --- output directory
if isfield(plotInfo, 'outDir')
    outDir = [plotInfo.outDir  filesep 'slices_coronal'];
else
    outDir = [params.storage.outputDir filesep 'slices_coronal'];
end
if ~exist(outDir, 'dir')
    mkdir(outDir);
end  

% --- save
set(f, 'PaperPositionMode','auto');
saveas(f, [outDir filesep figname '.fig']);
if params.plot_brainTopo.printResolution == 0
    print(f, '-dpng','-r0', [outDir filesep figname '.png']);
else
    print(f, '-dpng','-r600', [outDir filesep figname '.png']);
end
close(f);    
display(['Figure: ' figname ' stored in: ' outDir]);
end










%% test: delete
% ch = 1;
% thisCh = selCh_H(ch);                               % chnl index in H
% [ix,iy,iz] = mni2vox(H.channels(thisCh).MNI_x, H.channels(thisCh).MNI_y, H.channels(thisCh).MNI_z, xi, yi, zi); % index of MNI coor
% 
% figure;
% set(gcf, 'Position', [680 -375 1827 1353]);
% set(gcf, 'Colormap', clrmap.fig);
% 
% subplot(1,2,1); hold on;
% h_brain = image(-xi, yi, Z(:,:,iz));
% set(h_brain, 'CData', brain_cInds(:,:,iz)', 'CDataMapping','direct');
% %caxis([1 size(get(gcf,'ColorMap'),1)]);             % strech the climits for the whole figure's color map
% plot(H.channels(thisCh).MNI_x, H.channels(thisCh).MNI_y, 'ro', 'LineWidth',2,'MarkerSize',12);
% title('interp');
% xlabel('x-MNI');
% ylabel('y-MNI');
% axis image;
% colorbar;
% 
% subplot(1,2,2); hold on;
% h_brain = image(-xi, yi, Z(:,:,iz));
% set(h_brain, 'CData', brain_cInds(:,:,iz)', 'CDataMapping','direct');
% h_chnls = image(xi, yi, Z(:,:,iz));
% set(h_chnls, 'CData', chnls_cInds(:,:,iz)', 'AlphaData',chnls_aData(:,:,iz)', 'CDataMapping','direct');
% 
% %caxis([1 size(get(gcf,'ColorMap'),1)]);             % strech the climits for the whole figure's color map
% plot(H.channels(thisCh).MNI_x, H.channels(thisCh).MNI_y, 'ro', 'LineWidth',2,'MarkerSize',12);
% title('interp');
% xlabel('x-MNI');
% ylabel('y-MNI');
% axis image;
% colorbar;

%% useful? design convolution kernel: 3D sphere
% kernel.size_mm = 2;
% kernel.size_vox = abs(closestval(xi,kernel.size_mm) - closestval(xi,0));
% kernel.side = -kernel.size_vox:kernel.size_vox;             % side of the cube, in [indices]
% kernel.center = closestval(kernel.side,0);                  % index of cube center
% kernel.vals = zeros(size(kernel.side,2),size(kernel.side,2),size(kernel.side,2));        % 3D cube values
% for m = 1:size(kernel.vals,1)
%     for n = 1:size(kernel.vals,2)
%         for p = 1:size(kernel.vals,3)
%             if sqrt((m-kernel.center)^2+(n-kernel.center)^2+(p-kernel.center)^2) <= kernel.size_vox
%                 kernel.vals(m,n,p) = 1;                     % within sphere = 1
%             end
%         end
%     end
% end

%% convolve brain with electrode values in 3D: too slow ...
% vals_xyz = zeros(size(VI));
% for ch = 1:size(selCh_H,2)
%     thisCh = selCh_H(ch);
%     tmp = zeros(size(VI));
%     [ix,iy,iz] = mni2vox(H.channels(thisCh).MNI_x, H.channels(thisCh).MNI_y, H.channels(thisCh).MNI_z, xi, yi, zi);
%     tmp(ix,iy,iz) = 1;
%     tmp = convn(tmp, kernel.vals .* vals(ch), 'same');
%     vals_xyz = vals_xyz + tmp;
%     display(['channel: ' num2str(ch) '/' num2str(size(selCh_H,2)) ' done.']);
% end