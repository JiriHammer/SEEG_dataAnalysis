function plot_ieegChnlsInBrain(params, plotInfo)
% plots channels in MNI coordinates in slices (axial, sagittal, coronal)
% uses radiology convention for slices (mirrored axial slices)
% as default plots MRI, optionally CT (or post-implant MRI)
% (c) Jiri, Jan17

%% defaults
if ~isfield(plotInfo, 'labelSlices'), plotInfo.labelSlices = true; end
if ~isfield(plotInfo, 'fontSize'), plotInfo.fontSize = 14; end
if ~isfield(plotInfo, 'txt_pos'), plotInfo.txt_pos = [0.016, 0.98]; end

%% color maps
if ~isfield(plotInfo, 'clrmap')
    clrmap.mr = gray(128);                              % colormap for MRI brain (T1, T2, ...)
    inds_mr = [1,size(clrmap.mr,1)];                    % indices to figure's colormap
    clrmap.ct = copper(128);                            % colormap for CT
    inds_ct = size(clrmap.mr,1)+[1,size(clrmap.ct,1)];  % indices to figure's colormap
    clrmap.fig = cat(1, clrmap.mr, clrmap.ct);          % colormap of the figure
else
    clrmap.mr = plotInfo.clrmap.mr;                     % colormap for MRI brain (T1, T2, ...)
    inds_mr = plotInfo.clrmap.inds_mr;                  % indices to figure's colormap
    clrmap.ct = plotInfo.clrmap.ct;                     % colormap for CT
    inds_ct = plotInfo.clrmap.inds_ct;                  % indices to figure's colormap    
    clrmap.fig = plotInfo.clrmap.fig;                   % colormap of the figure
end
alphaVal = 0.8;                                     % max. transparency of the CT

%% pass info from loaded brain MRI (T1) (see getBrainData.m)
mr_V = plotInfo.brain_T1.VI;      % interpolated volume
mr_x = plotInfo.brain_T1.xi;      % interpolated x-axis, in [mm] of MNI coors
mr_y = plotInfo.brain_T1.yi;      % interpolated y-axis, in [mm] of MNI coors
mr_z = plotInfo.brain_T1.zi;      % interpolated z-axis, in [mm] of MNI coors

%% pass info from loaded brain CT (see getBrainData.m)
ct_V = plotInfo.brain_CT.VI;      % interpolated volume
ct_x = plotInfo.brain_CT.xi;      % interpolated x-axis, in [mm] of MNI coors
ct_y = plotInfo.brain_CT.yi;      % interpolated y-axis, in [mm] of MNI coors
ct_z = plotInfo.brain_CT.zi;      % interpolated z-axis, in [mm] of MNI coors

%% map values to colormaps
% T1 3D brain data: map to colormap indices
%mr_cInds = cVals2cInds(mr_V, [min(mr_V(:)),max(mr_V(:))], [1,size(clrmap.mr,1)]);
mr_cInds = cVals2cInds(mr_V, [prctile(mr_V(:),1),prctile(mr_V(:),99)], inds_mr);

% CT 3D brain data: map to colormap indices
%ct_cInds = cVals2cInds(ct_V, [min(ct_V(:)),max(ct_V(:))], [1,size(clrmap_fig,1)]);
ct_cInds = cVals2cInds(ct_V, [prctile(ct_V(:),1),prctile(ct_V(:),99)], inds_ct);

% CT 3D brain data: map to transparency indices (small values are invisible)
ct_aVals = linTransform(ct_V, [prctile(ct_V(:),1),prctile(ct_V(:),99.9)], [0,alphaVal]);

%% MNI coors -> voxel index values
mni = [plotInfo.channel.MNI_x, plotInfo.channel.MNI_y, plotInfo.channel.MNI_z];
[mr_ix,mr_iy,mr_iz] = mni2vox(-mni(1), mni(2), mni(3), mr_x, mr_y, mr_z);
[ct_ix,ct_iy,ct_iz] = mni2vox(-mni(1), mni(2), mni(3), ct_x, ct_y, ct_z);

%% figure
if ~isfield(plotInfo, 'fig')
    f = figure('visible','on');
    %set(f, 'Position', [1 -479 2880 1472]);
    %set(f, 'Position', [1 41 1920 963]);
    set(f, 'units','normalized','outerposition',[0 0 1 1]);
    
    nRows = 1;
    nCols = 3;
    nPlot = [1 2 3];
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
set(f, 'Colormap', clrmap.fig)

%% AXIAL (xy) slice
% subplot(1,3,1)
% subplot(nRows, nCols, nPlot(1));
% subtightplot(nRows, nCols, [1,4], gap, marg_h, marg_w);
subtightplot(nRows, nCols, nPlot(1), gap, marg_h, marg_w);
set(gca,'Position', get(gca,'Position')+[0,-0.10,0.05,0.05]);
hold on;    

% MRI
Z = ones(size(mr_V(:,:,mr_iz)));
h_t1 = image(mr_x, mr_y, Z);
set(h_t1, 'CData', mr_cInds(:,:,mr_iz)', 'CDataMapping','direct');

% CT
Z = ones(size(mr_V(:,:,ct_iz)));
h_ct = image(ct_x, ct_y, Z);
set(h_ct, 'CData', ct_cInds(:,:,ct_iz)', 'CDataMapping','direct', 'AlphaData',ct_aVals(:,:,ct_iz)');
axis image;

% MNI coor (= channel)
plot([-mni(1), -mni(1)],get(gca,'ylim'), 'y');
plot(get(gca,'xlim'),[mni(2), mni(2)], 'y');
%plot(-mni(1), mni(2), 'ro', 'LineWidth',2,'MarkerSize',12);

% left / right orientation
txt_L = mr_x(1) + (mr_x(end)-mr_x(1))/100*5;        % ~ 5% offset from left  side
txt_R = mr_x(end) - (mr_x(end)-mr_x(1))/100*5;      % ~ 5% offset from right side
txt_U = mr_y(end) - (mr_y(end)-mr_y(1))/100*5;      % ~ 5% offset from upper side
text(txt_L, txt_U, 'R', 'fontsize',ceil(1.5*fontSize), 'fontw','bold', 'HorizontalAlignment','center', 'Color','w');
text(txt_R, txt_U, 'L', 'fontsize',ceil(1.5*fontSize), 'fontw','bold', 'HorizontalAlignment','center', 'Color','w');

%axis image;
if plotInfo.labelSlices
    title(['AXIAL: MNI(z) = ' num2str(mni(3))], 'FontWeight','bold', 'FontSize', ceil(1.2*fontSize));
    xlabel('x-MNI', 'FontSize', fontSize);
    ylabel('y-MNI', 'FontSize', fontSize);
    set(gca, 'FontSize', fontSize);
else
    set(gca, 'XTick',[], 'YTick',[]);
end

% click on plot to see it bigger
set(gca, 'ButtonDownFcn','call_copy');   
set(h_t1, 'ButtonDownFcn','call_copy');
set(h_ct, 'ButtonDownFcn','call_copy');

%% SAGITTAL (yz) slice
% subplot(1,3,2);
% subplot(nRows, nCols, nPlot(2));
% subtightplot(nRows, nCols, [2,5], gap, marg_h, marg_w);
subtightplot(nRows, nCols, nPlot(2), gap, marg_h, marg_w);
set(gca,'Position', get(gca,'Position')+[0,-0.10,0.05,0.05]);
hold on;    

% MRI
Z = ones(size(squeeze(mr_V(mr_ix,:,:))));
h_t1 = image(mr_y, mr_z, Z);
set(h_t1, 'CData', squeeze(mr_cInds(mr_ix,:,:))');

% CT
Z = ones(size(squeeze(mr_V(ct_ix,:,:))));
h_ct = image(ct_y, ct_z, Z);
set(h_ct, 'CData', squeeze(ct_cInds(ct_ix,:,:))', 'AlphaData',squeeze(ct_aVals(ct_ix,:,:))');
axis image;

% MNI coor (= channel)
plot([mni(2), mni(2)],get(gca,'ylim'), 'y');
plot(get(gca,'xlim'),[mni(3), mni(3)], 'y');
%plot(mni(2), mni(3), 'ro', 'LineWidth',2,'MarkerSize',12);

% front / back orientation
txt_L = mr_y(1) + (mr_y(end)-mr_y(1))/100*5;        % ~ 5% offset from left  side
txt_R = mr_y(end) - (mr_y(end)-mr_y(1))/100*5;      % ~ 5% offset from right side
txt_U = mr_z(end) - (mr_z(end)-mr_z(1))/100*5;      % ~ 5% offset from upper side
text(txt_L, txt_U, 'B', 'fontsize',ceil(1.5*fontSize), 'fontw','bold', 'HorizontalAlignment','center', 'Color','w');
text(txt_R, txt_U, 'F', 'fontsize',ceil(1.5*fontSize), 'fontw','bold', 'HorizontalAlignment','center', 'Color','w');

%axis image;
if plotInfo.labelSlices
    title(['SAGITTAL: MNI(x) = ' num2str(mni(1))], 'FontWeight','bold', 'FontSize', ceil(1.2*fontSize));
    xlabel('y-MNI', 'FontSize', fontSize);
    ylabel('z-MNI', 'FontSize', fontSize);
    set(gca, 'FontSize', fontSize);
else
    set(gca, 'XTick',[], 'YTick',[]);
end

% click on plot to see it bigger
set(gca, 'ButtonDownFcn','call_copy');   
set(h_t1, 'ButtonDownFcn','call_copy');
set(h_ct, 'ButtonDownFcn','call_copy');

%% CORONAL (xz) slice
% subplot(1,3,3);
% subplot(nRows, nCols, nPlot(3));
% subtightplot(nRows, nCols, [3,6], gap, marg_h, marg_w);
subtightplot(nRows, nCols, nPlot(3), gap, marg_h, marg_w);
set(gca,'Position', get(gca,'Position')+[0,-0.10,0.05,0.05]);
hold on;    

% MRI
Z = ones(size(squeeze(mr_V(:,mr_iy,:))));
h_t1 = image(-mr_x, mr_z, Z);
set(h_t1, 'CData', squeeze(mr_cInds(:,mr_iy,:))');

% CT
Z = ones(size(squeeze(mr_V(:,ct_iy,:))));
h_ct = image(-ct_x, ct_z, Z);
set(h_ct, 'CData', squeeze(ct_cInds(:,ct_iy,:))', 'AlphaData',squeeze(ct_aVals(:,ct_iy,:))');
axis image;

% MNI coor (= channel)
plot([mni(1), mni(1)],get(gca,'ylim'), 'y');
plot(get(gca,'xlim'),[mni(3), mni(3)], 'y');
%plot(mni(1), mni(3), 'ro', 'LineWidth',2,'MarkerSize',12);

% left / right orientation
txt_L = mr_x(1) + (mr_x(end)-mr_x(1))/100*5;        % ~ 5% offset from left  side
txt_R = mr_x(end) - (mr_x(end)-mr_x(1))/100*5;      % ~ 5% offset from right side
txt_U = mr_z(end) - (mr_z(end)-mr_z(1))/100*5;      % ~ 5% offset from upper side
text(txt_L, txt_U, 'L', 'fontsize',ceil(1.5*fontSize), 'fontw','bold', 'HorizontalAlignment','center', 'Color','w');
text(txt_R, txt_U, 'R', 'fontsize',ceil(1.5*fontSize), 'fontw','bold', 'HorizontalAlignment','center', 'Color','w');

%axis image;
if plotInfo.labelSlices
    title(['CORONAL: MNI(y) = ' num2str(mni(2))], 'FontWeight','bold', 'FontSize', ceil(1.2*fontSize));
    xlabel('x-MNI', 'FontSize', fontSize);
    ylabel('z-MNI', 'FontSize', fontSize);
    set(gca, 'FontSize', fontSize);
else
    set(gca, 'XTick',[], 'YTick',[]);
end

set(gca, 'XDir','reverse');

% click on plot to see it bigger
set(gca, 'ButtonDownFcn','call_copy');   
set(h_t1, 'ButtonDownFcn','call_copy');
set(h_ct, 'ButtonDownFcn','call_copy');

%% text
if ~isfield(plotInfo, 'fig')
    if ~isfield(plotInfo.channel, 'neurologyLabel'), plotInfo.channel.neurologyLabel = 'n.a.'; end
    tx = axes('visible','off', 'position',[0 0 1 1]);
    mytitle = ['subject: ' params.storage.subjTag ', channel = ' plotInfo.channel.name '(' num2str(plotInfo.thisCh) ')'...
        ', MNI = [' num2str(mni,'% 3.0f') '], anatomical atlas = ' plotInfo.channel.ass_brainAtlas ...
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
        outDir = params.storage.outputDir;
    end
    if ~exist(outDir, 'dir')
        mkdir(outDir);
    end  

    % figure name
    if isfield(plotInfo, 'figName')
        figname = plotInfo.figName;
    else
        figname = 'notNamed';
    end

    % save
    set(f, 'PaperPositionMode','auto');
%     saveas(f, [outDir filesep figname '.fig']);
    if params.plot_brainTopo.printResolution == 0
        print(f, '-dpng','-r0', [outDir filesep figname '.png']);
    else
        print(f, '-dpng','-r600', [outDir filesep figname '.png']);
    end
    close(f);    
    display(['Figure: ' figname ' stored in: ' outDir]);
end




