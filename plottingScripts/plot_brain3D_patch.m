function h = plot_brain3D_patch(plotInfo)
% plots semi-transparent 3D brain model into existing axes
% plots colorcoded channels as circles
% highlights selected brain areas

% (c) Jiri, Apr17, Nov18
% based on previous function brain3D_model.m

if ~isfield(plotInfo, 'visible_axis'), plotInfo.visible_axis = true; end
if ~isfield(plotInfo, 'plot_chVals_asCircles'), plotInfo.plot_chVals_asCircles = true; end
if ~isfield(plotInfo, 'plot_chVals_asClrPatch'), plotInfo.plot_chVals_asClrPatch = false; end
if ~isfield(plotInfo, 'plot_text_orientation_LR'), plotInfo.plot_text_orientation_LR = true; end
if isscalar(plotInfo.circle_size)
    plotInfo.circle_size = repmat(plotInfo.circle_size, [size(plotInfo.vals,1),1]);
else
    assert(size(plotInfo.circle_size,1) == size(plotInfo.vals,1));
end
if ~isfield(plotInfo, 'alpha_val'), plotInfo.alpha_val = 0.9; end
if isscalar(plotInfo.alpha_val)
    plotInfo.alpha_val = repmat(plotInfo.alpha_val, [size(plotInfo.vals,1),1]);
else
    assert(size(plotInfo.alpha_val,1) == size(plotInfo.vals,1));
end

%% required values
clrmap = plotInfo.clrmap;
fv = plotInfo.brain.fv;
alphaVal = plotInfo.alphaVal;
vals = plotInfo.vals;
clims = plotInfo.clims;
xi = plotInfo.brain.xi;
yi = plotInfo.brain.yi;
zi = plotInfo.brain.zi;
axColor = 'w';

%% ----------------channel values (face patch colors)------------
facecolor = repmat([1 1 1], length(fv.faces), 1);
facealpha = repmat(alphaVal, length(fv.faces), 1);
if plotInfo.plot_chVals_asClrPatch
    distanceThreshold = 5;     % in [mm]
    for ch = 1:size(plotInfo.chnlsMNI,2)    
        
        % index of colormap & color
        i_clr = cVals2cInds(vals(ch), [clims(1),clims(2)], size(clrmap.brain,1)+[1,size(clrmap.chnls,1)]);
        clr = clrmap.fig(i_clr,:);
        
        % MNI
        mni_x = plotInfo.chnlsMNI(1,ch); 
        mni_y = plotInfo.chnlsMNI(2,ch); 
        mni_z = plotInfo.chnlsMNI(3,ch);
        
        % maybe useful for debugging
        %scatter3(ix,iy,iz, circle_size, 'MarkerFaceColor',clr, 'MarkerEdgeColor','none');

        % find nearest vertices & faces
        pos_vrtx = find(abs(fv.vertices(:,1) - mni_y) <= distanceThreshold & ...
                       abs(fv.vertices(:,2) - mni_x) <= distanceThreshold & ...
                       abs(fv.vertices(:,3) - mni_z) <= distanceThreshold  );
        pos_face = find(sum(ismember(fv.faces, pos_vrtx),2)>0);
        
        % update face color and transparency 
        facecolor(pos_face,:) = repmat(clr, length(pos_face), 1);
        facealpha(pos_face,:) = repmat(plotInfo.alpha_val(ch), length(pos_face), 1);
%         if mod(ch,10) == 1
%             disp([' - brain 3D patch: ' num2str(ch) '/' num2str(size(plotInfo.chnlsMNI,2)) ' done.']);
%         end

    end
end

%% ----------- >>> brain: patch plot <<< ----------------------
% -----old version (as of 30.10.2018)---
% cData = ones(size(fv.vertices,1),1);
% brain_cInds = cVals2cInds(cData, [0,1], [1,size(clrmap.brain,1)]);
% p = patch('Faces', [fv.faces(:,2),fv.faces(:,1),fv.faces(:,3)], ...
%           'Vertices',[fv.vertices(:,2),fv.vertices(:,1),fv.vertices(:,3)], ...
%           'EdgeColor','none', 'CData',brain_cInds, 'CDataMapping','direct', ...
%           'FaceColor','interp','FaceAlpha',0.1);   % 'FaceAlpha',alphaVal
% ---------------------------------------
p = patch('Faces', [fv.faces(:,2),fv.faces(:,1),fv.faces(:,3)], ...
          'Vertices',[fv.vertices(:,2),fv.vertices(:,1),fv.vertices(:,3)], ...
          'EdgeColor','none', 'FaceVertexCData',facecolor, 'FaceColor','flat', ...
          'FaceVertexAlphaData',facealpha, 'FaceAlpha','flat', 'AlphaDataMapping','none');   % 'FaceAlpha',alphaVal
%isonormals(x,y,z,V,p);

%% ----------------areas----------------------
if isfield(plotInfo, 'aarea')
    for aa = 1:size(plotInfo.aarea,1)
        p_aa = patch('Faces', [plotInfo.aarea{aa}.fv.faces(:,2),plotInfo.aarea{aa}.fv.faces(:,1),plotInfo.aarea{aa}.fv.faces(:,3)], ...
                  'Vertices',[plotInfo.aarea{aa}.fv.vertices(:,2),plotInfo.aarea{aa}.fv.vertices(:,1),plotInfo.aarea{aa}.fv.vertices(:,3)], ...
                  'EdgeColor','none', 'FaceAlpha',0.1);          
        p_aa.FaceColor = plotInfo.aarea{aa}.faceClr; ...'red';
        %isonormals(x,y,z,V,p);
    end
end

%% ---------------plot properties-----------------
lighting phong
axis tight;
axis equal
view(0,90)
camlight headlight
material dull;
set(gca, 'color', plotInfo.fig_bkgClr);
if plotInfo.visible_axis
    set(gca, 'XColor',plotInfo.txtClr,'YColor',plotInfo.txtClr,'ZColor',plotInfo.txtClr);
    xlabel('x');
    ylabel('y');
    zlabel('z');
else
    set(gca, 'XColor',plotInfo.fig_bkgClr,'YColor',plotInfo.fig_bkgClr,'ZColor',plotInfo.fig_bkgClr);
    set(gca, 'XTick',[], 'YTick',[], 'ZTick',[]);
end

%% ----------------channels as circles----------------------
if plotInfo.plot_chVals_asCircles
    circle_size = plotInfo.circle_size;
    ch_scatter = nan(1,size(plotInfo.chnlsMNI,2));
    for ch = 1:size(plotInfo.chnlsMNI,2)    
        i_clr = cVals2cInds(vals(ch), [clims(1),clims(2)], size(clrmap.brain,1)+[1,size(clrmap.chnls,1)]);
        clr = clrmap.fig(i_clr,:);
    %     [ix,iy,iz] = mni2vox(plotInfo.chnls(ch).MNI_x, plotInfo.chnls(ch).MNI_y, plotInfo.chnls(ch).MNI_z, xi, yi, zi); % index of MNI coor
        mni_x = plotInfo.chnlsMNI(1,ch); mni_y = plotInfo.chnlsMNI(2,ch); mni_z = plotInfo.chnlsMNI(3,ch);
        ch_scatter(ch) = scatter3(mni_x,mni_y,mni_z, circle_size(ch), 'MarkerFaceColor',clr, 'MarkerEdgeColor','none');
    end
end

%% --------------orientation (L/R)----------------
if plotInfo.plot_text_orientation_LR
    xLims = get(gca,'xlim');
    yLims = get(gca,'ylim');
    txt_L = xLims(1) + diff(xLims)/100*5;        % ~ 5% offset from left  side
    txt_R = xLims(end) - diff(xLims)/100*5;      % ~ 5% offset from right side
    txt_U = yLims(end) - diff(yLims)/100*5;      % ~ 5% offset from upper side
    txt_D = yLims(1) + diff(yLims)/100*5;        % ~ 5% offset from lower side
%     h_txt_L = text(txt_L, txt_U, txt_D, 'L', 'fontsize',12, 'fontw','bold', 'HorizontalAlignment','center', 'Color',plotInfo.txtClr);
%     h_txt_R = text(txt_R, txt_U, txt_D, 'R', 'fontsize',12, 'fontw','bold', 'HorizontalAlignment','center', 'Color',plotInfo.txtClr);
    h_txt_L = text(txt_L, txt_D, txt_D, 'L', 'fontsize',12, 'fontw','bold', 'HorizontalAlignment','center', 'Color',plotInfo.txtClr);
    h_txt_R = text(txt_R, txt_D, txt_D, 'R', 'fontsize',12, 'fontw','bold', 'HorizontalAlignment','center', 'Color',plotInfo.txtClr);    
end

%% --------------other text (e.g. subjTag)----------------
if isfield(plotInfo, 'plot_text_tag')
    xLims = get(gca,'xlim');
    yLims = get(gca,'ylim');
    txt_L = xLims(1) + diff(xLims)/100*5;        % ~ 5% offset from left  side
    txt_R = xLims(end) - diff(xLims)/100*5;      % ~ 5% offset from right side
    txt_U = yLims(end) - diff(yLims)/100*5;      % ~ 5% offset from upper side
    txt_D = yLims(1) + diff(yLims)/100*5;        % ~ 5% offset from lower side
    h_txt_tag = text(txt_L, txt_U, txt_D, plotInfo.plot_text_tag, 'fontsize',12, 'fontw','bold', 'HorizontalAlignment','center', 'Color',plotInfo.txtClr);
end

%% update plot
drawnow;

%% return handles
h = struct;     % handles to graphics objects
h.patchBrain = p;
% h.ch_scatter = ch_scatter;
% h.txt_L = h_txt_L;
% h.txt_R = h_txt_R;

