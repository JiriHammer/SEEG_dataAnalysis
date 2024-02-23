function h = brain3D_plot(brainVolume, chnls_MNI_VAL, plotInfo)
% plots semi-transparent 3D brain model into existing axes
% plots colorcoded channels as circles
% highlights selected brain areas
% input
%   - plotInfo = params.plot_brain3D + some other fields (TO DO)
%   - brainVolumes = 
%   - chnls_MNI_VAL = 
% (c) Jiri, Nov21, Nov18, Apr17
% based on: plot_brain3D_patch.m


%% default channel value limits
if ~isfield(plotInfo, 'thisVolume'), plotInfo.thisVolume = 1; end
if ~isempty(chnls_MNI_VAL)
    if ~isfield(plotInfo, 'chVals_lims')
        chVals_lims = [prctile(chnls_MNI_VAL(:,4),5),prctile(chnls_MNI_VAL(:,4),95)];   % 5th & 95th percentile 
        if chVals_lims(2) <= chVals_lims(1)
            chVals_lims = [min(chnls_MNI_VAL(:,4)),max(chnls_MNI_VAL(:,4))];   % 5th & 95th percentile 
        end
    else
        chVals_lims = plotInfo.chVals_lims;
    end
end
assert(chVals_lims(2) > chVals_lims(1));

%% figure / text colors
if strcmp(plotInfo.backgroundColor, 'w')
    plotInfo.txtClr = 'k';
elseif strcmp(plotInfo.backgroundColor, 'k')
    plotInfo.txtClr = 'w';
else
    plotInfo.txtClr = 'c';
end

%% figure specified?
if ~isfield(plotInfo, 'figHandle')
    f = fig_make;
    subplot(1,1,1);
    hold on;
else
    f = plotInfo.figHandle;
end
set(f, 'visible','on', 'Color',plotInfo.backgroundColor);
set(f, 'InvertHardcopy','off');                 % preserves black background
set(f, 'Renderer','OpenGL');    
opengl('hardware');

%% init (first volume - can have color-coded values)
if plotInfo.thisVolume == 1
    fv = brainVolume.fv;
    face_clr = brainVolume.clrmap;              % color of patch face
    facecolor = repmat(face_clr, length(fv.faces), 1);
    face_alf = brainVolume.alfmap(end);         % transparency of patch face
    facealpha = repmat(face_alf, length(fv.faces), 1);  
end

%% -------map channel values onto brain surface (face patch colors) ?------
if plotInfo.chVals_asPatches && plotInfo.thisVolume == 1
    disp(' - 3D brain plotting: channel values are mapped onto first brain volume!');
    distanceThreshold = 5;     % in [mm], good value = 5
    ch_clrMap = plotInfo.chVals_colorMap;
    trn = plotInfo.chVals_patchTransp;     % transparency of changed values
    for ch = 1:size(chnls_MNI_VAL,1)    
        
        % index of colormap & color
        i_clr = cVals2cInds(chnls_MNI_VAL(ch,4), chVals_lims, [1,size(ch_clrMap,1)]);
        clr = ch_clrMap(i_clr,:);
        
        % MNI
        mni_x = chnls_MNI_VAL(ch,1); mni_y = chnls_MNI_VAL(ch,2); mni_z = chnls_MNI_VAL(ch,3);
        
        % maybe useful for debugging
        %scatter3(ix,iy,iz, circle_size, 'MarkerFaceColor',clr, 'MarkerEdgeColor','none');

        % find nearest vertices & faces
        pos_vrtx = find(abs(fv.vertices(:,1) - mni_y) <= distanceThreshold & ...
                       abs(fv.vertices(:,2) - mni_x) <= distanceThreshold & ...
                       abs(fv.vertices(:,3) - mni_z) <= distanceThreshold  );
        pos_face = find(sum(ismember(fv.faces, pos_vrtx),2)>0);
        
        % update face color and transparency 
        facecolor(pos_face,:) = repmat(clr, length(pos_face), 1);   % set color
        facealpha(pos_face,:) = repmat(trn, length(pos_face), 1);   % set transparency
    end
end

%% ----------- >>> brain: patch plot <<< ----------------------
fv = brainVolume.fv;
if plotInfo.thisVolume == 1
    p = patch('Faces', [fv.faces(:,2),fv.faces(:,1),fv.faces(:,3)], ...
              'Vertices',[fv.vertices(:,2),fv.vertices(:,1),fv.vertices(:,3)], ...
              'EdgeColor','none', 'FaceVertexCData',facecolor, 'FaceColor','flat', ...
              'FaceVertexAlphaData',facealpha, 'FaceAlpha','flat', 'AlphaDataMapping','none');   % 'FaceAlpha',alphaVal
    %isonormals(x,y,z,V,p);
else
    p = patch('Faces', [fv.faces(:,2),fv.faces(:,1),fv.faces(:,3)], ...
              'Vertices',[fv.vertices(:,2),fv.vertices(:,1),fv.vertices(:,3)], ...
              'EdgeColor','none', 'FaceAlpha',0.1);          
    p.FaceColor = brainVolume.clrmap;
    p.FaceAlpha = brainVolume.alfmap(end);
    %isonormals(x,y,z,V,p);
end        
h_patch = p;

%% ---------------plot properties-----------------
lighting phong
axis tight;
axis equal
% view(0,90)
view(plotInfo.viewAngle);
camlight headlight
material dull;
set(gca, 'color', plotInfo.backgroundColor);
if plotInfo.visible_axis
    set(gca, 'XColor',plotInfo.txtClr,'YColor',plotInfo.txtClr,'ZColor',plotInfo.txtClr);
    xlabel('x');
    ylabel('y');
    zlabel('z');
else
    set(gca, 'XColor',plotInfo.backgroundColor,'YColor',plotInfo.backgroundColor,'ZColor',plotInfo.backgroundColor);
    set(gca, 'XTick',[], 'YTick',[], 'ZTick',[]);
end

%% ----------------channels as circles----------------------
ch_scatter = [];
if plotInfo.chVals_asCircles && ~isempty(chnls_MNI_VAL)
    if isfield(plotInfo, 'def_circle_size')
        circle_size = plotInfo.def_circle_size;     % defined "from above"
    else
        circle_size = linTransform(chnls_MNI_VAL(:,4), chVals_lims, plotInfo.circleSizeLims);
    end
    ch_clrMap = plotInfo.chVals_colorMap;
    ch_scatter = nan(1,size(chnls_MNI_VAL,1));
    for ch = 1:size(chnls_MNI_VAL,1)    
        i_clr = cVals2cInds(chnls_MNI_VAL(ch,4), chVals_lims, [1,size(ch_clrMap,1)]);
        clr = ch_clrMap(i_clr,:);
    %     [ix,iy,iz] = mni2vox(plotInfo.chnls(ch).MNI_x, plotInfo.chnls(ch).MNI_y, plotInfo.chnls(ch).MNI_z, xi, yi, zi); % index of MNI coor
        mni_x = chnls_MNI_VAL(ch,1); mni_y = chnls_MNI_VAL(ch,2); mni_z = chnls_MNI_VAL(ch,3);
        
        ch_scatter(ch) = scatter3(mni_x,mni_y,mni_z, circle_size(ch), 'MarkerFaceColor',clr, 'MarkerEdgeColor','none');
%         ch_scatter(ch) = scatter3(mni_x,mni_y,mni_z, circle_size(ch), 'MarkerFaceColor','k', 'MarkerEdgeColor','k');        % only black
%         ch_scatter(ch) = scatter3(mni_x,mni_y,mni_z, circle_size(ch), 'MarkerFaceColor',[0.5 0.5 0.5], 'MarkerEdgeColor',[0.2 0.2 0.2]);        % only grey
%         disp(['ch = ' num2str(ch) ': val = ' num2str(chnls_MNI_VAL(ch,4)) ', size = ' num2str(circle_size(ch)) ', clr = ' num2str(clr)]);
    end
end

%% --------------orientation (L/R)----------------
if plotInfo.text_LR && plotInfo.thisVolume == 1
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
if ~isempty(plotInfo.text_tag) && plotInfo.thisVolume == 1
    xLims = get(gca,'xlim');
    yLims = get(gca,'ylim');
    txt_L = xLims(1) + diff(xLims)/100*5;        % ~ 5% offset from left  side
    txt_R = xLims(end) - diff(xLims)/100*5;      % ~ 5% offset from right side
    txt_U = yLims(end) - diff(yLims)/100*5;      % ~ 5% offset from upper side
    txt_D = yLims(1) + diff(yLims)/100*5;        % ~ 5% offset from lower side
    h_txt_tag = text(txt_L, txt_U, txt_D, plotInfo.text_tag, 'fontsize',12, 'fontw','bold', 'HorizontalAlignment','center', 'Color',plotInfo.txtClr);
end

%% update plot
drawnow;

%% return handles
h = struct;     % handles to graphics objects
h.patchBrain = h_patch;
h.ch_scatter = ch_scatter;
% h.txt_L = h_txt_L;
% h.txt_R = h_txt_R;
