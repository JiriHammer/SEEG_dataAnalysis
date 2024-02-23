function [h_axes,h_image]  = clrBar_axes(plotInfo)
% plots "colorbar" as a new image into new axes of existing figure
% useful in case of two or more colormaps are catted together
% draws only specified range of the figure's colorbar

% (c) Jiri, Apr17

%% defaults
%if ~isfield(plotInfo, 'axPos'), plotInfo.axPos = [0.96 0.2 0.01 0.6]; end
if ~isfield(plotInfo, 'axPos'), plotInfo.axPos = [0.90 0.2 0.02 0.6]; end
if ~isfield(plotInfo, 'fig'), plotInfo.fig = gcf; end
if isfield(plotInfo, 'clims')
    clims = plotInfo.clims;
else
    error('no specification of color limits!');
end
if isfield(plotInfo, 'inds')
    inds = plotInfo.inds;
else
    inds = 1:size(get(plotInfo.fig, 'Colormap'),1);
    %error('no specification of colorbar indices to be used!');
end
if ~isfield(plotInfo, 'clrbar_label'), plotInfo.clrbar_label = ''; end
if ~isfield(plotInfo, 'txtClr'), plotInfo.txtClr = 'k'; end

%% colorbar values
clrbar_vals = [clims(1):1/1000*diff(clims):clims(2)]';
clrbar_xAx = 1:10;
clrbar_yAx = 1:size(clrbar_vals,1);
%clrbar_inds = cVals2cInds(repmat(clrbar_vals,[1,10]), [clims(1),clims(2)], size(clrmap.brain,1)+[1,size(clrmap.chnls,1)]);
clrbar_inds = cVals2cInds(repmat(clrbar_vals,[1,size(clrbar_xAx,2)]), [clims(1),clims(2)], [inds(1),inds(end)]);

%% axes position
h_axes = axes('visible','on', 'position',plotInfo.axPos);  % position
set(h_axes, 'Parent', plotInfo.fig)

%% image
h_image = image(clrbar_xAx, clrbar_yAx, ones(size(clrbar_inds)));
set(h_image, 'CData', clrbar_inds, 'CDataMapping','direct');

%% axes labels
minInd = min(clrbar_yAx); maxInd = max(clrbar_yAx); 
tickVals = round([minInd: (maxInd-minInd)/4 :maxInd]);                                                              % tick vals/labels
tickName = cell(1,length(tickVals));
for tick = 1:length(tickVals)
    x = clrbar_vals(tickVals(tick));
    if x < 0.005
        tickName{tick} = num2str(x, '%01.2e');    
    elseif x < 0.05
        tickName{tick} = num2str(x, '%01.3f');
    else
        tickName{tick} = num2str(x, '%01.1f');
    end
end
set(h_axes, 'YDir','normal');
set(h_axes, 'XTick',[], 'YTick',tickVals, 'YTickLabel',tickName);
set(h_axes, 'YAxisLocation','right');
ylabel(plotInfo.clrbar_label, 'Color',plotInfo.txtClr);
set(h_axes, 'XColor',plotInfo.txtClr,'YColor',plotInfo.txtClr,'ZColor',plotInfo.txtClr);
