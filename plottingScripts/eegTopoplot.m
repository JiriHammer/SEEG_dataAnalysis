function eegTopoplot(H_topo, mapVals, info)
% plots head topology of EEG values

% (c) Jiri, Aug16

% figure
if ~isfield(info.plot, 'figHandle')
    f = figure('visible','on', 'Position',[1, 1, 1920, 1200]);
    set(f, 'PaperUnits','inches', 'PaperPosition',[0 0 19 12]);
else
    f = info.plot.figHandle;
end

% axes
axPosition = info.plot.axPosition;
set(gcf,'color', 'w');
ax1 = axes('Parent',f, 'Position',axPosition, 'NextPlot','replacechildren');
axis square;
set(ax1, 'XTick',[], 'YTick',[]);

% color limits
% maxC = max(abs(mapVals));
% clim = [-maxC maxC];
clim = info.plot.clims;
colormap(redblue(256));

% topoplot: 'mapVals' = array with 1 power value for each channel
h = posplot(H_topo ,mapVals,'-delaunay','-clim',clim); % creates the topoplot

hold on;
axis off;
set(gca,'DataAspectRatio',[1 1 1])

% some optional things, head circle, nose, points for electrodes ...
drawCircle3d([mean(H_topo.cp(:,1)),mean(H_topo.cp(:,2))*1.05,1,.459],'LineWidth',2,'Color','k')  %head circle
plot3([0.4479, 0.5002, 0.5525],[0.9872, 1.1, 0.9872]*0.97,[1, 1, 1],'k','LineWidth',2) % nose

h = posplot(H_topo,mapVals,'-delaunay','-clim',clim); % somehow, I have to plot this again after the circle/nose... don't know why...

% plot all electrodes as black points
h1=plot3(H_topo.cp(:,1), H_topo.cp(:,2),ones(length(H_topo.cp),1),...
    '.k','LineWidth',4,'MarkerSize',10);
