function plotEC_spectrum2axes(EC, freqAxis, timeAxis, cLims, clz, EC_info)
% plots effective connecivity (EC) data as imagesc to axes

% (c) Jiri, Nov22

%% defualt settings
if ~isfield(EC_info, 'plot_xlabel'), EC_info.plot_xlabel = true; end
if ~isfield(EC_info, 'plot_ylabel'), EC_info.plot_ylabel = true; end

%% >>> plot: imagesc <<<
h = imagesc(timeAxis, freqAxis, EC, cLims);
% h = imagesc(timeAxis, freqAxis, EC);
colormap(gca,brewermap(256,'*RdBu')); 
axis tight;
colorbar;

%% paradigm-specific times
if isfield(EC_info.info, 'paraTimes')
    xLims = get(gca, 'xlim');
    yLims = get(gca, 'ylim');
    for t = EC_info.info.paraTimes(clz,:)
        if ~isnan(t)
            if t > xLims(1) && t < xLims(2)
                plot([t, t], yLims, 'Color',EC_info.info.colors(clz,:), 'LineWidth',2);
            end
        end
    end
end

%% title, labels, ticks
% title({[EC_info.src_groupName '(' num2str(EC_info.src_Nch) ')']; 
%        [EC_info.trg_groupName '(' num2str(EC_info.trg_Nch) ')']; ...
%        [EC_info.methodName ': clz = ' EC_info.info.clzNames{clz,1}]});
title(EC_info.str_title);
if EC_info.plot_xlabel
    xlabel('time (s)');
else
    set(gca, 'XTickLabel',{''});
end
if EC_info.plot_ylabel
    ylabel('freq (Hz)');
else
    set(gca, 'YTickLabel',{''});
end

% if EC_info.thisPlot > (EC_info.nRows-1)*EC_info.nCols   % ~ last row   
%     xlabel('time (s)');
% else
%     set(gca, 'XTick',[]);
% end    

% if mod(EC_info.thisPlot,EC_info.nCols) == 1     % ~ first column
%     ylabel('freq (Hz)');
% else
%     set(gca, 'YTick',[]);
% end

xticks = [-2:2];
for t = xticks
    plot([t t], ylim, '--k');   % time ticks
end
yticks = [10, 30, 50, 100, 120];
for t = yticks
    plot(xlim, [t t], '--k');   % freq ticks
end
