function plotEC_plotband2axes(y_avg, y_sem, timeAxis, EC_info)
% plots EC to axes as plotband

% (c) Jiri, Nov22


%% >>> plot G_src -> G_trg <<<
for clz = 1:size(y_avg,2)
    plotband(timeAxis, y_avg(:,clz), y_sem(:,clz), EC_info.info.colors(clz,:));
end
axis tight;

%% paradigm-specific times
if isfield(EC_info.info, 'paraTimes')
    xLims = get(gca, 'xlim');
    yLims = get(gca, 'ylim');    
    for clz = 1:size(y_avg,2)
        clr = EC_info.info.colors(clz,:);
        for t = EC_info.info.paraTimes(clz,:)
            if ~isnan(t)
                if t > xLims(1) && t < xLims(2)
                    plot([t, t], yLims, 'Color',clr, 'LineWidth',1);
                end
            end
        end
    end
end

%% plot significances
if isfield(EC_info, 'hVals')
    plot2axes_signif_filled(timeAxis, EC_info.hVals, EC_info.pVals);
end

%% title, labels, ticks
% title({[EC_info.src_groupName '(' num2str(EC_info.src_Nch) ') -> ' EC_info.trg_groupName '(' num2str(EC_info.trg_Nch) ')']; ...
%        [EC_info.methodName ' (' EC_info.fbName ') : clz = ' EC_info.info.clzNames{clz,1}]});
title({[EC_info.src_groupName '(' num2str(EC_info.src_Nch) ')']; 
       [EC_info.trg_groupName '(' num2str(EC_info.trg_Nch) ')']});   
if EC_info.thisPlot > (EC_info.nRows-1)*EC_info.nCols   % ~ last row   
    xlabel('time (s)');
else
    set(gca, 'XTick',[]);
end    
if mod(EC_info.thisPlot,EC_info.nCols) == 1     % ~ first column
    ylabel(EC_info.fbName);
else
%     set(gca, 'YTick',[]);
end

box on;
grid on;
