%% fig: subtighplot - examplary figure
% useful for copy-pasting

% (c) Jiri, Jul23

%% examplary figure
% f = fig_make;
f = fig_make('fig_W',20, 'fig_H',20); 
nRows = 6;
nCols = 6;
nPlot = 1;

marg_h = [0.10 0.05];   % margin from [bottom, top]
marg_w = [0.10 0.05];   % margin from [L, R] side
gap = [0.05, 0.01];     % between axes from [top, side]

% dummy data
t = -2:0.01:2;
w = 1:nRows*nCols;
txt_labels = {'axes 1', 'axes 2', 'axes 3', 'axes 4', 'axes 5', 'axes 6'};

fontSize = 10;
lw = 2;

for row = 1:nRows
    for col = 1:nCols
        % axes
%         subplot(nRows,nCols,nPlot);
        ax = subtightplot(nRows, nCols, nPlot, gap, marg_h, marg_w); 
        hold on;
        
        % plot
        plot(t, sin(w(nPlot)*t), 'b', 'LineWidth',lw);
        
        % axis props
        axis tight;
        xlim([-1.3,1.3]);
        set(ax, 'XTick',[-1 -0.5 0 0.5 1]);
        set(ax, 'XTickLabel',{'-1', '', '0', '', '1'});        
        set(ax, 'XTickLabelRotation',90);        
        
        set(gca, 'ylim',[-0.55, 1.1]);
        plot([0 0], ylim, '--k');
        plot(xlim, [0 0], '--k');
        
        % xlabel: only bottom row 
        if nPlot > nCols*(nRows-1)
            xlabel('time lag (s)', 'FontSize',fontSize);
        else
            set(gca,'XTickLabel', {''});
        end
        
        % ylabel: only L column 
        if mod(nPlot,nCols) == 1
            ylabel('signal', 'FontSize',fontSize);
        else
            set(gca,'YTickLabel', {''});
        end
        title(['sine at ' num2str(w(nPlot)) ' Hz'], 'FontSize',fontSize);
        
        % text (traj vars names)
        xLims = get(gca, 'xlim');
        yLims = get(gca, 'ylim');        
        if mod(nPlot,nCols) == 1
            h_txt = text(xLims(1)-1.5, yLims(1)+0.2*diff(yLims), txt_labels{row}, 'FontSize',10, 'FontWeight','bold', 'Color','k');   
            set(h_txt,'Rotation',90);
        end
        if nPlot<=nCols
            h_txt = text(xLims(1)+0.2*diff(xLims), yLims(2)+0.4*diff(yLims), txt_labels{col}, 'FontSize',10, 'FontWeight','bold', 'Color','k');   
        end        
        set(gca, 'FontSize',fontSize);
        box on;
        grid on;

        nPlot = nPlot+1;
    end
end

% text on figure
% tx = axes('visible','off', 'position',[0 0 1 1]);
% mytitle = 'XCORR: ALL';
% mytitle = strrep(mytitle, '_','\_');
% text(0.016, 0.97, mytitle, 'fonts', 16, 'fontw', 'bold');    

%% save
figName = 'test';
figDir = 'F:\dox\ms4_distance\figs\tp_tuning\v15_notch50\xcorr_trajs_v2';
fig_save(f, figName, figDir, 'res',600);
% close(f);             
