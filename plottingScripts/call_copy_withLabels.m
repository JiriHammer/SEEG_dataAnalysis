function call_copy_withLabels(xstring, ystring)

a = gca;
f = gcf;
cm=get(f,'colormap'); 
cf=get(f,'color');
% mousepos = get(gcf,'CurrentPoint');
set(f,'units','pixel')
set(a,'units','pixel')
posA = get(a,'Position');
posF = get(f,'Position');
set(f,'units','normalized')
set(a,'units','normalized')
posAn = get(a,'Position');
posFn = get(f,'Position');

 set(0, 'DefaultFigureMenu', 'figure');  

% fsizen = [0.2917, 0.3889]; % Default matlab normalized
fsize = [560, 420]; % Default matlab pixel
coord = posF(1:2)+posA(1:2);
% coordn = posAn(1:2).*posFn(3:4);

MainScreensize = get(0, 'ScreenSize');
MonitorPos = get(0,'MonitorPositions');
if coord(1) <=  MainScreensize(3)/2;        %subplot was on the left side of the main screen
    coord(1) = MainScreensize(3)/2 - fsize(1)*1.5;
    coord(2) = MainScreensize(4)/2 - fsize(2)/2;
    
elseif coord(1) <=  MainScreensize(3) && coord(1) >  MainScreensize(3)/2  %subplot was on the right side of the main screen
    
    coord(1) = MainScreensize(3) - fsize(1)*1.5;
    coord(2) = MainScreensize(4)/2 - fsize(2)/2;
    
elseif coord(1) >  MainScreensize(3) ||  coord(1) < 0   %subplot was on the second screen
    
    coord(1) = MonitorPos(1,3) - fsize(1)*1.5;
    coord(2) = MonitorPos(1,4)/2 - fsize(2)/2;
    
end

%prevent figures from overlapping over the screen edge
% coord(coordn+fsizen > 1) = coord(coordn+fsizen > 1) -fsize(1);

% f = figure('units','normalized','Position',[posA(1:2).*posF(3:4), 0.2917, 0.3889]);
% copyobj(a,f);
copyobj(gca,figure('units','pixel','Position',[coord, fsize]))
set(gca,'units','normalized','Position',[0.13 0.11 0.775 0.815]);

set(gca,'xcolor', 'black');
xlhand = get(gca,'xlabel');
set(xlhand,'string','X','fontsize',9);   

set(gcf,'KeyPressFcn',@react);
set(gca,'ButtonDownFcn','');
h_child=get(gca,'Children');
set(h_child,'ButtonDownFcn','');
set(gcf,'colormap',cm);
set(gcf,'color',cf);

if strcmp(get(h_child(1),'Type'),'surface');
    colorbar;
end;

% ylim('auto');

if nargin > 1
 xlabel(xstring,'FontSize',14,'FontName','Lucida Sans');
 ylabel(ystring,'FontSize',14,'FontName','Lucida Sans');
end
%axis on;
