function H=makeSignificanceBar_jiri(x,y,p,thr)
%makeSignificanceBar produces the bar and defines how many asterisks we get for a 
%given p-value & significance threshold 'thr'

% for example:
% x = [0.6923, 0.8462]; % x-axis location for the sigbar
% y = 0.9500;           % y-axis location for the sigbar
% p = 0.001;            % p-value

% (c) Jiri, Apr17, modified from Rob Campbell - CSHL 2013, see: sigstar.m (at MathWorks.com)

if nargin < 4
    thr = 0.01;
end

%% significance level
% if p<=1E-3
%     stars='***'; 
% elseif p<=1E-2
%     stars='**';
% elseif p<=0.05
%     stars='*';
% else isnan(p)
%     stars='n.s.';
%     %stars='';
% end

% if p<=thr(3)
%     stars='***'; 
% elseif p<=thr(2)
%     stars='**'; 
% elseif p<=thr(1)
%     stars='*'; 
% else
%     stars='n.s.';
%     %stars='';
% end

if p<=thr(1)
    stars='*'; 
else
    %stars='n.s.';
    stars='';
end

%% plot data
x=repmat(x,2,1);
y=repmat(y,4,1);

H(1)=plot(x(:),y,'-k','LineWidth',1.5,'Tag','sigstar_bar');

%starY=mean(y)+myRange(ylim)*offset;
starY=y(1)+0.02*diff(ylim);
H(2)=text(mean(x(:)),starY,stars,...
    'HorizontalAlignment','Center',...
    'FontSize', 16, ...
    'BackGroundColor','none',...
    'Tag','sigstar_stars');

% adjust ylims
Y=ylim;
if Y(2)<starY
    ylim([Y(1),starY+diff(Y)*0.05])
end


yd=diff(ylim)*0.01; %Ticks are 1% of the y axis range
y=get(H(1,1),'YData');
y(1)=y(1)-yd;
y(4)=y(4)-yd;   
set(H(1,1),'YData',y)

    
