function [xf,p] = fsgwindow(x,ord,fw,lag,step,positions)
% function sgwindow - fast windowed Savitzky-Golay filtering
% 
% syntax: xf = fsgwindow(x,ord,fw,lag)
% syntax: [xf,p] = fsgwindow(x,ord,fw,lag,step,positions)
%
% Input parameters:
% x: Input data, will be filtered column-wise.
% ord: Filter order.
% fw: For scalar input: total filter width. Kernel will be (a)symmetricly
%     centered around 'lag'.
%     If 'fw' is a 2-element vector, it is giving [left_with, right_width].
% lag: Samples are picked 'lag' samples away from the end of each window.
%      (i.e. going back to the past from the current position)
% step: Window will be advanced in steps of step.
% positions: If given, only picks values from windows at 'positions',
%            alligned to the end of the window.
%
% Output parameters:
% xf: Filtered x.
% p: Window positions. Equals 'positions', when given.

if ndims(x)==2 && min(size(x))==1
    x= x(:);
end

if nargin<4
    lag= 0;
end

if nargin<5
    step= 1;
end

if nargin<6
    p= lag+(1:step:size(x,1));
else
    p= positions+lag;
    if min(p)<1 || max(p)>size(x,1)
        error('position values out of range.')
    end
end

if length(fw)<2
    flength= fw;
    if lag<(fw+1)/2
        b= SavitzkyGolay([-lag,(fw-lag-1)],'PolynomialOrder',ord);
    else
        % this case is untested (but rarely used)
        b= SavitzkyGolay([-ceil((fw-2)/2),floor(fw/2)],'PolynomialOrder',ord);
        p= p- (lag-floor((fw-1)/2));
    end
else
    flength= -fw(1)+fw(2)+1;
    b= SavitzkyGolay(fw(1:2),'PolynomialOrder',ord);
end

% --- original ---
% if (size(x,1)<1000)
%     temp= filter(b,1,[x;zeros(flength,size(x,2))]);
% else
%     temp= fftfilt(fliplr(b),[x;zeros(flength,size(x,2))]);
% end
% -------------

% ---- test ----
temp= filter(b,1,[x;zeros(flength,size(x,2))]);
% --------------

xf= temp(p,:);
