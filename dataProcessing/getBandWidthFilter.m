function [loF, hiF, n] = getBandWidthFilter(F)
% look-up table. Based on carrier frequency F [Hz], returns freq. band boundaries
% and order n for butterworth filter

% (c) Jiri, Jul12

assert(F >= 0);

if F < 1.5
    loF = 2/3*F;
    hiF = 4/3*F;
    n = 2;
elseif F >=1.5 && F < 7    % 2 Hz bin
    loF = F-1;
    hiF = F+1;
    n = 3;  
elseif F >=7 && F < 15     % 4 Hz bin
    loF = F-2;
    hiF = F+2;
    n = 4;      
elseif F >=15 && F < 30     % 6 Hz bin
    loF = F-3;
    hiF = F+3;
    n = 4;
else                        % 10 Hz bin
    loF = F-5.0;
    hiF = F+5.0;
    n = 5;
end