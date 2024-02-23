%% filter's frequency response
% compares different filter responses (sgolay, butterworth lowpass)

%% settings
srate = 1200;

%% sgolay
order=2;
winLen=1.5*srate + 1;  % in [s]
B=sgolay(order,winLen); 
sgB = B(round(winLen/2),:);
sgA = 1;

%% butterworth
hiF = 1;  % in [Hz]
Wn = hiF/(srate/2);           % normalized bandpass frequencies
n = 3;                          % butterworth order
[bwB,bwA] = butter(n, Wn, 'low');          % returns polynoms of Butterw. filter

%% Filter Visualization Tool
fvtool(sgB, sgA, bwB, bwA)

%% gaussian smoothing kernel
srate = 500;
tb = 0.125;
nb = round(tb*srate);

%wvtool(gausswin(nb));

b = gausswin(nb);
a = 1;
fhandle = fvtool(b,a, 'Analysis','freq');
set(fhandle, 'Fs',srate, 'NormalizedFrequency','off');

%% sgolay
srate = 500;
tb = 0.5;
nb = round(tb*srate)+1;
b = sgolay(2, nb);
a = 1;
fhandle = fvtool(b(:,round(end/2)),a, 'Analysis','freq');
set(fhandle, 'Fs',srate, 'NormalizedFrequency','off');

%% rad
srate = 500;
%tb = 0.333;
tb = 0.100;
nb = round(tb*srate)+1;
b = ones(nb,1)*1/nb;
%b(round(end/2)) = b(round(end/2))+1;
a = 1;
fhandle = fvtool(b,a, 'Analysis','freq');
set(fhandle, 'Fs',srate, 'NormalizedFrequency','off');

