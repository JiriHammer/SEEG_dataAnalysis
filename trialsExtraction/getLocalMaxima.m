function y = getLocalMaxima(params, x)
% extracts local maxima from vector x
% returns y = 2D: samples x values of the local maxima

% (c) Jiri, Jun12

%% define sampling rate
load(params.storage.cacheFile, 'sampleRateAfterFiltering');
assert(exist('sampleRateAfterFiltering','var') == 1);
srate = sampleRateAfterFiltering;

%% determine smoothing filter after derivation
Wn = 5/(srate/2);                          % normalized bandpass frequencies
n = 3;                                      % butterworth order
[b,a] = butter(n, Wn, 'low');               % returns polynoms of Butterw. filter

%% estimate derivatives & smooth
i_array = 1:size(x,1);
der1 = estFirstDerivative(x, 1/srate);
der1 = filtfilt(b, a, der1);
der2 = estFirstDerivative(der1, 1/srate);
der2 = filtfilt(b, a, der2);

%% de-select first & last seconds from triggering (filter artifacts)
i_tCut1 = abs(round(params.triggering.time2cut(1)*srate));      % in [samples], cut edges of data
i_tCut2 = abs(round(params.triggering.time2cut(2)*srate));      % in [samples], cut edges of data

i_array = i_array((i_tCut1:end-i_tCut2));                       % holds the indices
der1 = der1(i_tCut1:end-i_tCut2);
der2 = der2(i_tCut1:end-i_tCut2);

%% detect local maxima
i_der1 = find(abs(der1) < 8e-2);                    % 1st derivation = 0, indices in der1
i_der1 = i_der1(diff(cat(1,-666,i_der1)) > 1);      % take only first samples that satisfied the condition
i_der2 = find(der2 < 0);                            % 2nd derivation < 0, indices in der2
i_locMax = intersect(i_array(i_der1), i_array(i_der2));

%% return output 'y'
y = nan(length(i_locMax), 2);
for n = 1:length(i_locMax)
    y(n,1) = i_locMax(n);
    y(n,2) = x(i_locMax(n));
end

%% debug: detection, figure, x-axis in samples
% t = [1:size(x,1)]./srate;
% f = figure('visible', 'on', 'Position', [1, 1, 1920, 1200]);
% set(f,'PaperUnits','inches','PaperPosition',[0 0 19 12]);
% hold on;
% plot(t, x, 'r');
% plot(i_array./srate, der1./max(der1), 'm');
% plot(i_array./srate, der2./max(der2), 'c');
% plot([min(t), max(t)], [0, 0], '--k');
% for n = 1:length(i_locMax)
%     plot([i_locMax(n), i_locMax(n)]./srate, get(gca,'ylim'), '--g');        % triggers
% end    
% legend('origSignal', 'der 1', 'der 2');
% plot([min(t), max(t)], [0, 0], '--k');
% xlabel('time [s]');
% close(f);       
   
%% debug: detection, figure, x-axis in samples
%     f = figure('visible', 'on', 'Position', [1, 1, 1920, 1200]);
%     set(f,'PaperUnits','inches','PaperPosition',[0 0 19 12]);
%     hold on;
%     plot(x, 'r');
%     plot(der1./std(der1,0,1), 'm');
%     plot(der2./std(der2,0,1), 'c');
%     plot([min(t), max(t)], [0, 0], '--k');
%     for n = 1:length(i_locMax)
%         plot([i_locMax(n), i_locMax(n)], get(gca,'ylim'), '--g');        % triggers
%     end    
%     legend('origSignal', 'der 1', 'der 2');  