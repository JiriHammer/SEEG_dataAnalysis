function sLevel = getSignificanceLevel(freqTag)
% returns significance level based on empirical distribution
% sLevel = [sLo, sHi] ~ [min,max] of the empirical distribution

% (c) Jiri, Jul13

%% settings
fileName = '/export/homer/hammer/Documents/paper3_pauseResponse/data/pauseTriggering_tStep5/surrog_sLevel/sLevel_m1.mat';

%% load empirical data distribution
assert(exist(fileName,'file') == 2);
load(fileName, 'surrPauseTime');          % format: surrPauseTime{freq}.yAvgs(pause,time,surr)

%% find the right tag
i_tag = [];
for freq = 1:size(surrPauseTime,1)
    if strcmp(surrPauseTime{freq}.freqTag{freq}, freqTag)
        i_tag = freq;
    end
end
assert(~isempty(i_tag));
dd = surrPauseTime{i_tag}.yAvgs(:);     % surrogate data distribution

%% significance level
sLevel = [min(dd), max(dd)];
