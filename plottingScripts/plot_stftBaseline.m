function plot_stftBaseline(params, freqAxis, timeAxis, B_avg, B_std, figName)
% plots STFT baseline for all channels
% input: B_avg, B_std = 2D: [freq x channels]
% prepapres struct 'trialsData'
% plotting: plotTrials.m

% (c) Jiri, Aug23

if nargin < 6
    figName = 'baseline_RS';
end

%% data
trialsData = struct;
trialsData.xVals = freqAxis;
trialsData.yVals = B_avg;
trialsData.yErrs = B_std;

%% channel names
params.nCh = size(trialsData.yVals,2);
trialsData.info.chNames = getChannelNames(params, lower(params.response.signalType), 'resp');

%% other info
trialsData.info.xlabel = 'freq [Hz]';
trialsData.info.verLines = [10, 30, 50, 100];
trialsData.info.text = ['subj = ' params.storage.subjTag ', baseline from resting state, duration = ' num2str(timeAxis(end)/60) ' min'];
trialsData.info.outDir = [params.storage.outputDir filesep 'baseline_RS'];
trialsData.info.figName = figName;

%% plot
plotTrials(params, trialsData);
