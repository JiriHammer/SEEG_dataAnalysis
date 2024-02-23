function plotSingleTrials(params, trials)
% plots single trials

% (c) Jiri, Nov16

%% aux. info
clear trialsData;
trialsData = struct;
trialsData.info.verLines = [-2:2];
trialsData.info.horLines = 0;
trialsData.info.yLims = getYLims(trials.data);
% trialsData.info.outDir = [params.storage.outputDir filesep 'singleTrialPlots'];

%% storage info
if isfield(trials, 'outDir')
    trialsData.info.outDir = trials.outDir;
else
    trialsData.info.outDir = [params.storage.outputDir filesep 'singleTrialPlots'];  % default
end

%% channel names
params.nCh = size(trials.data,2);
trialsData.info.chNames = getChannelNames(params, lower(params.response.signalType));

%% selected time
i_t = closestval(trials.time,params.plot_triggering.time2plot(1)):closestval(trials.time,params.plot_triggering.time2plot(2));

%% outliers of the disctribution
outliers = getOutliers(trials.data(i_t,:,:));

%% plot
for tr = 1:size(trials.data,3)
    trialsData.yVals = trials.data(i_t,:,tr);
    trialsData.xVals = trials.time(i_t);
    trialsData.info.outliers = outliers(:,tr);
    trialsData.info.rejected = trials.rejected(i_t,:,tr);
    trialsData.info.figName = ['trial_' num2str(tr) '_data'];
    [tf, i_c] = ismember(trials.labels(tr),cell2mat(trials.clzNames(:,2)));
    assert(tf);          
    trialsData.info.text = [trials.text ' TRIAL: ' num2str(tr) ',  subject: ' params.storage.subjTag ...
        ', triggered: ' params.triggering.cutPoint ...
        ', class: ' trials.clzNames{i_c,1}];        

    plotTrials(params, trialsData);
end
