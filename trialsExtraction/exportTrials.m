function exportTrials(params, trials2export)
% exports trials from trigger analysis
% exports selected channels (selCh_H_resp) & header structure H

% (c) Jiri, Mar21 - for bachelors at CVUT

%% trials
trials = struct;
trials.data = trials2export.data;
trials.time = trials2export.time;
trials.labels = trials2export.labels;
trials.label_names_values = trials2export.clzNames;
trials.info = trials2export.info;

%% update H.channels in header for Yeo7 atlas (saves H to cacheFile)
load(params.storage.cacheFile, 'H', 'selCh_H_resp');
updateHeaders_atlasInfo(params, selCh_H_resp, H, params.storage.cacheFile);

%% header
clear H selCh_H_resp
load(params.storage.cacheFile, 'H', 'selCh_H_resp', 'srate');
trials.samplingFrequency = srate;
channels = H.channels;
selected_channels = selCh_H_resp;

%% save
outdir = params.saveResults.exportTriggering_trials;
if ~exist(outdir, 'dir')
    mkdir(outdir);
end 
disp('saving trials ...');
fileName = [outdir filesep params.storage.subjTag '_trials.mat'];
save(fileName, 'trials', 'channels', 'selected_channels', 'params', '-v7.3');
disp(['saved trials into: ' fileName]);



