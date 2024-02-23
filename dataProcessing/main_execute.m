function main_execute(params, subj)
% analyzes the iEEG experiments for a single subject 'subj' in subject list
% runs analyses based on configuration

% (c) Jiri, 2016 - ...
% based on: whitenoise/mainEcog.m (c) Jiri, 2011 - 2015

%% subject tag
params.storage.subjTag = params.storage.subjList{subj,1};
disp(['Processing job: subj = ' params.storage.subjTag ', outDir = ' params.storage.outName ]);
disp([' - results dir = ' params.storage.dir_results]);

%% load data from D-struct & H-struct
params = load_main(params);
if isempty(params.storage.analyzedSess)
    display(['WARNING: paradigm type = ' params.paradigm.usedParadigm ': No sesssions found for subject: ' params.storage.subjTag]);
    return;
end

%% filtering of predictor & response datasets
loadAndFilter(params, params.predictor);
loadAndFilter(params, params.response);

%% rejection of bad data
rejectionAnalysis(params, params.predictor);
rejectionAnalysis(params, params.response);

%% check the data quality
signalQualityAnalysis(params);

%% trajectory analysis
%trajQualityAnalysis(params);

%% decoding / prediction based on multiple linear regression
decoding_linRegression(params);

%% tuning analysis
tuningAnalysis(params);

%% mutual information
mutualInfoAnalysis(params);

%% triggering analysis
triggerAnalysis(params);

%% save params & clean-up HDD (should be the last function)
clearStorage(params);


