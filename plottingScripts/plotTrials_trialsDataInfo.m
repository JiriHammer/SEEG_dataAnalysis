function trialsData = plotTrials_trialsDataInfo(data_info)
% returns auxiliary info in trialsData.info

% (c) Jiri, Oct20

%% decoding: 'trialsData' auxiliary info
clear trialsData;
trialsData = struct;
trialsData.info.verLines = [-1, 0, 1];

trialsData.info.xlabel = 'time offset (s)';
trialsData.info.ylabel = 'DA (CC)';

% y-limits for each channel (if commented out, the program estimates common y-lims for all channels)
trialsData.info.yLims = [-0.25, 0.45];      % defined common y-limits for all channels
trialsData.info.horLines = [-0.2:0.2:0.4];
trialsData.info.outliers = zeros(data_info.nCh,1);

trialsData.info.outDir = data_info.figDir;
trialsData.info.figName = [data_info.figName '_' data_info.subjTag];
trialsData.info.text = ['Decoding. Subject: ' data_info.subjTag ', ' data_info.text_nt];     

% channel names
params = data_info.params;
params.nCh = data_info.nCh; 
trialsData.info.chNames = getChannelNames(params,params.predictor.signalType, 'pred');
