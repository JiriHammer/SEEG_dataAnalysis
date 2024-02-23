function covStruct = getRespCovariance(dataStructDir, cacheFile, signalType, chPos)
% loads corresponding response data
% computes their autocovariance

% (c) Jiri, Feb12

%% load params from cache file -> used sessions
assert(exist(cacheFile,'file') == 2);
load(cacheFile, 'params');
assert(exist('params','var') == 1);

%% load response data
dataStructFile = [dataStructDir filesep 'expData_srate' num2str(params.amp.srate) '.mat'];
assert(exist(dataStructFile,'file') == 2);

load(dataStructFile, 'S');
assert(exist('S', 'var') == 1);
assert(size(S,1) == length(params.storage.sessionCacheFiles));
assert(S{1}.srate == params.amp.srate);

resp = cell(length(params.simulation.trainSession),1);
for sess = 1:length(params.simulation.trainSession)
    
    thisSess = params.simulation.trainSession(sess);
    assert(isfield(S{thisSess}, signalType));
    eval(['resp{sess} = S{thisSess}.' signalType ';' ]);
end
    
%% select channel
selCh = getDecodingChannels(params.connectionTable, signalType, chPos);

%% compute mean autocovariance
maxTime = 10;  % in [s]
timeAxis = [-maxTime:1/params.amp.srate:maxTime];
meanRespCov = [];
for sess = 1:size(resp,1)
    respCov = xcov(resp{sess}(:,selCh), maxTime*params.amp.srate, 'coef');
    meanRespCov= cat(2, meanRespCov, respCov);
    %plot(timeAxis, chCov, palette{sess});
end   
meanRespCov = mean(meanRespCov, 2);

% attach as handle to statData
covStruct.respCov = meanRespCov;
covStruct.timeCov = timeAxis;

%% check if scaling alters the shape? - yes, do not scale!
% figure;
% hold on;
% plot(covStruct.timeCov, covStruct.respCov, 'b');
% plot(covStruct.timeCov, 0.5*covStruct.respCov, 'r');