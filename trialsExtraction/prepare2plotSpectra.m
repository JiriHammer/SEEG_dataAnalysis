function spectralData = prepare2plotSpectra(params, spectra, clz)
% converts 'spectra' struct to 'spectralData' struct for plotting
% plotting function: plotSpectra(params, spectralData)

% (c) Jiri, Sep22


%% aux. info
spectralData = struct;
spectralData.info.verLines = params.plot_triggering.verLines;
spectralData.info.horLines = [13 50 100 120];
spectralData.info.xlabel = 'time (s)';
spectralData.info.ylabel = 'freq (Hz)';
%spectralData.info.outDir = [params.storage.outputDir filesep 'classMeansPlots'];
spectralData.info.outDir = [params.storage.dir_results filesep params.storage.outName];
    
%% color limits
%     spectralData.info.cLims = [-4, 4];     % manual !
% spectralData.info.cLims = [-0.5, 0.5];     % z-score, manual !
% spectralData.info.cLims = [-0.1, 0.1];     % z-score, manual !
c_tmp = spectra.data;                   % automatic ...
% c_tmp(spectra.rejected == 1) = NaN;
spectralData.info.cLims = getYLims(nanmean(c_tmp,4));
clear c_tmp;
    
%% channel names
params.nCh = size(spectra.data,3);
spectralData.info.chNames = getChannelNames(params, 'ieeg');  
    
%% time selection from plot_triggering.time2plot
i_t = closestval(spectra.time,params.plot_triggering.time2plot(1)):closestval(spectra.time,params.plot_triggering.time2plot(2));

%% spectralData vals (class means)
spectralData.xVals = spectra.time(i_t);
spectralData.yVals = spectra.freq';

%% class labels & text on figure
clzLabels = unique(spectra.labels); 
[tf, i_c] = ismember(clzLabels(clz),cell2mat(spectra.clzNames(:,2)));
assert(tf);
spectralData.info.figName = ['spectra_' params.storage.subjTag '_clz' num2str(i_c)];      % use real class names
i_clz = find(spectra.labels == clzLabels(clz));                     % selected trials
spectralData.info.text = ['subject: ' params.storage.subjTag ...
    ', triggered: ' params.triggering.cutPoint ...
    ', class: ' spectra.clzNames{i_c,1} '(' num2str(length(i_clz)) ')'];     

%% not-rejected indices -> vals
if isfield(spectra, 'rejected')
    inds = ~spectra.rejected(:,i_t,:,:);                                           % not rejected indices = 1, rejected indices = 0
    i_otherTrials = setdiff(1:size(spectra.data,4),i_clz);              % all other trials
    inds(:,:,:,i_otherTrials) = false;                                  % not selected trials are set to 0 (i.e. rejected)
    vals = spectra.data(:,i_t,:,:);
    vals(inds == 0) = NaN;       
else
    vals = spectra.data(:,i_t,:,i_clz);
end

%% spectralData vals = color-coded values (class means)
spectralData.cVals = nanmean(vals,4);       % 3D: [freq x time x chnls]
%spectralData.cVals = nanmedian(vals,4);      % 3D: [freq x time x chnls]

%% paradigm times
if ~isempty(spectra.info.paraTimes)
    spectralData.info.paraTimes = spectra.info.paraTimes(clz,:);
else
    spectralData.info.paraTimes = [];
end
spectralData.info.clr_clz = spectra.info.colors(clz,:);
    
   

