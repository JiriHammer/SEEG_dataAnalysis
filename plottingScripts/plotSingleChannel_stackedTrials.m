function plotSingleChannel_stackedTrials(params, spectra, spectralData)
% plot single stracked single trials for all freq. bands (ch-by-ch)
% spectra = 
%         data: [129×85×57×188 double] = 4D: [freq x time x chnl x trials]
%         freq: [1×129 double]
%         time: [85×1 double]
%       labels: [1×188 double]
%     clzNames: {2×2 cell}
%     rejected: [129×85×57×188 double]

% (c) Jiri, Jun20

%% color limits
cLims = 2*spectralData.info.cLims;

%% load H & selected channels
load(params.storage.cacheFile, 'H', 'selCh_H_resp');
selCh_H = selCh_H_resp;
assert(size(selCh_H,2) == size(spectra.data,3));

%% load RT

%% for each channel:
n = 1;
for ch = 1:size(spectra.data,3)
    
    %% compute freq.bands activations for all trials
    ch_vals = [];
    nFreq = size(params.triggering.freqBands,1);
    spectralData.info.chNames = cell(1, nFreq);
    for fb = 1:nFreq
        
        % select frequencies
        freqBand = params.triggering.freqBands{fb,2};
        freqBandName = params.triggering.freqBands{fb,1};
        spectralData.info.chNames{fb} = freqBandName;
        i_fr = closestval(spectra.freq,freqBand(1)):closestval(spectra.freq,freqBand(2));
        
        % select time
        %selTime = params.plot_triggering.time2plot;
        selTime = [-100, 100];      % select all samples to see also baseline
        i_t = closestval(spectra.time,selTime(1)):closestval(spectra.time,selTime(2));
        
        % sort trials: by RTs
        
        ch_vals = cat(3, ch_vals, squeeze(mean(spectra.data(i_fr,i_t,ch,:),1))');          % 3D = trials x time x freq
    end
    
    %% sort trials: by class label
    tmp_vals = nan(size(ch_vals));
    k = 0;
    spectralData.info.horLines = [];
    for clz = 1:size(spectra.clzNames,1)
        i_clz = find(spectra.labels == spectra.clzNames{clz,2});
        tmp_vals(k+[1:size(i_clz,2)],:,:) = ch_vals(i_clz,:,:);
        k = k + size(i_clz,2);
        spectralData.info.horLines = cat(2, spectralData.info.horLines, k); 
    end
    assert(isempty(find(isnan(tmp_vals),1)));
    
    %% modify 'spectralData'
    spectralData.cVals = tmp_vals;
    spectralData.yVals = [1:size(ch_vals,1)]';
    spectralData.xVals = spectra.time(i_t);
    spectralData.info.ylabel = 'trials';
    spectralData.info.figName = ['ch_' num2str(ch)];
    spectralData.info.outDir = [params.storage.outputDir filesep 'stackedTrials'];
    spectralData.info.marg_h = [0.05 0.05];
    spectralData.info.marg_w = [0.05 0.05];
    spectralData.info.gap = [0.050, 0.050];
    spectralData.info.verLines = [params.triggering.baseline, params.plot_triggering.time2plot];
    
    % text: chnl info
    thisCh = selCh_H(n,ch);                 % chnl index in H  
    plotInfo.thisCh = thisCh;
    plotInfo.channel = H.channels(thisCh);          
    %plotInfo.txt_pos = [0.016, 0.98-0.03*(n-1)];
    spectralData.info.text = textOnFig_chnlInfo(params, plotInfo);    
        
    %% plot stacked trials
    spectralData.info.cLims = cLims;
    plotSpectra(params, spectralData);
    
end
