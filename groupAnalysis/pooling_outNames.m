function pooling_outNames(params, dirsList, figDir)
% pools across results in outNames (different analysis settings)
% result = from triggerAnalysis -> trialData
% plots one figure per channel for all subjects

% (c) Jiri, May22

%% settings: activations
freqBands = cat(1, {'erp', [nan, nan]}, params.triggering.freqBands);   % add ERP to freq. bands

%% pool over subj

for subj = 1:size(params.storage.subjList,1)
    subjTag = params.storage.subjList{subj,1};
    params.storage.subjTag = subjTag;
    plotInfo.subjTag = subjTag;

    %% cache file (based on params)
    cacheFile = [dirsList{1,1} filesep subjTag filesep 'cacheFile.mat'];
    assert(exist(cacheFile,'file') == 2);
    plotInfo.cacheFile = cacheFile;
    
    %% load structure H & selected channels
    clear H selCh_H_resp;
    load(cacheFile, 'H', 'selCh_H_resp');    
    selCh_H = selCh_H_resp;

    %% anatomical assignements from isarg_atlas
    [ass_isargAtlas, H] = updateHeaders_atlasInfo(params, selCh_H_resp, H, cacheFile);
    
    %% load trialsData over outNames & freqBands -> D
    D = cell(size(freqBands,1),size(dirsList,1));
    for dirs = 1:size(dirsList,1)
        
        %% cache file
        cacheFile = [dirsList{dirs,1} filesep subjTag filesep 'cacheFile.mat'];
        assert(exist(cacheFile,'file') == 2);
        plotInfo.cacheFile = cacheFile;

        %% load freq. band activations
        for freq = 1:size(freqBands,1)
            varName = ['trialsData_' freqBands{freq,1}];  
            % try to load freq. band activation from cache file
            if ~ismember(varName, who('-file', cacheFile))   % if the variable 'varName' exists
                error(['WARNING: subject = ' subjTag ', variable = ' varName ' was not found in cacheFile: ' cacheFile]);
            else
                load(cacheFile, varName);
                assert(exist(varName,'var') == 1);
                eval(['trialsData = ' varName ';']); 
                clear(varName); 
                D{freq,dirs} = trialsData;
                clear trialsData;
            end
        end
    end
    
    %% figure: plot for each channel
    for ch = 1:size(selCh_H,2)
        f = fig_make;
        nCols = size(dirsList,1);
        nRows = size(freqBands,1);    
        nPlot = 1;
        
        plotInfo.fig = f;
        plotInfo.nRows = nRows;
        plotInfo.nCols = nCols;  
        % plotInfo.marg_h = marg_h;
        % plotInfo.marg_w = marg_h;
        % plotInfo.gap = gap+0.03;
        plotInfo.fontSize = 10;        
        
        % plot for each freq & outname
        for freq = 1:size(D,1)
            for dirs = 1:size(D,2)
                h_axes = subplot(nRows, nCols, nPlot);
%                 h_axes = subtightplot(nRows, nCols, v_subs(nPlot), 2*gap, 2*marg_h, 2*marg_w);
                hold on;
                plotInfo.chnlsAxes = h_axes;
                trialsData = D{freq,dirs};
                colors = trialsData.info.colors;
                trialsData.info.chNames{ch} = [];
                plotChannel2Axes(trialsData, ch, plotInfo);

                if nPlot > (nRows-1)*nCols
                    xlabel('time [s]');
                end
                if mod(nPlot,nCols) == 1
                    ylabel(freqBands{freq,1});
                end
                i_wordClasses = strfind(trialsData.info.text, ', classes: ');
                titleStr = trialsData.info.text(i_wordClasses+size(', classes: ',2):end);
                if nPlot <= nCols
                    i_word = strfind(dirsList{dirs,1}, '\');
                    outName = dirsList{dirs,1}(i_word(end)+1:end);
                    titleStr = {strrep(outName, '_','\_');, titleStr};
                end
                title(titleStr);
                nPlot = nPlot+1;
            end
        end   
        
        % text on top
%         tx = axes('visible','off', 'position',[0 0 1 1]);
%         mytitle = [];
%         mytitle = strrep(mytitle, '_','\_');
%         text(0.016, 0.97, mytitle, 'fontsize', 14, 'fontw', 'bold');
        plotInfo.txt_pos = [0.016, 0.97];
        plotInfo.thisCh = selCh_H(ch);
        plotInfo.channel = H.channels(selCh_H(ch));
        textOnFig_chnlInfo(params, plotInfo);    
        
        % save figure
        outDir = [figDir filesep subjTag];
        figname = ['ch_' num2str(ch)];
        fig_save(f, figname, outDir);
        close(f);            
        
    end     % of chnls
end     % of subj
    
    
  
    