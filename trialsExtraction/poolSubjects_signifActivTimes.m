function poolSubjects_signifActivTimes(params, AA_type, SA_params)
% plots histograms of times of significant activations in anatomical areas 
% activation time = when activations in given freq. band were significantly different
% !!! works for only 2 classes !!!
% pools together subject processed by a given job
% settings (subjects, directory names, ...) are stored in 'params'
% data stored in struct:
% trialsData = 
%           info: [1×1 struct]
%          xVals: [50×1 double]
%          yVals: [50×72×2 double]
%          yErrs: [50×72×2 double]
%          pVals: [50×72 double]
%          hVals: [50×72 double]
%     pVals_base: [50×72×2 double]
%     hVals_base: [50×72×2 double]
    
% (c) Jiri, Feb19

disp('Histogram of significant activations in anatomical areas ...');

%% default activation times
if nargin < 3
    SA_params = struct;
    SA_params.whichActivationTimes = 'first_significant';   
    SA_params.signifValsNames = {'pVals','hVals'}; % choices: {'pVals','hVals'} or {'pVals_base','hVals_base'}
end

%% default anatomy area
if nargin < 2 || isempty(AA_type)
    AA_type = 'neurologists_anatomy';   
end

%% define anatomy areas
AA_list = anatomicalAreas_getList(params, AA_type);

%% frequency bands
freqBands = cat(1, {'erp', [nan, nan]}, params.triggering.freqBands);   % add ERP to freq. bands

%% TO DO: min. number of channels to show in histograms
nChnls_cutoff = 2;

%% allocate time vector for binning: load trialsData_hiGamma (last of freq. band)
cacheFile = [params.storage.dir_results filesep params.storage.outName filesep params.storage.subjList{1} filesep 'cacheFile.mat'];
assert(exist(cacheFile,'file') == 2);
varName = ['trialsData_' freqBands{end,1}];     % last of freq. band
load(cacheFile, varName);
assert(exist(varName,'var') == 1);
eval(['trialsData = ' varName ';']); 
clear(varName);
  
%% time vector binning
if isfield(SA_params, 'histogramBinSize')
    tBins_stepSize = SA_params.histogramBinSize;
else
    tBins_stepSize = mean(diff(trialsData.xVals),1);    % default (estimated from data)
end
tBins_vector = trialsData.xVals(1)-tBins_stepSize/2:tBins_stepSize:trialsData.xVals(end)+tBins_stepSize/2;

nClz = size(trialsData.(SA_params.signifValsNames{1}),3);

%% output directory
outDir =[params.storage.dir_results filesep params.storage.outName filesep 'timesOfActivations_anatomAreas' filesep AA_type '_' SA_params.whichActivationTimes '_' SA_params.signifValsNames{1}];
if ~exist(outDir, 'dir')
    mkdir(outDir);
end          

%% plot frequency bands + ERP (over subjects and anatomy areas)
for freq = 1:size(freqBands,1)
    for c = 1:nClz
        varName = ['trialsData_' freqBands{freq,1}]; 
        t = cell(size(AA_list,1),2);            % times of activation: 2D = AA x 2 x clz, where 2 => bigger (or smaller)
        sbj = cell(size(AA_list,1),2);          % counts subjects
        chnls = zeros(size(AA_list,1),2);       % counts channels

        % go thru all subjects
        for subj = 1:size(params.storage.subjList,1)
            subjTag = params.storage.subjList{subj,1};

            % cache file
            cacheFile = [params.storage.dir_results filesep params.storage.outName filesep subjTag filesep 'cacheFile.mat'];
            assert(exist(cacheFile,'file') == 2);

            % load structure H & selected channels
            clear H selCh_H_resp;
            load(cacheFile, 'H', 'selCh_H_resp');
            selCh_H = selCh_H_resp;

            % try to load freq. band activation from cache file -> trialsData
            if ~ismember(varName, who('-file', cacheFile))   % if the variable 'varName' exists
                disp(['WARNING: subject = ' subjTag ', variable = ' varName ' was not found in cacheFile: ' cacheFile]);
            else
                % load trialsData
                clear(varName);
                clear trialsData;
                load(cacheFile, varName);
                assert(exist(varName,'var') == 1);
                eval(['trialsData = ' varName ';']); 
                clear(varName);

                % only 2 classes selected?
                if size(trialsData.yVals,3) ~= 2 && strcmp(SA_params.signifValsNames{1},'pVals')
                    disp('WARNING: p-values are computed only for 2 classes');
                    disp('WARNING: No significances were computed.');
                    disp(['WARNING: Number of selected classes = ' num2str(size(trialsData.yVals,3))]);
                    disp('WARNING: Check your params.triggering.classes selection. Aborting!');
                    return;
                end            
                assert(size(trialsData.yVals,2) == size(selCh_H_resp,2));
                assert(size(trialsData.yVals,3) == size(trialsData.info.clzNames,1));
                assert(isfield(trialsData, SA_params.signifValsNames{1}));   % ~ trialsData.pVals
                assert(isfield(trialsData, SA_params.signifValsNames{2}));   % ~ trialsData.hVals
                assert(~any(isnan(trialsData.(SA_params.signifValsNames{1})(:))));

                % sort used channels to AA
                for ch = 1:size(selCh_H,2)
                    chVals_h = trialsData.(SA_params.signifValsNames{2})(:,ch,c);        % ~ trialsData.hVals(:,ch,c);
                    chVals_p = trialsData.(SA_params.signifValsNames{1})(:,ch,c);        % ~ trialsData.pVals(:,ch,c);

                    % was significant activation?
                    if any(chVals_h)

                        % get anatomy area name of each channel
                        thisAA = [];
                        for n = 1:size(selCh_H,1)
                            thisCh = selCh_H(n,ch);                               % chnl index in H  
                            thisAA = cat(2, thisAA, anatomicalArea_getName(AA_type, AA_list, H.channels(thisCh)));
                        end

                        % was assigned?
                        wasAssigned = false;
                        for a = 1:size(thisAA,2)
                            [aa_found,i_aa] = ismember(lower(thisAA{a}), lower(AA_list(:,1)));
                            if aa_found  

                                % >>> get significant activation time <<<
                                SA_params.thisClz = c;
                                t_SA = getSignificantActivationTimes(chVals_h, chVals_p, trialsData.xVals, squeeze(trialsData.yVals(:,ch,:)), SA_params);
                                for clz = 1:size(t_SA,2)
                                    t{i_aa,clz} = cat(2, t{i_aa,clz}, t_SA{clz}');
                                    if ~isempty(t_SA{clz})
                                        chnls(i_aa,clz) = chnls(i_aa,clz) + 1;
                                        sbj{i_aa,clz} = cat(2, sbj{i_aa,clz}, subj);                
                                    end
                                end
                                wasAssigned = true;

                            end
                        end
                        if ~wasAssigned
                            for n = 1:size(selCh_H,1)
                                thisCh = selCh_H(n,ch);                               % chnl index in H 
                                disp([subjTag ', ch = ' H.channels(thisCh).name '(' num2str(thisCh) ') - not assigned: neurologist = ' H.channels(thisCh).neurologyLabel ', cytoarch = ' H.channels(thisCh).ass_cytoarchMap ', toolbox = ' H.channels(thisCh).ass_brainAtlas]);
                            end
                        end
                    end
                end
            end
        end

        %% find number of subplots
        nSubs = 0;
        for aa = 1:size(t,1)
            if chnls(aa,1) >= nChnls_cutoff || chnls(aa,2) >= nChnls_cutoff    % at least 'nChnls_cutoff' (e.g. 2) different channels contribute
                nSubs = nSubs+1;
            end
        end

        %% find max. bin count
        binCounts = [];
        f = figure;
        histVals = nan(size(tBins_vector,2)-1,size(t,1),2);   % binned values, 3D: time x AA x clz (=bigger/smaller)
        for aa = 1:size(t,1)
            for clz = 1:size(t,2)
                if size(t{aa,clz},2) >= nChnls_cutoff    % at least 'nChnls_cutoff' (e.g. 2) different channels contribute            
                    hh = histogram(t{aa,clz}, tBins_vector);
                    binCounts = cat(1, binCounts, max(hh.Values(:)));
                    histVals(:,aa,clz) = hh.Values;
                end
            end
        end
        close(f);

        %% -----------plot histogram distribution for each AA------------------
        % figure settings
        f = figure('units','normalized','outerposition',[0 0 1 1]);
        [nRows, nCols] = getSubplotLayout(nSubs);
        marg_h = [0.05 0.05];
        marg_w = [0.04 0.04];
        gap = [0.006, 0.006];
        yLims = [0, ceil(1.05*max(binCounts))];

        % text
        tx = axes('visible','off', 'position',[0 0 1 1]);
        mytitle = ['time histograms of  selected  anatomy areas, frequency = ' freqBands{freq,1} ': '];
        if strcmp(SA_params.signifValsNames{2},'hVals')
            title_signs = {'>','<'};
            for clz = 1:size(params.triggering.classes,1)
                clr = trialsData.info.colors(clz,:);
                mytitle = [mytitle '\color[rgb]{' num2str(clr(1)) ' ' num2str(clr(2)) ' ' num2str(clr(3)) '}' params.triggering.classes{1,1} ' ' title_signs{clz} ' ' params.triggering.classes{2,1} ', '];
            end    
            mytitle(end-1:end) = [];    % removes ', '
        elseif strcmp(SA_params.signifValsNames{2},'hVals_base')
            title_signs = {'>','<'};
            for clz = 1:size(title_signs,2)
                clr = trialsData.info.colors(clz,:);
                mytitle = [mytitle '\color[rgb]{' num2str(clr(1)) ' ' num2str(clr(2)) ' ' num2str(clr(3)) '}' params.triggering.classes{c,1} ' ' title_signs{clz} ' baseline, '];
            end    
            mytitle(end-1:end) = [];    % removes ', '
        end
            
        mytitle = strrep(mytitle, '_','\_');
        text(0.016, 0.98, mytitle, 'fontsize', 16, 'fontw', 'bold');

        % >>> plot histogram of times distribution <<<
        nPlot = 1;
        for aa = 1:size(t,1)
            if size(t{aa,1},2) >= nChnls_cutoff || size(t{aa,2},2) >= nChnls_cutoff    % at least 'nChnls_cutoff' (e.g. 2) different channels contribute
                h = subtightplot(nRows, nCols, nPlot, gap, marg_h, marg_w);
                hold on;

                % significant times when: class_1 > class_2 or vice versa
                for clz = 1:size(t,2)
                    if size(t{aa,clz},2) >= nChnls_cutoff
                        h_hist = histogram(t{aa,clz}, tBins_vector);
                        h_hist.FaceColor = trialsData.info.colors(clz,:);

                        % median of distribution
                        m = median(t{aa,clz},2);
                        plot([m, m], yLims, '--', 'LineWidth',2, 'Color', trialsData.info.colors(clz,:));
                    end
                end
                set(h, 'ylim', yLims);

                % plot ticks
                for xtck = [-1.0, -0.5, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]
                    if xtck > min(xlim) && xtck < max(xlim)
                        plot([xtck, xtck], yLims, '--', 'Color', [0.4,0.4,0.4]);
                    end
                end            
                plot([0, 0], yLims, ':k', 'LineWidth',1.5);

                % labels
                xlabel('time (s)');
    %             str_title = [AA_list{aa} ': nCh = ' num2str(size(t{aa},2)) ', nSubj = ' num2str(numel(unique(sbj{aa})))];
                str_title = cat(2,{AA_list{aa}}, {['nCh = ' num2str(chnls(aa,1)+chnls(aa,2)) ', nSubj = ' num2str(numel(unique(cat(2,sbj{aa,1},sbj{aa,2}))))]});
                %title(str_title);
                xLims = xlim;
                pos_text = [xLims(1)+0.5*diff(xLims), yLims(1)+0.9*diff(yLims)];
                text(pos_text(1),pos_text(2), str_title, ...
                    'FontSize',12, 'FontWeight','bold', 'HorizontalAlignment','center', 'Color','k', 'Interpreter','none');     

                nPlot = nPlot+1;
            end
        end

        %% save figure
        figname=['tHist_AA_clz' num2str(c) '_' freqBands{freq,1}];
        set(f, 'PaperPositionMode','auto');
        saveas(f, [outDir filesep figname '.fig']);
        print(f, '-dpng','-r0', [outDir filesep figname '.png']);
        %print(f, '-dpng','-r600', [outDir filesep figname '.png']);
        %print(f, '-dtiff', '-noui', '-r0', [outDir filesep figname '.tif']);
        disp(['Figure: ' figname ' stored in: ' outDir]);        
        close(f);

        %% ------plot sorted anatomic areas (AAs) on their activations times----
        f = figure('units','normalized','outerposition',[0 0 1 1], 'Color','w');
        title_signs = {'>','<'};

        % text
        tx = axes('visible','off', 'position',[0 0 1 1]);
        mytitle = ['times sorted on maximal bin count of selected  anatomy areas, frequency = ' freqBands{freq,1}];    
        mytitle = strrep(mytitle, '_','\_');
        text(0.016, 0.98, mytitle, 'fontsize', 16, 'fontw', 'bold');

        % plot both classes (bigger/smaller)
        for clz = 1:size(t,2)

            % cat all histogram values of used AAs (might differ for each class)
            histVals_used = [];
            anatArea_used = [];
            n = 1;
            for aa = 1:size(t,1)
                if ~any(isnan(histVals(:,aa,clz)))
                    histVals_used = cat(2, histVals_used, histVals(:,aa,clz));
                    anatArea_used{n} = [AA_list{aa,1} ' (' num2str(chnls(aa,clz)) '/' num2str(numel(unique(sbj{aa,clz}))) ')'];
                    n = n+1;
                end
            end
            histVals_used = histVals_used'; % used AAs (= rows) & time with binned counts (= cols)

            % normalize histogram values
            histVals_used = histVals_used./repmat(max(histVals_used,[],2), [1,size(histVals_used,2)]);

            % sort according to time with the maximal count
            [m,i_max] = max(histVals_used,[],2);
            [b,i_sorted] = sort(i_max);
            histVals_sorted = histVals_used(i_sorted,:);

            %  >>> plot sorted & normalized bin counts <<<
            %subplot(1,2,clz);
            subtightplot(1, 2, clz, gap+0.15, marg_h+0.08, marg_w+0.08);
            hold on;
            %imagesc(histVals_sorted);
            imagesc(getBinCenters(tBins_vector), [0.5:size(histVals_sorted,1)], histVals_sorted, [0 1]);
            h_clrbar = colorbar;
            ylabel(h_clrbar, 'normalized count of significant differences', 'FontSize',11, 'FontWeight','bold');

            axis tight
            if clz == 1
                colormap(gca,flip(hot(256),1));   % colormap(gca,brewermap(256,'*RdBu'));
            else
                colormap(gca, cat(1, [1 1 1], brewermap(256,'PuBu'), [0 0 0]));
            end

            % plot x-ticks
            for xtck = [-1.0, -0.5, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]
                if xtck > min(xlim) && xtck < max(xlim)
                    plot([xtck, xtck], ylim, '--', 'Color', [0.6,0.6,0.6]);
                end
            end            
            plot([0, 0], ylim, ':k', 'LineWidth',1.5);

            for ytck = 0:size(histVals_sorted,1)
                plot(xlim,[ytck, ytck], '-', 'Color', [0.6,0.6,0.6]);
            end

            set(gca, 'Ytick',[1:size(histVals_sorted,1)]-0.5);
            set(gca, 'YTickLabels', anatArea_used(i_sorted));
            box on;
             % title
            if strcmp(SA_params.signifValsNames{2},'hVals')
                mytitle = [params.triggering.classes{1,1} ' ' title_signs{clz} ' ' params.triggering.classes{2,1}];
            elseif strcmp(SA_params.signifValsNames{2},'hVals_base')
                mytitle = [params.triggering.classes{c,1} ' ' title_signs{clz} ' baseline'];
            end
            clr = trialsData.info.colors(clz,:);
            title(mytitle, 'Interpreter','none', 'Color',clr);
            xlabel('time (s)');
        end

        %% save figure
        figname=['sortedMaxBinCount_AA_clz' num2str(c) '_' freqBands{freq,1}];
        set(f, 'PaperPositionMode','auto');
        saveas(f, [outDir filesep figname '.fig']);
        print(f, '-dpng','-r0', [outDir filesep figname '.png']);
        %print(f, '-dpng','-r600', [outDir filesep figname '.png']);
        %print(f, '-dtiff', '-noui', '-r0', [outDir filesep figname '.tif']);
        disp(['Figure: ' figname ' stored in: ' outDir]);        
        close(f);
    end
end
