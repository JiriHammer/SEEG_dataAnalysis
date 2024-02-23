function plotEC_chGroups_imagesc(params, M, EC_info, selCh_groups, subjTag)
% plots effective connectivity (EC): a-la spectra
% works only for 2 groups !!!
% tbd ?

% (c) Jiri, Nov22

if ~isfield(EC_info, 'methodName')
    EC_info.methodName = params.connectivity.connectivityMeasure{1};
end

%% channels group selection
[selCh_all, selCh_groups] = ch2roi_selCh_atlasRoi(params, size(M,1));

%% init struct spectra (for connectivity plotting)
freqAxis = EC_info.freq';     % assumes row vector, 1 x N
clzLabels = EC_info.clzLabels;

disp([EC_info.methodName ' connectivity analysis plotting: subj = ' subjTag ' ...']);
assert(size(M,5) == size(clzLabels,2));

%% figure (only for 2 groups)
f = fig_make;
% nCols = numel(find(L == 1));  % each group to each group
% nRows = size(clzLabels,2);  % classes below each other

nRows = 2;  % each group to each other(!) group
nCols = size(clzLabels,2);  % classes
nPlot = 1;

% colormap
% colormap(f,brewermap(256,'*RdBu'));  
% colormap(f,hot(256));
colormap(f,jet(256));

% color limits
cVals = [];
for clz = 1:size(clzLabels,2)
    tmp = M(selCh_groups{1},selCh_groups{2},:,:,clz);
    cVals = cat(1, cVals, tmp(:));
    tmp = M(selCh_groups{2},selCh_groups{1},:,:,clz);
    cVals = cat(1, cVals, tmp(:));    
end
cLims = getYLims(cVals, 1);

%  ----------------------- plot EC: group1 -> group2 ----------------------
for clz = 1:size(clzLabels,2)
    g1 = 1;
    g1_ch = selCh_groups{g1};         % indices of channels in COH matrix
    groupName_1 = params.connectivity.selectedChnls{g1,2};
    if ~isempty(g1_ch)
        g2 = 2;
        g2_ch = selCh_groups{g2}; % indices of channels in M matrix
        groupName_2 = params.connectivity.selectedChnls{g2,2};
        if ~isempty(g2_ch)

            % data (mean over group2 & group1): G1 -> G2
            C = squeeze(nanmean(nanmean(M(g1_ch,g2_ch,:,:,clz),2),1));       % 2D: F x T
%             cLims = getYLims(C(:), 3);
            
            % axes
            subplot(nRows, nCols, nPlot);
            hold on;
            box on;

            EC_info.src_groupName = params.connectivity.selectedChnls{g1,2};
            EC_info.trg_groupName = params.connectivity.selectedChnls{g2,2};
            plotEC_spectrum2axes(C, freqAxis, EC_info.time, cLims, clz, EC_info);

            nPlot = nPlot+1;
        end
    end
end

%  ----------------------- plot EC: group2 -> group1 ----------------------
for clz = 1:size(clzLabels,2)
    g1 = 1;
    g1_ch = selCh_groups{g1};         % indices of channels in COH matrix
    groupName_1 = params.connectivity.selectedChnls{g1,2};
    if ~isempty(g1_ch)
        g2 = 2;
        g2_ch = selCh_groups{g2}; % indices of channels in COH matrix
        groupName_2 = params.connectivity.selectedChnls{g2,2};
        if ~isempty(g2_ch)

            % data (mean over group2 & group1): G2 -> G1
            C = squeeze(nanmean(nanmean(M(g2_ch,g1_ch,:,:,clz),2),1));       % 2D: F x T
%             cLims = getYLims(C(:), 3);
            
            % axes
            subplot(nRows, nCols, nPlot);
            hold on;
            box on;

            EC_info.src_groupName = params.connectivity.selectedChnls{g2,2};
            EC_info.trg_groupName = params.connectivity.selectedChnls{g1,2};
            plotEC_spectrum2axes(C, freqAxis, EC_info.time, cLims, clz, EC_info);
            
            nPlot = nPlot+1;
        end
    end
end

%% save figure
figname = [subjTag '_' EC_info.methodName '_' groupName_1 '_' groupName_2];
outDir = [params.storage.dir_results filesep EC_info.methodName];
fig_save(f, figname, outDir);
close(f); 



