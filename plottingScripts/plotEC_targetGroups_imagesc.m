function EC_avg = plotEC_targetGroups_imagesc(params, M, EC_info, selCh_groups, subjTag)
% plots effective connectivity (EC) of selected channel group (= params.grup_source) to all other groups
% plots = imagesc (a-la spectra) for each class: src->trg & trg->src
% input: M = 5D: ch_src, ch_trg,freq, time, clz
% channel groups indices (indexing M) in cell: selCh_groups{group}
% output: EC_avg{group,clz} = averaged EC (src->trg) for each target group (or empty if no channels found)

% (c) Jiri, Nov22

if ~isfield(EC_info, 'methodName')
    EC_info.methodName = params.connectivity.connectivityMeasure{1};
end

%% number of groups to plot -> N_groups
N_groups = 0;
for g2 = 1:size(selCh_groups,1)
    if g2 ~= params.grup_source
        if ~isempty(selCh_groups{g2})
            N_groups = N_groups+1;
        end
    end
end
if N_groups == 0
    disp('no channel groups found.');
    return;
end

%% apply z-score?
if params.connectivity.plotEC_zscore
    M = zscore(M, 0, 4);    % dim 4 = time
end

%% init struct spectra (for connectivity plotting)
freqAxis = EC_info.freq';     % assumes row vector, 1 x N
clzLabels = EC_info.clzLabels;

disp([EC_info.methodName ' connectivity analysis plotting: subj = ' subjTag ' ...']);
assert(size(M,5) == size(clzLabels,2));
g1 = params.grup_source;

%% figure (only for 2 groups)
f = fig_make;

% nRows = N_groups;  % each group to each other(!) group
% nCols = 2*size(clzLabels,2);  % classes
nCols = N_groups;  % each group to each other(!) group
nRows = 2*size(clzLabels,2);  % classes
nPlot = 1;

EC_info.nCols = nCols;
EC_info.nRows = nRows;

marg_h = [0.1 0.1];
marg_w = [0.04 0.04];
gap = [0.06, 0.02];

% colormap
% colormap(f,brewermap(256,'*RdBu'));  
% colormap(f,hot(256));
colormap(f,jet(256));

% color limits
cVals = [];
for clz = 1:size(clzLabels,2)
    for g2 = 1:size(selCh_groups,1)
        if g2 ~= params.grup_source
            if ~isempty(selCh_groups{g2})
                tmp = M(selCh_groups{g1},selCh_groups{g2},:,:,clz);
                cVals = cat(1, cVals, tmp(:));
                tmp = M(selCh_groups{g2},selCh_groups{g1},:,:,clz);
                cVals = cat(1, cVals, tmp(:));    
            end
        end
    end
end
cLims = getYLims(cVals, 1);

%% EC plot
EC_avg = cell(size(selCh_groups,1),size(clzLabels,2));  % 2D cell: groups x clz
g1_ch = selCh_groups{g1};       % source channels
for clz = 1:size(clzLabels,2)
    for g2 = 1:size(selCh_groups,1)
        if g2 ~= params.grup_source
            if ~isempty(selCh_groups{g2})
                g2_ch = selCh_groups{g2};   % target channels
            
                % ---------------------G1 -> G2----------------------------
                % data (mean over group2 & group1): G1 -> G2
                C = squeeze(nanmean(nanmean(M(g1_ch,g2_ch,:,:,clz),2),1));       % 2D: F x T
    %             cLims = getYLims(C(:), 3);
    
                % averaged EC (src->trg)
                EC_avg{g2,clz} = C;     % for each target group (but not trg->src !)

                % axes
                if nCols > 20
                    subplot(nRows, nCols, nPlot);
                else
                    subtightplot(nRows, nCols, nPlot, gap, marg_h, marg_w);
                end       
                EC_info.thisPlot = nPlot;
                hold on;
                box on;

                % >>> plot G1 -> G2 <<<
                EC_info.src_groupName = params.connectivity.selectedChnls{g1,2};
                EC_info.trg_groupName = params.connectivity.selectedChnls{g2,2};
                EC_info.src_Nch = size(selCh_groups{g1},2);
                EC_info.trg_Nch = size(selCh_groups{g2},2);
                plotEC_spectrum2axes(C, freqAxis, EC_info.time, cLims, clz, EC_info);

                % ---------------------G2 -> G1----------------------------
                % data (mean over group2 & group1): G2 -> G1 ("reversing" order of source & target)
                C = squeeze(nanmean(nanmean(M(g2_ch,g1_ch,:,:,clz),2),1));       % 2D: F x T
    %             cLims = getYLims(C(:), 3);

                % axes
                if nCols > 20
                    subplot(nRows, nCols, nPlot+nCols);
                else
                    subtightplot(nRows, nCols, nPlot+nCols, gap, marg_h, marg_w);
                end      
                EC_info.thisPlot = nPlot+nCols;
                hold on;
                box on;

                % >>> plot G2 -> G1 <<<
                EC_info.src_groupName = params.connectivity.selectedChnls{g2,2};
                EC_info.trg_groupName = params.connectivity.selectedChnls{g1,2};
                EC_info.src_Nch = size(selCh_groups{g2},2);
                EC_info.trg_Nch = size(selCh_groups{g1},2);                
                plotEC_spectrum2axes(C, freqAxis, EC_info.time, cLims, clz, EC_info);
                
                nPlot = nPlot+1;
            end
        end
    end
    nPlot = 2*nCols+1;
end
          
%% colorbar
hCb = colorbar;
if params.connectivity.plotEC_zscore
    clrbar_title = [EC_info.methodName ' (zscore)'];
else
    clrbar_title = EC_info.methodName;
end
ylabel(hCb, clrbar_title, 'FontSize',10);
set(hCb, 'position', [0.96 0.2 0.01 0.6]);

%% save figure
figname = [subjTag '_' EC_info.methodName '_' num2str(g1) '_' params.connectivity.selectedChnls{g1,2}];
if params.connectivity.plotEC_zscore
    outDir = [params.storage.dir_results filesep EC_info.methodName filesep params.atlasName '_g2all_zscored'];
else
    outDir = [params.storage.dir_results filesep EC_info.methodName filesep params.atlasName '_g2all'];
end
fig_save(f, figname, outDir);
close(f); 


