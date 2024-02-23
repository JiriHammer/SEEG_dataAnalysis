function M_subj = plotFC_spectral(M, EC_info, chGroups)
% plots results of functional connectivity (FC) - e.g. DTF
% input:
%   - M = FC matrix, 5D: M(ch_trg,ch_src,freq,time,clz)
%   - EC_info = freq & time + other info (see: connectivityAnalysis.m)
%   - chGroups = struct with channel indices (see: ch2roi_getSelChnls.m)
% output:
%   - M_subj = FC subj mean, 5D: M_subj(grp_trg,grp_src,freq,time,clz)

% (c) Jiri, Oct23

nGroups = size(chGroups.ch_inds_trials,1);
nClz = size(EC_info.info.clzNames,1);

%% init output (mean over source & target channels)
M_subj = nan(nGroups,nGroups,size(M,3),size(M,4),size(M,5));    % 5D: grp_trg x grp_src x freq x time x clz x groups

%% figure
f = fig_make;
nCols = nGroups*nClz;
nRows = nGroups;
nPlot = 1;

%% plot DTF: source (src) -> target (trg)
for src = 1:nGroups
    if isfield(chGroups, 'ch_inds_M')
        ch_src = chGroups.ch_inds_M{src};       % DTF on sel chnls, source channel indices to matrix M
    else
        ch_src = chGroups.ch_inds_trials{src};  % DTF on all chnls, source channel indices to matrix M
    end    
    for trg = 1:nGroups
        if isfield(chGroups, 'ch_inds_M')
            ch_trg = chGroups.ch_inds_M{trg};       % DTF on sel chnls, target channel indices to matrix M
        else
            ch_trg = chGroups.ch_inds_trials{trg};  % DTF on all chnls, target channel indices to matrix M
        end
        if ~isempty(ch_src) && ~isempty(ch_trg)
            % src->trg from M -> FC_src_trg
            FC_src_trg = squeeze(mean(mean(M(ch_trg,ch_src,:,:,:),1),2));   % 3D: freq x t x clz
            
            for clz = 1:nClz
                subplot(nRows,nCols,nPlot);
                hold on;
                if clz > 1
                    set(gca, 'Position', get(gca,'Position')+[-0.01 0 0 0]);
                end
                
                % plotting into -> EC_info
                EC_info.thisPlot = nPlot;
                EC_info.nRows = nRows;
                EC_info.nCols = nCols;
                EC_info.str_title = plotFC_getTitle(chGroups, src, trg, EC_info.info.clzNames{clz});
                
                % plot FC
%                 cLims = [-0.2,0.2]; 
                cLims = getYLims(FC_src_trg,3);
                plotEC_spectrum2axes(FC_src_trg(:,:,clz), EC_info.freq, EC_info.time, cLims, clz, EC_info);
                
                % FC subj mean -> M_subj = 5D: grp_trg x grp_src x freq x time x clz x groups
                M_subj(trg,src,:,:,clz) = FC_src_trg(:,:,clz);
                
                nPlot = nPlot+1;
            end
            
        else
            for clz = 1:nClz
                nPlot = nPlot+1;
            end
        end
    end
end
                
% save
fig_save(f, EC_info.figName, EC_info.outDir, 'format','png', 'res',600);
close(f);   

            
            
                
        

return;
%% OLD, MSU? (now=Oct23)
i_ch_g1 = EC_info.selCh_groups{1};
i_ch_g2 = EC_info.selCh_groups{2};
FC_g1_g2 = squeeze(mean(mean(M(i_ch_g1,i_ch_g2,:,:,:),1),2));   % 3D: freq x t x clz
FC_g2_g1 = squeeze(mean(mean(M(i_ch_g2,i_ch_g1,:,:,:),1),2));   % 3D: freq x t x clz
M_subj = cat(4, FC_g1_g2, FC_g2_g1);                            % 4D: freq x t x clz x grp
if do_zscore
    M_subj = zscore(M_subj, 0, 2);
    cLims = [0 0.2];
end

%% interaction labels (hardcoded for i_ch_g1)
g_src = 1;
g_trg = 2;
FC_roiInteraction{1,1} = [FC_roi{g_src,2} '(' num2str(size(i_ch_g1,2)) ') -> ' FC_roi{g_trg,2} '(' num2str(size(i_ch_g2,2)) ')'];
g_src = 2;
g_trg = 1;        
FC_roiInteraction{2,1} = [FC_roi{g_src,2} '(' num2str(size(i_ch_g2,2)) ') -> ' FC_roi{g_trg,2} '(' num2str(size(i_ch_g1,2)) ')'];

%% plot as spectra for 2 channels
f = fig_make;
nCols = 2;
nRows = 2;
nPlot = 1;
for clz = 1:nClz
    for grp = 1:nGroups
        % axes
        subplot(nRows, nCols, nPlot);
        hold on;

        % plot: imagesc
        h = imagesc(EC_info.time, EC_info.freq, M_subj(:,:,clz,grp));
%                 h = imagesc(EC_info.time, EC_info.freq, M_subj(:,:,clz,grp), cLims);
        colormap(gca,brewermap(256,'*RdBu')); 
        axis tight
        colorbar;

        % time t = 0
        plot([0, 0], ylim, '--k');

        % labels
        xlabel('time (s)');
        ylabel('freq (Hz)');
        title(['DTF: subj = ' subjTag ', ROI = ' FC_roiInteraction{grp} ', clz = ' EC_info.info.clzNames{clz}]);
        nPlot = nPlot+1;
    end
end               

% save
figName = ['DTF_raw_subj_' subjTag];
outDir = [FC_outDir filesep 'DTF_fig'];
fig_save(f, figName, outDir, 'format','png', 'res',600);
close(f);   