function M_subj = plotFC2axes_5DM_spectral(M, EC_info, chGroups, thr_nChPerGroup)
% plots results of functional connectivity (FC) to existing axes - e.g. DTF
% works for two neural networks (NN)
% input:
%   - M = FC matrix, 5D: M(ch_trg,ch_src,freq,time,clz)
%   - EC_info = freq & time + other info (see: connectivityAnalysis.m)
%   - chGroups = struct with channel indices (see: ch2roi_getSelChnls.m)
% output:
%   - M_subj = FC subj mean, 5D: M_subj(grp_trg,grp_src,freq,time,clz)

% (c) Jiri, Dec23

if nargin < 4
    thr_nChPerGroup = 0;
end

if ~isfield(EC_info, 'plot_src_equalTo_trg')
    EC_info.plot_src_equalTo_trg = false;
end

nGroups = size(chGroups.ch_inds_trials,1);
% nClz = size(EC_info.info.clzNames,1);
nClz = size(M,5);

%% init output (mean over source & target channels)
M_subj = nan(nGroups,nGroups,size(M,3),size(M,4),size(M,5));    % 5D: grp_trg x grp_src x freq x time x clz

%% figure
% f = fig_make;
% nCols = nGroups;    % +1 => difference between NN 
% nRows = nClz;
% nPlot = 1;

nCols = EC_info.nCols;
nRows = EC_info.nRows;
marg_h = [0.08 0.08];   % margin from [bottom, top]
marg_w = [0.03 0.03];   % margin from [L, R] side
gap = [0.02, 0.02];     % between axes from [top, side]

%% plot DTF: source (src) -> target (trg)
nPlot = 1;
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
        if size(ch_src,2) >= thr_nChPerGroup && size(ch_trg,2) >= thr_nChPerGroup
            if (src ~= trg) || EC_info.plot_src_equalTo_trg
                % src->trg from M -> FC_src_trg
                FC_src_trg = squeeze(mean(mean(M(ch_trg,ch_src,:,:,:),1),2));   % 3D: freq x t x clz
                FC_trg_src = squeeze(mean(mean(M(ch_src,ch_trg,:,:,:),1),2));   % 3D: freq x t x clz
                
                if any(~isnan(FC_trg_src(:))) || any(~isnan(FC_src_trg(:)))
                    % cLims
                    cLims = getYLims(FC_src_trg,3);

                    % ----------- plot each clz -------------
                    for clz = 1:nClz
                        thisPlot = EC_info.thisCol + (nPlot-1)*EC_info.nCols;   % plot clz below each other
%                         subplot(EC_info.nRows,EC_info.nCols,thisPlot);
                        ax = subtightplot(nRows, nCols, thisPlot, gap, marg_h, marg_w); 
                        hold on;
        %                 if clz > 1
        %                     set(gca, 'Position', get(gca,'Position')+[-0.01 0 0 0]);
        %                 end

                        % plotting into -> EC_info
                        EC_info.thisPlot = thisPlot;
                        EC_info.str_title = plotFC_getTitle(chGroups, src, trg, EC_info.info.clzNames{clz}, EC_info.fcName);
                        
                        % xlabel: only bottom row 
                        if thisPlot > nCols*(nRows-1)
                            EC_info.plot_xlabel = true;
                        else
                            EC_info.plot_xlabel = false;
                        end
                        
                        % ylabel: only L column
                        if mod(thisPlot,nCols) == 1
                            EC_info.plot_ylabel = true;
                        else
                            EC_info.plot_ylabel = false;
                        end
        
                        % >>> plot FC <<<
                        plotEC_spectrum2axes(FC_src_trg(:,:,clz), EC_info.freq, EC_info.time, cLims, clz, EC_info);
                        nPlot = nPlot+1;
                        
                        % FC subj mean -> M_subj = 5D: grp_trg x grp_src x freq x time x clz x groups
                        M_subj(trg,src,:,:,clz) = FC_src_trg(:,:,clz);
                    end           
                end
%                 nPlot = nPlot+1;
            end
            
            % ----------- plot TRG->SRC - SRC->TRG -------------
%             thisPlot = nPlot + (clz-1)*nCols;   % plot clz below each other
            
        else
            M_subj = [];
%             close(f);
            return;
        end
    end
end
