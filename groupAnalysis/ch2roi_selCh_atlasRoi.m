function [selCh_all, selCh_groups] = ch2roi_selCh_atlasRoi(params, nCh, useOnlySignif)
% TBD (now=Oct23)
% returns all selected channels 'selCh_all' & channel groups 'selCh_groups'
% returns channel indices 'ch' indexing the selected channels: selCh_H(ch)
% specified in: params.connectivity.selectedChnls
% based on atlasName & areaName in the atlas, see: anatomicalAreas_getList.m
% for example: 
% params.connectivity.selectedChnls = {...    % each row = tuple: atlasName, areaName (must match names in anatomicalAreas_getList.m) 
%     'Yeo7','Default'; ...
%     'Yeo7','Dorsal Attention'; ...  
%     };  

% (c) Jiri, Sep22

if nargin < 3
    useOnlySignif = false;
end

%% atlas & ROI (default)
if ~isfield(params.connectivity, 'selectedChnls')
    params.connectivity.selectedChnls = {...    % each row = tuple: atlasName, areaName (must match names in anatomicalAreas_getList.m) 
        'Yeo7','Default'; ...
        'Yeo7','Dorsal Attention'; ...  
        };  
end

%% init
selCh_all = [];
selCh_groups = cell(size(params.connectivity.selectedChnls,1),1);

%% go for each ROI (= 1 row in params.connectivity.selectedChnls)
for roi = 1:size(params.connectivity.selectedChnls,1)
    brainAtlas_name = params.connectivity.selectedChnls{roi,1};
    sel_ROI = params.connectivity.selectedChnls{roi,2};
    if strcmp(brainAtlas_name, 'all')       % special case: all channels
        selCh_groups{roi} = 1:nCh;
        selCh_all = 1:nCh;
    else    
        % load H & selCh_H
        load(params.storage.cacheFile, 'H', 'selCh_H_resp');
        selCh_H = selCh_H_resp;
        assert(nCh == size(selCh_H,2));
        
        % anatomical assignements from isarg_atlas
        [ass_isargAtlas, H] = updateHeaders_atlasInfo(params, selCh_H_resp, H, params.storage.cacheFile);

        % --- significance of activations (soa) from trialsData: all channels  ---
        if useOnlySignif
            soa_allCh = ch2roi_soa_trialsData(params);
        end
        
        % assign each channel
        for ch = 1:size(selCh_H,2)
            % --- significance of activations: selected channel -> ch_wasSignif
            if useOnlySignif
                ch_wasSignif = soa_allCh(ch);
            else
                ch_wasSignif = true;    % significance plays no role -> set to true
            end
            
            % --- anatomic area assignment -> thisAA  
            thisAA = ch2roi_assignAA(brainAtlas_name, ass_isargAtlas(ch), params.connectivity);
            
            % --- add to selected channel group? ---
            if strcmp(thisAA, sel_ROI) && ch_wasSignif  % is it the selected ROI?
                selCh_groups{roi} = cat(2, selCh_groups{roi}, ch);
                selCh_all = cat(2, selCh_all, ch);
            end
        end

    end
end
selCh_all = sort(unique(selCh_all));

% if isempty(selCh_all)
%     disp('No channels found for connectivity computation. Check selection in: params.connectivity.selectedChnls');
%     return;
% end
