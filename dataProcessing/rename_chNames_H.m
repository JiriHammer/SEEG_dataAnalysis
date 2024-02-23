function rename_chNames_H(params, filterName, dataTag)
% renames selCh_names (selCh_names_pred or selCh_names_resp)
% needed after spatial filter is applied

% (c) Jiri, Nov18

%% load 'H' & selCh
varName = ['selCh_H_' dataTag];              % = 'selCh_H_pred' or 'selCh_H_resp'
load(params.storage.cacheFile, 'H', varName);
selCh_H = eval(varName);

%% load selCh_names (default, in case spatial filter = nan)
varName = ['selCh_names_' dataTag];          % = 'selCh_names_pred' or 'selCh_names_resp'
load(params.storage.cacheFile, 'H', varName);
selCh_names = eval(varName);

%% CAR: rename selected channels
if strfind(filterName, 'car')        % ~ common average re-reference (CAR)
    selCh_names = cell(size(selCh_H));
    for ch = 1:size(selCh_H,2)
        assert(strcmp(selCh_names{1,ch}, H.channels(selCh_H(ch)).name));
        selCh_names{2,ch} = 'CAR';
    end
end

%% BIP: rename selected channels
if strcmp(filterName, 'bip')        % ~ bipolar reference (BIP)
    selCh_names = cell(size(selCh_H));
    for ch = 1:size(selCh_H,2)
        selCh_names{1,ch} = H.channels(selCh_H(1,ch)).name;
        selCh_names{2,ch} = H.channels(selCh_H(2,ch)).name;
        disp([' - BIP: ' H.channels(selCh_H(1,ch)).name ' - ' H.channels(selCh_H(2,ch)).name]);
    end
end

%% save
varName = ['selCh_names_' dataTag];          % = 'selCh_names_pred' or 'selCh_names_resp'
eval([varName '=selCh_names;' ]);         
save(params.storage.cacheFile, varName, '-append');


