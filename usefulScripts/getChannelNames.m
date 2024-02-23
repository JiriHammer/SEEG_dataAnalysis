function chNames = getChannelNames(params, signalType, dataTag, isShort)
% determines channel names based on signal type or naming style (specified in 'signal_style') 
% for example: signal_style = eeg, ecog, seeg
% or anatomical assignment in structure H

% (c) Jiri, Apr16, Jun17

chNames = [];
if nargin < 3
    dataTag = 'resp';
end

if nargin < 4
    isShort = true;
end

%% trajectory?
if ismember(lower(signalType),{'xpos','xvel','xacc','abspos','absvel','absacc','dddt_pos','dddt_neg'})
    chNames = {signalType};
    return;
end

%% copy+noise model?
if ismember(lower(signalType),{'cnmodel'})
    assert(params.nCh == 2);    % ch1 = cnModel & ch2 = noiseOnly
    chNames = cell(1,params.nCh);
    chNames{1} = ['cnModel = ' params.model_copyNoise.copy{1}];
    for ch = 2:params.nCh-1
        chNames{1} = [chNames{1} ' + ' params.model_copyNoise.copy{ch}];
    end    
    chNames{2} = 'cnModel = noise';
    return;
end

%% sua simulation?
if ismember(lower(signalType),{'suasim'})
    if ismember('suaSimulation_populations', params.response.dataProcessing)
        assert(size(params.suaSimulation.popSizes,2) == params.nCh);
        chNames = cell(1, params.nCh);
        for ch = 1:params.nCh
            chNames{ch} = ['N = ' num2str(params.suaSimulation.popSizes(ch))];
        end
    else
        for ch = 1:params.nCh
            chNames{ch} = ['ch ' num2str(ch)];
        end
    end
    return;
end

%% load: 'H', 'selCh_H', 'selCh_names'
if ~isfield(params, 'H')
    % load struct H
    load(params.storage.cacheFile, 'H');
    
    % load selCh_H
    varName = ['selCh_H_' dataTag];
    load(params.storage.cacheFile, varName);
    selCh_H = eval(varName);
    
    % load selCh_names
%     varName = ['selCh_names_' dataTag];
%     load(params.storage.cacheFile, varName);
%     selCh_names = eval(varName);    
else
    selCh_H = params.selCh_H;
    H = params.H;
end
assert(params.nCh == size(selCh_H,2));

%% determine iEEG channel names: anatomical atlas
if strcmpi(signalType,'ieeg')
    chNames = cell(1, size(selCh_H,2));
    for ch = 1:size(selCh_H,2)
        for n = 1:size(selCh_H,1)
            thisCh = selCh_H(n,ch);
            %chNames{ch} = [chNames{ch} H.channels(thisCh).name '(' num2str(thisCh) ')-'];
            chNames{ch} = [chNames{ch} H.channels(thisCh).name '-'];
        end
        chNames{ch}(end) = ':';     % = 'P3(3):' or 'P3(3)-P4(4):' (for example)
        
        for n = 1:size(selCh_H,1)
            thisCh = selCh_H(n,ch);
            maxLngth = 5;
            tmp = [];

            % laterality
            if H.channels(thisCh).MNI_x >= 0
                laterality = 'R';
            elseif H.channels(thisCh).MNI_x < 0
                laterality = 'L';
            end       
            
            % preferred ESM label?
            if isfield(params, 'prefer_ESM')
                if isfield(H.channels(thisCh), 'esm')
                    tmp = H.channels(thisCh).esm;
                end
            end 

            % (A) neurologist's label
            if strcmp(tmp, 'n.a.') || strcmp(tmp, 'n.s.r.') || isempty(tmp)
                if isfield(H.channels(thisCh), 'neurologyLabel')
                    tmp = H.channels(thisCh).neurologyLabel;
                end
            end

            % (B) cytoarchitectonic map
            if strcmp(tmp, 'n.a.') || strcmp(tmp, 'n.s.r.') || isempty(tmp)
                if isfield(H.channels(thisCh), 'ass_cytoarchMap')
                    tmp = H.channels(thisCh).ass_cytoarchMap;
                end
            end

            % (C) brain atlas = macro anatomy label
            if strcmp(tmp, 'n.a.') || strcmp(tmp, 'n.s.r.') || isempty(tmp)
                if isfield(H.channels(thisCh), 'ass_brainAtlas')
                    tmp = H.channels(thisCh).ass_brainAtlas;
                    maxLngth = 4;
                end
            end

            % shorten the name ?
            if isShort
                tmp = shortenChannelName(tmp, maxLngth);
            end
            chNames{ch} = [chNames{ch} ' ' laterality ' ' tmp '-'];
            %chNames{ch} = [num2str(thisCh) ': ' tmp];
        end
        chNames{ch}(end) = '';
    end
end

%% determine channel names: EEG
if strcmpi(signalType,'eeg')
    chNames = cell(1, size(selCh_H,2));
    for ch = 1:size(selCh_H,2)
        thisCh = selCh_H(ch);
        assert(strcmp(H.channels(thisCh).signalType,'EEG'));
        chNames{ch} = [num2str(thisCh) ': ' H.channels(thisCh).name];
    end
end

%% remove special format chars
%chNames = strrep(chNames, '_','\_');

