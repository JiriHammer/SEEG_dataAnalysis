function chGroups = getChannelGroups(params, groupSpecification)
% returns channel groups (cell array) based on group specification

% (c) Jiri, Apr16

chGroups = [];

%% load 'H' & selCh
name_selCh = ['selCh_H_' params.thisData];
load(params.storage.cacheFile, 'H', name_selCh);
selCh_H = eval(name_selCh);
selCh_raw = 1:size(selCh_H,2);      % corresponds to channels in 'pred_raw' or 'resp_raw'

%% CAR: channel groups = headboxes of amplifier with different REFs & GNDs
if strcmp(groupSpecification, 'perHeadbox')
    if isfield(H.channels(1), 'headboxNumber')      % re-referencing per headbox
        hdbxAll = [];
        for ch = 1:size(H.channels,2)
            hdbxAll = cat(2, hdbxAll, H.channels(ch).headboxNumber);        % headbox numbers from all channels
        end
        hdbxSel = hdbxAll(selCh_H);                                         % headbox numbers from selected channels (e.g. iEEG channels)
        assert(size(selCh_raw,2) == size(hdbxSel,2));

        hdbxUnq = unique(hdbxSel);                                          % headbox numers
        chGroups = cell(1,size(hdbxUnq,2));
        for grp = 1:size(hdbxUnq,2)
            thisGrp = hdbxUnq(grp);
            for ch = 1:size(selCh_raw,2)
                if hdbxSel(ch) == thisGrp
                    chGroups{grp} = cat(2, chGroups{grp}, ch);
                end
            end
        end
    else                                            % only 1 headbox in recording
        chGroups{1,1} = selCh_raw;
    end
    
end

%% CAR: channel groups = SEEG electrode shanks
if strcmp(groupSpecification, 'perElectrode')
    elsh_all = [];
    for ch = 1:size(H.channels,2)
        elsh_all{ch} = extractFromString(H.channels(ch).name, 'string');    % string part of all channel names
    end    
    elsh_sel = elsh_all(selCh_H);                                           % string part of selected channels
    elsh_unq = unique(elsh_sel);                                            % electrode shank names
    
    chGroups = cell(1,size(elsh_unq,2));
    for grp = 1:size(elsh_unq,2)
        thisGrp = elsh_unq{grp};
        for ch = 1:size(selCh_raw,2)
            if strcmp(extractFromString(H.channels(ch).name, 'string'), thisGrp)
                chGroups{grp} = cat(2, chGroups{grp}, ch);
            end
        end
    end    
end

%% BIP: channel groups = neighboring SEEG channels on same electrode shank
if strcmp(groupSpecification, 'bip')
    %%
    selCh_H_bip = [];           % needed??? tbd?
    chGroups = [];
    grp = 1;                    % indexing of ch-pairs ("groups"), from 1:N
    channels_bip = struct;      % re-define H.channels
    ch = 2;
    while ch <= size(selCh_H,2) 
%     for ch = 2:size(selCh_H,2)
        chPair_1 = selCh_H(ch-1);
        chPair_2 = selCh_H(ch);
        elName_1 = extractFromString(H.channels(chPair_1).name, 'string');
        elName_2 = extractFromString(H.channels(chPair_2).name, 'string');
        elNum_1 = str2double(extractFromString(H.channels(chPair_1).name, 'number'));
        elNum_2 = str2double(extractFromString(H.channels(chPair_2).name, 'number'));        
        if strcmp(elName_1, elName_2)       % same electrode shank
            if abs(elNum_1 - elNum_2) == 1     % direct neighbors (TO DO: nearest neighbor distance)
                chGroups{grp} = [ch-1, ch];
%                 selCh_H_bip = cat(2, selCh_H_bip, [chPair_1; chPair_2]);    % old, tbd
                selCh_H_bip = cat(2, selCh_H_bip, grp);
                
                % update channels_bip struct
                channels_bip(grp).name = [H.channels(chPair_1).name '-' H.channels(chPair_2).name];
                channels_bip(grp).numberOnAmplifier = [H.channels(chPair_1).numberOnAmplifier; H.channels(chPair_2).numberOnAmplifier];
                channels_bip(grp).signalType = H.channels(chPair_1).signalType; % taken only from chPair 1
                channels_bip(grp).MNI_x = mean([H.channels(chPair_1).MNI_x, H.channels(chPair_2).MNI_x],2);
                channels_bip(grp).MNI_y = mean([H.channels(chPair_1).MNI_y, H.channels(chPair_2).MNI_y],2);
                channels_bip(grp).MNI_z = mean([H.channels(chPair_1).MNI_z, H.channels(chPair_2).MNI_z],2);
                channels_bip(grp).ass_brainAtlas = [H.channels(chPair_1).ass_brainAtlas '-' H.channels(chPair_2).ass_brainAtlas];
                channels_bip(grp).ass_cytoarchMap = [H.channels(chPair_1).ass_cytoarchMap '-' H.channels(chPair_2).ass_cytoarchMap];
                if isfield(H.channels(chPair_1), 'neurologyLabel')
                    channels_bip(grp).neurologyLabel = [H.channels(chPair_1).neurologyLabel '-' H.channels(chPair_2).neurologyLabel];
                else
                    channels_bip(grp).neurologyLabel = 'n.a.';
                end
                channels_bip(grp).esm = [H.channels(chPair_1).esm '-' H.channels(chPair_2).esm];
                channels_bip(grp).seizureOnset = any([H.channels(chPair_1).seizureOnset,  H.channels(chPair_2).seizureOnset]);
                channels_bip(grp).interictalOften = any([H.channels(chPair_1).interictalOften,  H.channels(chPair_2).interictalOften]);
                if isfield(H.channels(chPair_1), 'lesion')
                    channels_bip(grp).lesion = any([H.channels(chPair_1).lesion,  H.channels(chPair_2).lesion]);
                else
                    channels_bip(grp).lesion = false;
                end
                               
                disp([' - BIP ' num2str(grp) ' = ' H.channels(chPair_1).name '-' H.channels(chPair_2).name]);
                grp = grp+1;
            end
        end
        
        % which neighboring contacts to inlcude?
        if strcmp(params.bip.chGroups, '12_23')
            ch = ch+1;      % inlcudes all bip. pairs (e.g. A2-A1 & A3-A2, so that A2 is included twice)
        elseif strcmp(params.bip.chGroups, '12_34')
            ch = ch+2;      % includes one bip. pair! (e.g. A2-A1 & A4-A3, so that A2 is included once!)
        end
    end     % of while
    
    %% save re-defined pointers of channels in H struct
    if params.thisSess == size(params.storage.sessionCacheFiles,2)
        eval([name_selCh '=selCh_H_bip;' ]);        % = 'selCh_H_pred' or 'selCh_H_resp'
        H.channels = channels_bip;
        save(params.storage.cacheFile, name_selCh, 'H', '-append');    
        disp(' - BIP referencing: re-definition of H.channels !!!');
    end
end

