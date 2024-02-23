%% plot values (SNR, DA) of different anatomical areas
% from different:
% - neuronal features (data_info.neurVars)
% - kinematic variables (data_info.trajVars)
% - subjects (subjList)

% (c) Jiri, Jun17

%% settings: SNR
pathBeg = 'G:\dox\proj_carGameAnalysis\data\kinVars_tuning';
figDir  = 'G:\dox\proj_carGameAnalysis\figs\kinVars_tuning\anatomicAreas';
if ~exist(figDir, 'dir')
    mkdir(figDir);
end  
trajVars = {'xpos','xvel','absPos','absVel'};
neurVars = {'lfc','hiGamma','beta'};
var2show = 'snr_vals';      % choices: 'snr_vals','snr_rank'
%var2show = 'snr_woCorrTerm';      % old implementation

signalTypes = {'ecog', 1;       % = signalType, signalNumber
               'seeg', 2};
AA_type = 'largeAA';        % 'individualAA', 'largeAA'

%% load SNR or DA values
dataStruct = struct;
dataStruct.var2show = var2show;             % determines if tuning or decoding results are loaded
dataStruct.subtractSNR_x_abs = true;        % for SNR: subtracts xPos - absPos & xVel - absVel
dataStruct.pathBeg = pathBeg;
dataStruct.figDir = figDir;
dataStruct.trajVars = trajVars;
dataStruct.neurVars = neurVars;
[d_nts, dataStruct] = loadPooledData(dataStruct);
subjList = dataStruct.subjList;
t_lag = dataStruct.t_lag;
H_all = dataStruct.H_all;

%% individual cytoarchitectonic (anatomical) areas
%cytoarchAreas = getLargerAnatomicalArea([], 'list_individualAA');

%% larger anatomical areas (e.g. cytoarchitectonic)
anatAreas = getLargerAnatomicalArea([], ['list_' AA_type]);

%% sort channels to anatomical areas
n_ass = 0;
n_tot = 0;
table_ch2aa = cell(size(subjList,1),1);     % each cell has rows: subj;chnl;aa_index
for subj = 1:size(subjList,1)
    for ch = 1:size(H_all{subj}.selCh,2)
        thisCh = H_all{subj}.selCh(ch);
        i_aa = assign_ch2aa(H_all{subj}.channels(thisCh), anatAreas, subjList{subj}, thisCh);
        table_ch2aa{subj} = cat(2, table_ch2aa{subj}, [subj; ch; i_aa]);  % each cell has rows: subj;chnl;i_aa
        if isnan(i_aa), n_ass = n_ass+1; end
        n_tot = n_tot+1;
    end
    display('+-------------------------------------------------------------+');
end
display(['Total number of channels = ' num2str(n_tot)]);
display([' - number of channels not assigned = ' num2str(n_ass) ' = ' num2str(n_ass/n_tot*100) ' %']);

%% assign signal type: ecog or seeg ?
table_ecog_seeg = cell(size(subjList,1),1);
for subj = 1:size(subjList,1)
    for ch = 1:size(H_all{subj}.selCh,2)
        thisCh = H_all{subj}.selCh(ch);
        was_assigned = false;
        for sig = 1:size(signalTypes,1)
            if ~isempty(strfind(lower(H_all{subj}.channels(thisCh).signalType),signalTypes{sig,1}))
                table_ecog_seeg{subj} = cat(2, table_ecog_seeg{subj}, [ch; signalTypes{sig,2}]);
                was_assigned = true;
            end
        end
        if ~was_assigned
            error([subjList{subj} ', ch = ' num2str(thisCh) ', unknown signal type: ' H_all{subj}.channels(thisCh).signalType]);
        end
    end
end

%% erase not used AAs (at least 2 diff. subjects contribute)
used_aa = [];
nSubj_aa_signal = zeros(size(anatAreas,1),size(signalTypes,1));
for aa = 1:size(anatAreas,1)
    for sig = 1:size(signalTypes,1)
        for subj = 1:size(subjList,1)   
            is_thisSIG = table_ecog_seeg{subj}(end,:) == signalTypes{sig,2};
            is_thisAA = table_ch2aa{subj}(end,:) == aa;
            if sum(is_thisAA & is_thisSIG) > 0
                nSubj_aa_signal(aa,sig) = nSubj_aa_signal(aa,sig) + 1;
            end
        end
    end
end
            

