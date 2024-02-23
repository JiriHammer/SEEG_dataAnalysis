function [d_nts, dataStruct] = loadPooledData(dataStruct)
% loads data (eg. SNR or DA), typically as 2D matrix = [lag x chnl]
% from different:
% - neuronal features (data_info.neurVars)
% - kinematic variables (data_info.trajVars)
% - subjects (subjList)
% returns 'd_nts' = 3D cell: neur x traj x subj
% where: d_nts{neur,traj,subj} = 2D: lag x chnl

% (c) Jiri, Jun17

%% defaults: neuronal features
if ~isfield(dataStruct, 'neurVars') 
    neurVars = {'lfc','hiGamma','beta'};
else
    neurVars = dataStruct.neurVars; 
end

%% defaults: kinematic variables
if ~isfield(dataStruct, 'neurVars') 
    trajVars = {'xpos','xvel','absPos','absVel'};
else
    trajVars = dataStruct.trajVars; 
end

%% defaults: subject list
if ~isfield(dataStruct, 'subjList') 
    subjList = {'11_FR1_day1'};
else
    subjList = dataStruct.subjList;
end

%% default variable to load. choices: 'snr_vals','snr_rank', 'da'
if ~isfield(dataStruct, 'var2show') 
    var2show = 'snr_vals'; 
else
    var2show = dataStruct.var2show;
end
if ~isfield(dataStruct, 'subtractSNR_x_abs') 
    dataStruct.subtractSNR_x_abs = false;
end

%% defaults settings
if ~isfield(dataStruct, 'pathBeg')
    pathBeg = 'G:\dox\proj_carGameAnalysis\data\kinVars_tuning';
else
    pathBeg = dataStruct.pathBeg;
end
if ~isfield(dataStruct, 'outName_prefix'), dataStruct.outName_prefix = 'tuning'; end
if ~isfield(dataStruct, 'plotFigs') 
    data_info.plotFigs = false; 
end
if ~isfield(dataStruct, 'figDir')
    figDir  = ['G:\dox\proj_carGameAnalysis\figs\kinVars_decoding\' data_info.neurVar '\AVG'];
else
    figDir = dataStruct.figDir;
end

if ~isfield(dataStruct, 'loadHeaderInfo') 
    dataStruct.loadHeaderInfo = true; 
end

%% alloction of result
disp('loading pooled data ...');
d_nts = cell(size(neurVars,2), size(trajVars,2),size(subjList,1));    % SNR or DA: 3D = [neur x traj x subj]
t_lag = [];
H_all = cell(size(subjList,1),1);

%% AVG: load over all neurs, trajs & subjects
if strfind(var2show, 'avg')
    bins_all = cell(size(neurVars,2), size(trajVars,2), size(subjList,1),1);
    for neur = 1:size(neurVars,2)   
        display(['loading tuning results for: ' neurVars{neur} ' ...']);
        for traj = 1:size(trajVars,2)
            for subj = 1:size(subjList,1)

                % AVG: single-subject & single-traj
                data_info = struct;
                data_info.trajVar = trajVars{traj};
                data_info.neurVar = neurVars{neur};
                data_info.pathBeg = pathBeg;
                data_info.outName_prefix = dataStruct.outName_prefix;
                data_info.figDir = figDir;
                data_info.subjTag = subjList{subj,1};
                data_info.var2show = var2show;
                %data_info.selLags = [-1, 1];
                data_info.plotFigs = false;            
                [trialsData,params] = tuning_fig_AVG(data_info);
                d_nts{neur,traj,subj} = trialsData.cVals;
                bins_all{neur,traj,subj} = trialsData.yVals;
                dataStruct.tunInfo = trialsData.info.tunInfo;   % binning info
                
                % lag values
                if isempty(t_lag)
                    t_lag = trialsData.xVals;
                else
                    %assert(size(t_lag,1) == size(trialsData.xVals,1));  % should be same lags for all
                    if size(t_lag,1) ~= size(trialsData.xVals,1)
                        display(['WARNING: ' subjList{subj,1} ', srate = ' num2str(params.srate) ', N lags = ' num2str(size(trialsData.xVals,1)) ', N expected = ' num2str(size(t_lag,1))]);
                    end
                end
                
                % update header info
                if dataStruct.loadHeaderInfo
                    if isempty(H_all{subj})
                        load(params.storage.cacheFile, 'H', 'selCh_H_resp');
                        if size(trialsData.info.tunInfo,2) == 1
                            assert(size(trialsData.cVals,3) == size(selCh_H_resp,2)); 
                        elseif size(trialsData.info.tunInfo,2) == 2
                            assert(size(trialsData.cVals,4) == size(selCh_H_resp,2)); 
                        end
                        H.selCh = selCh_H_resp;
                        params.path2others.isarg_atlas = 'G:\shareData\visualization_normBrains\isarg_atlas_v8';
                        [ass_isargAtlas, H] = updateHeaders_atlasInfo(params, selCh_H_resp, H, params.storage.cacheFile);
                        H_all{subj} = H;
                    end
                end
            end
        end
    end
    dataStruct.bins_all = bins_all;     % bin centers
end       
                
%% SNR: load over all neurs, trajs & subjects
if strfind(var2show, 'snr')
    for neur = 1:size(neurVars,2)   
        display(['loading tuning results for: ' neurVars{neur} ' ...']);
        for traj = 1:size(trajVars,2)
            for subj = 1:size(subjList,1)

                % SNR: single-subject & single-trajectory
                data_info = struct;
                data_info.trajVar = trajVars{traj};
                data_info.neurVar = neurVars{neur};
                data_info.pathBeg = pathBeg;
                data_info.outName_prefix = dataStruct.outName_prefix;
                data_info.figDir = figDir;
                data_info.subjTag = subjList{subj,1};
                data_info.var2show = var2show;
                %data_info.selLags = [-1, 1];
                data_info.plotFigs = false;
                [trialsData,params] = tuning_fig_SNR(data_info);  % !!! loading !!!
                d_nts{neur,traj,subj} = trialsData.yVals;
                
                % lag values
                if isempty(t_lag)
                    t_lag = trialsData.xVals;
                else
                    %assert(size(t_lag,1) == size(trialsData.xVals,1));  % should be same lags for all
                    if size(t_lag,1) ~= size(trialsData.xVals,1)
                        display(['WARNING: ' subjList{subj,1} ', srate = ' num2str(params.srate) ', N lags = ' num2str(size(trialsData.xVals,1)) ', N expected = ' num2str(size(t_lag,1))]);
                    end
                end
                
                % update header info
                if dataStruct.loadHeaderInfo
                    if isempty(H_all{subj})
                        load(params.storage.cacheFile, 'H', 'selCh_H_resp');
                        assert(size(trialsData.yVals,2) == size(selCh_H_resp,2)); 
                        H.selCh = selCh_H_resp;
                        params.path2others.isarg_atlas = 'G:\shareData\visualization_normBrains\isarg_atlas_v8';
                        [ass_isargAtlas, H] = updateHeaders_atlasInfo(params, selCh_H_resp, H, params.storage.cacheFile);
                        H_all{subj} = H;
                    end
                end
                
                % significance & type of tuning
                if isfield(dataStruct, 'getSignificance')
                    if dataStruct.getSignificance
                        dataStruct.thisNeur = neurVars{neur};
                        dataStruct.thisTraj = trajVars{traj};
                        dataStruct.subjTag = subjList{subj,1};
                        
                        % significance of tuning values
                        [dataStruct.hVals{neur,traj,subj},dataStruct.sVals{neur,traj,subj}] = tuning_significance_get(trialsData, dataStruct);
                    end
                end

                if isfield(dataStruct, 'getTuningType')
                    if dataStruct.getTuningType                
                        % load AVG tuning values
                        data_info.var2show = 'avg';
                        trialsData_AVG = tuning_fig_AVG(data_info);
                        
                        % classify types of tuning (velocity-speed, position-distance)
                        [dataStruct.sType{neur,traj,subj}, dataStruct.typeNames{traj}] = tuning_significance_type(trialsData_AVG, dataStruct.hVals{neur,traj,subj}, dataStruct);
                    end
                end
                
            end       
        end
    end

    % subtract xVel - absVel & xPos - absPos
    if dataStruct.subtractSNR_x_abs
        for neur = 1:size(neurVars,2)   
            % subtract xVel - absVel
            tf_xvel = false; tf_absVel = false;
            if ismember('xvel', trajVars)
                [tf_xvel,i_xvel] = ismember('xvel', trajVars);
                [tf_absVel,i_avel] = ismember('absVel', trajVars);
            elseif ismember('vel_2D', trajVars)
                [tf_xvel,i_xvel] = ismember('vel_2D', trajVars);
                [tf_absVel,i_avel] = ismember('absVel_2D', trajVars);
            end                
            if tf_xvel && tf_absVel
                for subj = 1:size(subjList,1)
                    d_nts{neur,i_xvel,subj} = d_nts{neur,i_xvel,subj} - d_nts{neur,i_avel,subj};
                end
            end

            % subtract xPos - absPos
            tf_xpos = false; tf_absPos = false;
            if ismember('xpos', trajVars)
                [tf_xpos,i_xpos] = ismember('xpos', trajVars);
                [tf_absPos,i_apos] = ismember('absPos', trajVars); 
            elseif ismember('pos_2D', trajVars)
                [tf_xpos,i_xpos] = ismember('pos_2D', trajVars);
                [tf_absPos,i_apos] = ismember('absPos_2D', trajVars);
            end                       
            if tf_xpos && tf_absPos
                for subj = 1:size(subjList,1)
                    % hack, tbd
                    if size(d_nts{neur,i_xpos,subj},2) ~= size(d_nts{neur,i_apos,subj},2)
                        nCh = size(d_nts{neur,i_apos,subj},2);
                        d_nts{neur,i_xpos,subj} = d_nts{neur,i_xpos,subj}(:,nCh) - d_nts{neur,i_apos,subj};
                    else
                        d_nts{neur,i_xpos,subj} = d_nts{neur,i_xpos,subj} - d_nts{neur,i_apos,subj};
                    end
%                     d_nts{neur,i_xpos,subj} = d_nts{neur,i_xpos,subj} - d_nts{neur,i_apos,subj};
                end   
            end
            
            % subtract xAcc - absAcc
            tf_xacc = false; tf_absAcc = false;
            if ismember('xacc', trajVars)
                [tf_xacc,i_xacc] = ismember('xacc', trajVars);
                [tf_absAcc,i_aacc] = ismember('absAcc', trajVars); 
            elseif ismember('acc_2D', trajVars)
                [tf_xacc,i_xacc] = ismember('acc_2D', trajVars);
                [tf_absAcc,i_aacc] = ismember('absAcc_2D', trajVars); 
            end                      
            if tf_xacc && tf_absAcc
                for subj = 1:size(subjList,1)
                    d_nts{neur,i_xacc,subj} = d_nts{neur,i_xacc,subj} - d_nts{neur,i_aacc,subj};
                end
            end
            
        end    
    end
end

%% DA: load over all neurs, trajs & subjects
if strfind(var2show, 'da')
    for neur = 1:size(neurVars,2)   
        display(['loading decoding results for: ' neurVars{neur} ' ...']);
        for traj = 1:size(trajVars,2)
            for subj = 1:size(subjList,1)
                
                % DA: single-subject & single-trajectory
                data_info = struct;
                data_info.trajVar = trajVars{traj};
                data_info.neurVar = neurVars{neur};
                data_info.pathBeg = pathBeg;
                data_info.figDir = figDir;
                data_info.subjTag = subjList{subj,1};
                data_info.plotFigs = false;
                [trialsData,params] = decoding_fig_AVG(data_info);
                d_nts{neur,traj,subj} = trialsData.yVals;
                
                % lag values
                if isempty(t_lag)
                    t_lag = trialsData.xVals;
                else
                    assert(size(t_lag,1) == size(trialsData.xVals,1));  % should be same lags for all 
                end

                % header info
                if dataStruct.loadHeaderInfo
                    if isempty(H_all{subj})
                        load(params.storage.cacheFile, 'H', 'selCh_H_pred');
                        assert(size(trialsData.yVals,2) == size(selCh_H_pred,2)); 
                        H.selCh = selCh_H_pred;
                        H_all{subj} = H;
                    end
                end
                
%                 % from other script: cat DA over trajs
%                 data_traj_avg = cat(3, data_traj_avg, trialsData.yVals);    % cat over trajs
%                 data_traj_sem = cat(3, data_traj_sem, trialsData.yErrs);    % cat over trajs
%                 text_traj = [text_traj, '\color[rgb]{' num2str(clrs(traj,1)) ' ' num2str(clrs(traj,2)) ' ' num2str(clrs(traj,3)) '}' trajVars{traj} ' '];
                
            end       
        end
    end
end
        
%% MI (mutual information): load over all neurs, trajs & subjects
if strcmp(var2show, 'MI')
    for neur = 1:size(neurVars,2)   
        display(['loading mutual information results for: ' neurVars{neur} ' ...']);
        for traj = 1:size(trajVars,2)
            for subj = 1:size(subjList,1)

                % SNR: single-subject & single-trajectory
                data_info = struct;
                data_info.trajVar = trajVars{traj};
                data_info.neurVar = neurVars{neur};
                data_info.pathBeg = pathBeg;
                data_info.figDir = figDir;
                data_info.subjTag = subjList{subj,1};
                data_info.var2show = var2show;
                %data_info.selLags = [-1, 1];
                data_info.plotFigs = false;
                [trialsData,params] = fig_mi_one(data_info);
                d_nts{neur,traj,subj} = trialsData.yVals;
                
                % lag values
                if isempty(t_lag)
                    t_lag = trialsData.xVals;
                else
                    %assert(size(t_lag,1) == size(trialsData.xVals,1));  % should be same lags for all
                    if size(t_lag,1) ~= size(trialsData.xVals,1)
                        display(['WARNING: ' subjList{subj} ', srate = ' num2str(params.srate) ', N lags = ' num2str(size(trialsData.xVals,1)) ', N expected = ' num2str(size(t_lag,1))]);
                    end
                end
                
                % header info
                if dataStruct.loadHeaderInfo
                    if isempty(H_all{subj})
                        load(params.storage.cacheFile, 'H', 'selCh_H_resp');
                        assert(size(trialsData.yVals,2) == size(selCh_H_resp,2)); 
                        H.selCh = selCh_H_resp;
                        H_all{subj} = H;
                    end
                end

            end       
        end
    end
end

%% lag values & channel info
dataStruct.t_lag = t_lag;
if dataStruct.loadHeaderInfo
    dataStruct.H_all = H_all;
end

display('loading done.');
