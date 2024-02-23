%% plot 3D brain model: test script 
% - define 3D brain structures (should be coregistered & normalized to MNI space)
%       -> sees truct: params.plot_brain3D
% - define MNI coordinates and values to plot
%       -> chnls_MNI_VAL = 2D: [ch x 4], where 4 = 3 MNIs + 1 VAL
% - load brain volumes -> surfaces
%       -> brainVolumes_fixed.m
% - plot 3D brain model
%       -> brain3D_plot.m

% (c) Jiri, Nov21
% based on: plot_brain3D_patch.m & plot_brain3D_wrapper.m

%% output directory
outDir = 'G:\dox\proj_switching_EI\brainNetworks_Yeo7N';

%% params settings
params_default;
params.plot_brain3D = struct(...        % user interface structure: holds most (but not all!) of the user settings, see also inside brain3D_plot.m
    'plot_brain3D', true, ...           % 3D brain model
    'volumes2plot', {{...               % volumes to plot (coregistered .nii files)
        'icbm152', [1 1 1], 0.3;  ...     % 1x3 cell: name, color, transparency, brain choice: icbm152, wc1T1_colin27
...        'wDefault', [0 1 0], 0.2; ...
...        'wDorsal Attention', [1 1 0.5], 0.2; ...
        }}, ...          
    'chVals_asCircles', true, ...       % plots channel values as color circles
    'chVals_asPatches', true, ...       % each channel value is mapped onto surface (but may be overwritten by latter points!), see also distanceThreshold for focality inside brain3D_plot.m
    'chVals_colorMap', jet(256), ...    % colormap of channel values
    'chVals_patchTransp', 0.2, ...      % transparency of channel values
    'circleSizeLims', [5 60], ...       % in points, maps chnl values to circle size
    'size_interpolate', 1.0, ...        % in [mm], voxel size to which the brain is interpolated
    'backgroundColor', 'w', ...         % figure / axes background color, good choices: 'w', 'k'
    'text_LR', true, ...                % marks L / R side
    'text_tag', 'P1', ...    % string
    'viewAngle', [0, 90], ...           % axial = 0,90; sagittal = 100,10; coronal front = 0,180; coronal back = 0,0
    'visible_axis', false ...           % if axis & their labels are visible
); 

%% M1 (red) & S1 (blue) areas
% params.plot_brain3D = struct(...        % user interface structure: holds most (but not all!) of the user settings, see also inside brain3D_plot.m
%     'plot_brain3D', true, ...           % 3D brain model
%     'volumes2plot', {{...               % volumes to plot (coregistered .nii files)
%         'wc1T1_colin27', [1 1 1], 0.1;  ...     % 1x3 cell: name, color, transparency
%         'wMotor_4a', [1 0 0], 1.0; ...
%         'wMotor_4p', [1 0 0], 1.0; ...
%         'wPSC_3a', [0 0 1], 1.0; ...
%         'wPSC_3b', [0 0 1], 1.0; ...
%         'wPSC_1', [0 0 1], 1.0; ...
%         'wPSC_2', [0 0 1], 1.0; ...
%         }}, ...          
%     'chVals_asCircles', false, ...       % plots channel values as color circles
%     'chVals_asPatches', false, ...       % each channel value is mapped onto surface (but may be overwritten by latter points!), see also distanceThreshold for focality inside brain3D_plot.m
%     'chVals_colorMap', jet(256), ...    % colormap of channel values
%     'chVals_patchTransp', 0.2, ...      % transparency of channel values
%     'circleSizeLims', [5 60], ...       % in points, maps chnl values to circle size
%     'size_interpolate', 1.0, ...        % in [mm], voxel size to which the brain is interpolated
%     'backgroundColor', 'w', ...         % figure / axes background color, good choices: 'w', 'k'
%     'text_LR', true, ...                % marks L / R side
%     'text_tag', 'P1', ...    % string
%     'viewAngle', [0, 90], ...           % axial = 0,90; sagittal = 100,10; coronal front = 0,180; coronal back = 0,0
%     'visible_axis', false ...           % if axis & their labels are visible
% ); 

%% ============================ IMPLEMENTATION ============================

%% -------------------------- channels: MNI coors & vals ------------------

%% variables: traj, neur, subj
trajVars = {'absVel'};   
trajVars_legend = {'speed'};
neurVars = {'hiGamma'};
subjList = msDist_getSubjList('CG_motorSubj');
selTime = 0;

%% paths to data (v14)
pathBeg_tune = 'G:\dox\ms4_distance\data\kinVars_tuning\v14_dirNondir';  % tuning

%% tuning SNR: load iEEG SNR values (all channels)
disp('loading data ...');
var2show = 'snr_vals';   % choices: 'snr_vals','snr_rank'
dataStruct = struct;
dataStruct.var2show = var2show;
dataStruct.pathBeg = pathBeg_tune;
dataStruct.figDir = outDir;
dataStruct.trajVars = trajVars;
dataStruct.neurVars = neurVars;
dataStruct.subjList = subjList;
dataStruct.subtractSNR_x_abs = false;       % subtract xVel - absVel OR xPos - absPos
[d_tune, dataStruct] = loadPooledData(dataStruct);   % d_tune = 3D cell: neur x traj x subj

%% cat all channels: MNIs & VALs
ch_MNI = [];
ch_VAL = [];
for subj = 1:size(subjList,1)
    assert(size(d_tune{1,1,subj},2) == size(dataStruct.H_all{subj}.selCh,2));
    i_t = closestval(dataStruct.t_lag, selTime);
    for ch = 1:size(dataStruct.H_all{subj}.selCh,2)
        thisCh = dataStruct.H_all{subj}.selCh(ch);
        ch_MNI = cat(1, ch_MNI, ...
            [dataStruct.H_all{subj}.channels(thisCh).MNI_x, ...
             dataStruct.H_all{subj}.channels(thisCh).MNI_y, ...
             dataStruct.H_all{subj}.channels(thisCh).MNI_z]);
        ch_VAL = cat(1, ch_VAL, d_tune{1,1,subj}(i_t,ch)); 
    end
end
    
% channels MNI coors
chnls_MNI_VAL = [ch_MNI, ch_VAL];      % ch x 4, where 4 = 3 MNIs + 1 VAL

%% ------------------ brain volumes -> surfaces ---------------------------
[brainVols, clrmap, inds_vol] = brainVolumes_fixed(params, params.plot_brain3D);
% plotInfo.clrmap = clrmap;   % DEL ?

%% --------------------- PLOT 3D brain topo -------------------------------
f = figure('units','normalized','outerposition',[0 0 1 1]);
ax = subplot(1,1,1); 
hold on;

% ch-spec info
chVals_lims = [prctile(chnls_MNI_VAL(:,4),5),prctile(chnls_MNI_VAL(:,4),95)];

% plot
plotInfo = params.plot_brain3D;
plotInfo.axHandle = ax;
plotInfo.figHandle = f; ...get(gcf);
plotInfo.chVals_lims = chVals_lims;     % computed automatically if left out
for vol = 1:size(brainVols,2)
    if brainVols{vol}.loaded
        plotInfo.thisVolume = vol;
        
        if vol == 1
            plotInfo.chVals_asPatches = false;
        else
            plotInfo.chVals_asPatches = false;
        end
        
        % plot all hnl values at once
        h = brain3D_plot(brainVols{vol}, chnls_MNI_VAL, plotInfo);
        
        % plot one third of chnl values at a time
%         i_ch = (vol-1)*floor(size(chnls_MNI_VAL,1)/3) + [1:floor(size(chnls_MNI_VAL,1)/3)]; % select indices
%         h = brain3D_plot(brainVols{vol}, chnls_MNI_VAL(i_ch,:), plotInfo);        
    end
end
disp('done.');

%% save
figName = 'topoBrain';
fig_save(f, figName, outDir, 'format',{'png','fig'}, 'res',600);
%close(f);  

%%
fileName = ['G:\dox\proj_switching_EI\brainNetworks_Yeo7N\chnls_MNI_VAL.mat'];
save(fileName, 'chnls_MNI_VAL', '-v7.3');
disp('saved');
