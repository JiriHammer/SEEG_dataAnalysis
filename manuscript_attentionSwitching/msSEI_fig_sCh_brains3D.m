%% single trials: brain + DMN or DAN + selCh
% plots 3D brain with highlighted NN + selCh

% (c) Jiri, Feb24

%% settings
sel_AA = {
    'Default';
    'Dorsal Attention';
};
sel_viewAngles =  {
%     'sagittal L', [-90, 0];          % sagittal R = 100,10 (90, 0)
%     'axial',      [0, 90];          % axial = 0,90    
    'sagittal R', [90, 0];         % sagittal L = -100,10 (-90, 0)
    };

%% settings
subjTag = '20_PR11';
selected_channels = {'Rd1-Rd2', 'Ld6-Ld7'}; % DMN, DAN !!! order must match: sel_AA !!!
selected_networks = {'DMN', 'DAN'};
pathBeg = 'F:\dox\ms_switch_EI\data\v31_stft_session_bip\switchin_EI_IE_bip';

%% load trials
% --- cache file
cacheFile = [pathBeg filesep subjTag filesep 'cacheFile.mat'];
assert(exist(cacheFile,'file') == 2);

% load H
clear H selCh_H_resp
load(cacheFile, 'H', 'selCh_H_resp');  
selCh_H = selCh_H_resp;

%% select channels
i_ch = [];
ch_MNI = [];
for c = 1:size(selected_channels,2)
    for ch = 1:size(selCh_H,2)
        if strcmp(H.channels(selCh_H(ch)).name, selected_channels{c})
            i_ch = cat(2, i_ch, ch);
            ch_MNI = cat(1, ch_MNI, [H.channels(selCh_H(ch)).MNI_x, H.channels(selCh_H(ch)).MNI_y, H.channels(selCh_H(ch)).MNI_z]);
        end
    end
end

%% figure
params_default;
for nn = 1:size(sel_AA,1) 
    nCols = size(sel_viewAngles,1);   % ~ freq
    nRows = 1;
    fig_W = 6;         % in [cm]
    fig_H = 6;    % in [cm], 3 cm per row
    f = fig_make('fig_W',fig_W, 'fig_H',fig_H);
    marg_h = [0.05 0.05];   % margin from [bottom, top]
    marg_w = [0.05 0.10];   % margin from [L, R] side
    gap = [-0.01, -0.01];     % between axes from [top, side]  

    %  ------------------------ PLOT: 3D brain for both NN ---------------
    %     marg_w = [0.02 0.10];   % margin from [L, R] side
    % settings
    params.plot_brain3D = struct(...        % user interface structure: holds most (but not all!) of the user settings, see also inside brain3D_plot.m
        'plot_brain3D', true, ...           % 3D brain model
        'volumes2plot', {{...               % volumes to plot (coregistered .nii files)
            'colin27_t1_tal_lin masked', [1 1 1], 0.20;  ...     % 1x3 cell: name, color, transparency (other choices for norm brain = wc1T1_colin27)
            }}, ...          
        'chVals_asCircles', true, ...       % plots channel values as color circles
        'chVals_asPatches', false, ...       % each channel value is mapped onto surface (but may be overwritten by latter points!), see also distanceThreshold for focality inside brain3D_plot.m
        'chVals_colorMap', [0 0 0], ...    % colormap of channel values
        'chVals_patchTransp', 0.2, ...      % transparency of channel values
        'circleSizeLims', [20 20], ...       % in points, maps chnl values to circle size
        'size_interpolate', 1.0, ...        % in [mm], voxel size to which the brain is interpolated
        'backgroundColor', 'w', ...         % figure / axes background color, good choices: 'w', 'k'
        'text_LR', false, ...                % marks L / R side
        'text_tag', '', ...    % string
        'viewAngle', [0, 90], ...           % axial = 0,90; sagittal = 100,10; coronal front = 0,180; coronal back = 0,0
        'visible_axis', false ...           % if axis & their labels are visible
    ); 

    % 3D brain: include Yeo7 NNs
    chnls_MNI_VAL = [];
    chnls_clrs_nn = [
        [139, 0, 0]./255;  % DMN = dark red
        [0, 80, 0]./255;  % DAN = dark green
    ];

%     [aa_found, aa] = ismember(sel_AA{nn}, label_AA(:,1));   % aa = index to GRUP !!!
%     assert(aa_found);
%     assert(strcmp(label_AA{aa,1}, sel_AA{nn}));

    % include Yeo7 NNs
    params.plot_brain3D.volumes2plot{end+1,1} = sel_AA{nn};
    params.plot_brain3D.volumes2plot{end,2} = msSwitch_get_clr4aa(sel_AA{nn});
    params.plot_brain3D.volumes2plot{end,3} = 0.2;

    % chnls MNI
%     tmp_ch = EC_info.ch_MNI{nn};               % sel chnls MNI
%     tmp_ch = cat(2, tmp_ch, nn*ones(size(tmp_ch,1),1)); % add value = nn (was significant)
%     chnls_MNI_VAL = cat(1, chnls_MNI_VAL, tmp_ch);      % add (again) the significant channels, (note that they are duplicate)
% %     params.plot_brain3D.chVals_colorMap = cat(1, params.plot_brain3D.chVals_colorMap, msSwitch_get_clr4aa(sel_AA{nn}));
%     params.plot_brain3D.chVals_colorMap = cat(1, params.plot_brain3D.chVals_colorMap, chnls_clrs_nn(nn,:));
    chnls_MNI_VAL = [ch_MNI(nn,:), 1];

    % brain volumes -> surfaces
    [brainVols, clrmap, inds_vol] = brainVolumes_fixed(params, params.plot_brain3D);

    % plot: 3D brains: different view angles
    for v = 1:size(sel_viewAngles,1)   

        % project channels MNI
        chnls_MNI_VAL_toPlot = chnls_MNI_VAL;
        if strcmp(sel_viewAngles{v,1}, 'sagittal L')
            chnls_MNI_VAL_toPlot(chnls_MNI_VAL(:,1)>=0,:) = [];  % only channels on L side
            chnls_MNI_VAL_toPlot(:,1) = min(brainVols{1}.xi);   % all channels on L side - up front!
        elseif strcmp(sel_viewAngles{v,1}, 'axial')
            chnls_MNI_VAL_toPlot(:,3) = max(brainVols{1}.zi);  % all channels on top - up front!
        elseif strcmp(sel_viewAngles{v,1}, 'sagittal R')
            chnls_MNI_VAL_toPlot(chnls_MNI_VAL(:,1)<0,:) = [];  % only channels on R side
            chnls_MNI_VAL_toPlot(:,1) = max(brainVols{1}.xi);   % all channels on R side - up front!
        end     

        % axes
    %     nPlot = fb + (v-1)*nCols;
        nPlot = v;
        ax = subtightplot(nRows, nCols, nPlot, gap, marg_h, marg_w);
        hold on;

        % set: plot_brain3D & plotInfo
        plotInfo = params.plot_brain3D;
        plotInfo.axHandle = ax;
        plotInfo.figHandle = f;
        plotInfo.chVals_lims = [0, 2];

        % PLOT 3D brain vols
        for vol = 1:size(brainVols,2)
            if brainVols{vol}.loaded
                plotInfo.thisVolume = vol;
                % plot all chnl values at once: uncomment line 137 in brain3D_plot.m
                h = brain3D_plot(brainVols{vol}, chnls_MNI_VAL_toPlot, plotInfo);   
            end
        end    

        % title: FB
    %     if v == 1
    %         title(label_FB{fb,3}, 'FontSize',fontSize);  % , 'FontWeight','normal'
    %     end

        % set view angle
        view(sel_viewAngles{v,2});
    end % view angles

    %% save
    figName = [selected_networks{nn} '_brain3D_v1'];
    outDir = [pathBeg filesep 'msSEI_fig_sCh_singleTrials' filesep 'selCh_brain3D'];
    fig_save(f, figName, outDir, 'format','png', 'res',600);
    % close(f); 
    
end % of nn
