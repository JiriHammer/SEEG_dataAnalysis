%% plots transparent brains with highlighted neural networks from Yeo7 atlas

% (c) Jiri, Aug21

%% output directory
outDir = 'G:\dox\proj_switching_EI\brainNetworks_Yeo7N';
if ~exist(outDir, 'dir')
    mkdir(outDir);
end 

%% settings: paths to source data
% dataDir = 'F:\dox\ms5_deepDecoding\data\v9_realTime';

%% params settings
params_default;
% params.plot_brainTopo.plot_slices = false;
% params.plot_brainTopo.plot_projections = true;
params.plot_brainTopo = struct(...           % user interface structure: holds most (but not all!) of the user settings
    'time2plot', [-0.1:0.1:1.0], ...                    % in [s], w.r.t. 0 = cutting point (for example)
    'plot_slices', false, ...                            % 2D brain slices
    'plot_projections', false, ...                       % 3D brain model projections (axial, sagittal, coronal)
    'plot_animation', false, ...                        % 3D brain model GIF animation (takes longer time)
    'printResolution', 0, ...                           % choices: 0 (= screen resolution) or 600 (= dpi). Resolution of the figures.                          
    'colorMap', jet(128), ...                           % colormap for channel values
    'size_interpolate', 1.0, ...                        % in [mm], voxel size to which the brain is interpolated
    'size_coloredCube', 3.0, ...                        % in [mm], "voxel" size of the colored channel values
    'slicePlanes', {{'axial','sagittal','coronal'}},... % planes to plot, options: 'axial','sagittal','coronal'
    'model_views', []  ...           % camera view angles at which the snapshot of the 3D model are taken. For some odd reason (camera light?), view(0,90) makes grey background...
);

params.paradigm.usedParadigm = 'carDriving_motorSubjects';
params.paradigm.specificType = 'continuous turns';
[params.storage.pathBeg, params.storage.subjList] = get_subjectList(params);
% params.storage.subjList = subjList;

%% areas to plot
% anatomy areas in: G:\shareData\visualization_normBrains\areas_normalized2spm
mriFile_AnatAreas = {
    'wVisual', 'r'; ...
    'wSomatomotor', 'r'; ...
    'wDorsal Attention', 'r'; ...
    'wVentral Attention', 'r'; ...
    'wLimbic', 'r'; ...
    'wFrontoparietal', 'r'; ...
    'wDefault', 'r'; ...
    };     

%% --------------------- PLOT 3D brain topo -------------------------------
for nn = 1:size(mriFile_AnatAreas,1) ...:size(subjList,1)
    f = figure('units','normalized','outerposition',[0 0 1 1]);
    set(f, 'visible','on', 'Color','k');
    set(f, 'InvertHardcopy','off');                 % preserves black background
    %set(f, 'Colormap', clrmap.fig);
    set(gcf, 'Renderer','OpenGL');

    nRows = 1;
    nCols = 1;
    nPlot = 1;
    subjTag = '20_PR4';
    params.storage.subjTag = subjTag;
    
    % brain info
    plotInfo.plottingStyle = '3D_model';
    plotInfo.file2load = 'wc1T1_colin27';  % choices: subject-specific = 'wc1T1', normalized brain = 'wc1T1_colin27';
    plotInfo.brain = getBrainData(params, plotInfo);
    
    % L,R labels
    if nPlot == 1
        plotInfo.plot_text_orientation_LR = true;
    else
        plotInfo.plot_text_orientation_LR = false;
    end
    plotInfo.plot_text_tag = mriFile_AnatAreas{nn,1};
    
    plotInfo.aarea = cell(1,1);
    plotInfo.file2load = mriFile_AnatAreas{nn,1};      % choices: see: TO DO
    plotInfo.aarea{1} = getBrainData(params, plotInfo);
    plotInfo.aarea{1}.faceClr = mriFile_AnatAreas{nn,2};
    
    % channels MNI coors
     plotInfo.chnlsMNI = [];
    
    % axes
    ax = subplot(nRows,nCols,nPlot);
    if nPlot <= nCols
        set(gca,'Position', get(gca,'Position') + [0 -0.1 0 0]);
    else
        set(gca,'Position', get(gca,'Position') + [0 0.1 0 0]);
    end
    %plotInfo.axPos = ax.Position;  %  get(gca, 'Position');
    plotInfo.axHandle = ax;
    plotInfo.figHandle = f; ...get(gcf);
    plotInfo.visible_axis = false;    
    plotInfo.saveFigLater = true;
    nPlot = nPlot+1;
    
    % plot
    plotInfo.plot_chVals_asClrPatch = false;
    plotInfo.chnl_clims = [-1 1];
    plotInfo.circle_size = 20;
    plotInfo.fig_bkgClr = 'k';                          % choices: black ('k') or white ('w')
    plotInfo.text = ['subject: ' subjTag];   
    plotInfo.outDir = outDir;
    plotInfo.figName = ['topoBrain_' subjTag];
    plotInfo.viewAngle = [0 90];
    chnls_VAL = [];
    plot_brain3D_wrapper(params, chnls_VAL, plotInfo);
    1+1;
    
    % save
    figname = mriFile_AnatAreas{nn,1};
    fig_save(f, figname, outDir, 'format',{'png','fig'}, 'res',600);
    %close(f);  
end
