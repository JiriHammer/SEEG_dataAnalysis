function brain = brainVolumes_load(params, file2load, struct_plotBrain)
% loads 3D brain volume (MRI, CT, anatomic areas, ...) sepcified in file2load
% interpolates if different vosel size is specified
% returns brain = 
% SPM package MUST be included in Matlab path

% (c) Jiri, Nov21, Jan17
% based on: getBrainData.m

%% add SPM toolbox to paths
spm_dir = params.path2others.spm;
assert(exist(spm_dir,'dir') == 7);
addpath(spm_dir);

%% load normalized brain MRI (wT1.nii)
disp(['Loading brain file: ' file2load '.nii ...']);
if ~isempty(strfind(file2load, 'colin27')) || ~isempty(strfind(file2load, 'icbm152'))
    % normalized brain from colin27
    path2file = [params.storage.dir_neuroHammer filesep 'plot_3Dbrain' filesep 'normBrains_MRI'];
    separationThreshold = 0.5;                  % (value 0.5 separates gray matter from the dark background). Other values 0 - 1 may work also fine
elseif strfind(file2load, 'ch-spec')
    brain.loaded = false;
    return;
elseif ismember(file2load, {'wT1','wT2','wc1T1','wc2T1','wc3T1','wc4T1','wc5T1','wrCT'})
    % subject-specific brain
    % first, try header in the paradigm specific directory
    path2file = [params.storage.pathBeg filesep params.storage.subjTag filesep 'brainData' filesep 'output_spm'];
    if exist(path2file,'dir') ~= 7
        % second, try header in the header directory
        path2file =[params.storage.dir_shareData filesep 'subjHeaders' filesep params.storage.subjTag filesep 'brainData' filesep 'output_spm'];
        assert(exist(path2file,'file') == 7);
    end    
    separationThreshold = 0.5;                  % (value 0.5 separates gray matter from the dark background). Other values 0 - 1 may work also fine
else    
    % anatomy areas?
    atlasName = anatomicalAreas_atlas4AA(file2load);
    path2file = [params.path2others.normBrains filesep atlasName];
    separationThreshold = 0.1;                  % (value 0.5 separates gray matter from the dark background). Other values 0 - 1 may work also fine
end
assert(exist(path2file,'dir') == 7);
fileName = [path2file filesep file2load '.nii'];
if ~exist(fileName,'file')
    fileName = [path2file filesep file2load '.hdr'];
end
assert(exist(fileName,'file') == 2);
brain.hdr = spm_vol(fileName);
[brain.vol, brain.xyz] = spm_read_vols(brain.hdr);

%% MNI axis (in mm)
voxSize_old = [nan nan nan];
tmp = reshape(brain.xyz(1,:),size(brain.vol));
x = -squeeze(tmp(:,1,1));
x = sort(x);                    % force ascending order (01.08.2018)
assert(x(1) < x(end));          % assure ascending order
voxSize_old(1) = x(2) - x(1);

tmp = reshape(brain.xyz(2,:),size(brain.vol));
y = squeeze(tmp(1,:,1))';
assert(y(1) < y(end));
voxSize_old(2) = y(2) - y(1);

tmp = reshape(brain.xyz(3,:),size(brain.vol));
z = squeeze(tmp(1,1,:));
assert(z(1) < z(end));
voxSize_old(3) = z(2) - z(1);

%% interpolation (if needed)
voxSize_new = struct_plotBrain.size_interpolate;     % in [mm]
thr_smallNum = 0.001;
if abs(voxSize_new - voxSize_old(1)) < thr_smallNum && ...
   abs(voxSize_new - voxSize_old(2)) < thr_smallNum && ...
   abs(voxSize_new - voxSize_old(3)) < thr_smallNum
    VI = brain.vol;
    xi = x';
    yi = y';
    zi = z';
else
    [X,Y,Z] = meshgrid(x,y,z);
    xi = x(1):voxSize_new:x(end);
    yi = y(1):voxSize_new:y(end);
    zi = z(1):voxSize_new:z(end);
    [XI,YI,ZI] = meshgrid(xi,yi,zi);
    V = permute(brain.vol,[2,1,3]);                        % permutation needed by interp3 (for some reason)
    VI = interp3(X,Y,Z,V,XI,YI,ZI,'spline');
    VI = permute(VI,[2,1,3]);                                   % re-arrange back to [x,y,z] dimensions
end

%% return interpolated volume and axis
assert(size(VI,1) == size(xi,2));
assert(size(VI,2) == size(yi,2));
assert(size(VI,3) == size(zi,2));
brain.VI = VI;
brain.xi = xi;
brain.yi = yi;
brain.zi = zi;
brain.voxSize_new = voxSize_new;

%% isosurface for 3D model
if isfield(struct_plotBrain, 'plot_brain3D')
    disp('computing isosurface for 3D brain model ...');
    vLims = [prctile(VI(:),1), prctile(VI(:), 99)];
    if vLims(1) == vLims(2)
        vLims = [min(VI(:)), max(VI(:))];
        if vLims(1) < vLims(2)
            disp(['WARNING: only zero values found in: ' file2load]);
        end
    end
    V = linTransform(VI, vLims, [0, 1]);
    %brain.fv = isosurface(V, separationThreshold);    % surface, vertex 
    
    % try to add axes
    [XI,YI,ZI] = meshgrid(yi,xi,zi);
    brain.fv = isosurface(XI,YI,ZI,V, separationThreshold);    % surface, vertex 
end

%% was the brain structure loaded?
if isfield(struct_plotBrain, 'plot_brain3D')
    if isempty(brain.fv.vertices)
        brain.loaded = false;
    else
        brain.loaded = true;
    end
else
    brain.loaded = true;
end

