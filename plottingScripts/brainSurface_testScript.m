%% iEEG activity to .nii

% (c) Jiri, May22

%% load activations in MNI -> chnls_MNI_VAL
% variables: traj, neur, subj
trajVars = {'absVel'};   
trajVars_legend = {'speed'};
neurVars = {'hiGamma'};
subjList = msDist_getSubjList('CG_motorSubj');
selTime = 0;
outDir = 'G:\dox\proj_switching_EI\brainNetworks_Yeo7N';

% paths to data (v14)
pathBeg_tune = 'G:\dox\ms4_distance\data\kinVars_tuning\v14_dirNondir';  % tuning

% tuning SNR: load iEEG SNR values (all channels)
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

% cat all channels: MNIs & VALs
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

%% load .nii example -> VOL
fileName = 'C:\Users\Jiri\Documents\MATLAB\mni2fs\examples\AudMean.nii';
hdr=spm_vol(fileName);
[VOL, xyz] = spm_read_vols(hdr);

%% MNI axis -> x,y,z
tmp = reshape(xyz(1,:),size(VOL));
x = -squeeze(tmp(:,1,1));
x = sort(x);                    % force ascending order (01.08.2018)
assert(x(1) < x(end));          % assure ascending order
voxSize_old(1) = x(2) - x(1);

tmp = reshape(xyz(2,:),size(VOL));
y = squeeze(tmp(1,:,1))';
assert(y(1) < y(end));
voxSize_old(2) = y(2) - y(1);

tmp = reshape(xyz(3,:),size(VOL));
z = squeeze(tmp(1,1,:));
assert(z(1) < z(end));
voxSize_old(3) = z(2) - z(1);

%% modify its volume information
V = zeros(size(VOL));                           % init
for n = 1:size(chnls_MNI_VAL,1)
    ix = closestval(x, chnls_MNI_VAL(n,1));     % index of x-coor
    iy = closestval(y, chnls_MNI_VAL(n,2));
    iz = closestval(z, chnls_MNI_VAL(n,3));
    V(ix,iy,iz) = chnls_MNI_VAL(n,4);           % value at voxel
end
C_1 = ones(5,5,5);     % conv. kernel, 1 cubic cm, 1 cm = 5 * 2 mm (voxel size)
V_conv = convn(V,C_1,'same');

C_2 = ones(10,10,10);     % conv. kernel, 2 cubic cm, 2 cm = 10 * 2 mm (voxel size)
V_conv_2 = convn(V,C_2,'same');

C_10 = ones(50,50,50);     % conv. kernel, 10 cubic cm, 2 cm = 50 * 2 mm (voxel size)
V_conv_10 = convn(V,C_10,'same');

%% save as new .nii file
fileName = 'C:\Users\Jiri\Documents\MATLAB\mni2fs\examples\myMotor_wiConv_1.nii';
V_new = V_conv;

% fileName = 'C:\Users\Jiri\Documents\MATLAB\mni2fs\examples\myMotor_wiConv_2.nii';
% V_new = V_conv_2;

% fileName = 'C:\Users\Jiri\Documents\MATLAB\mni2fs\examples\myMotor_wiConv_10.nii';
% V_new = V_conv_10;

% fileName = 'C:\Users\Jiri\Documents\MATLAB\mni2fs\examples\myMotor_noConv.nii';
% V_new = V;

hdr_new = hdr;
hdr_new.fname=fileName;

spm_create_vol(hdr_new);
spm_write_vol(hdr_new, V_new);
disp('done');

%% plot using mni2fs
figure;
toolboxpath = fileparts(which('mni2fs'));
mni2fs_auto(fullfile(toolboxpath, 'examples/myMotor_wiConv_1.nii'),'rh');
% mni2fs_auto(fullfile(toolboxpath, 'examples/myMotor_wiConv_2.nii'),'rh');
% mni2fs_auto(fullfile(toolboxpath, 'examples/myMotor_wiConv_10.nii'),'rh');

