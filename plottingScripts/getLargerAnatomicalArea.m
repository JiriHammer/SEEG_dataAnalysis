function aaName = getLargerAnatomicalArea(ass_cytarch, mode)
% based on the cytoarchitectonic assignment 'ass_cytarch' by SPM Anatomy toolbox,
% returns larger anatomical structure 'aaName'.
% Useful for getting more channels into anatomical areas.

% (c) Jiri, May17
% current version of SPM Anatomy toolbox: 22c or 18 (ECoG)

%% default mode
if nargin < 2
    mode = 'aa_one';
end

%% definition of anatomical areas: motor, sensory
% AA = {
% 'FrontalPole',  {'Area Fp1','Area Fp2'};
% 'BasalForebrain',{'BF (Ch 1-3)','BF (Ch 4)'};
% 'OFC',          {'Area Fo1', 'Area Fo2','Area Fo3'};
% 'Broca',        {'Area 44','Area 45'};
% 'Motor',        {'Area 4a','Area 4p',    'Area 6'};
% 'Sensory',      {'Area 1','Area 2','Area 3a','Area 3b'};
% 'Insula',       {'Area Id1 ','Area Ig2 ','Area Ig1 '};
% 'Cingulum',     {'Area 25','Area s24','Area s32','Area 33'};
% 'Operculum',    {'Area OP1 [SII]','Area OP2 [PIVC]','Area OP3 [VS]','Area OP4 [PV]'};
% 'Auditory',     {'Area TE 1.0','Area TE 1.1','Area TE 1.2','Area TE 3',    'TE 3'};
% 'Amygdala',     {'Amygdala (CM)','Amygdala (SF)','Amygdala (AStr)','Amygdala (LB)'};
% 'Hippocampus',  {'DG (Hippocampus)','CA1 (Hippocampus)','CA2 (Hippocampus)','CA3 (Hippocampus)','HATA Region','Entorhinal Cortex','Subiculum'};
% 'IPL',          {'Area PFop (IPL)','Area PFm (IPL)','Area PF (IPL)','Area PFcm (IPL)','Area PGp (IPL)','Area PGa (IPL)','Area PFt (IPL)',     'IPC (PF)'};
% 'SPL/IPS',      {'Area 5M (SPL)','Area 5L (SPL)','Area 5Ci (SPL)','Area 7M (SPL)','Area 7A (SPL)','Area 7PC (SPL)','Area 7P (SPL)','SPL (7PC)',   'Area hIP1 (IPS)','Area hIP3 (IPS)','Area hIP2 (IPS)'};
% 'Visual',       {'Area FG1','Area FG2','Area  FG3','Area FG4','Area hOc1 [V1]','Area hOc2 [V2]','Area hOc4v [V4(v)]', 'Area hOc5 [V5/MT]','Area hOc3v [V3v]','Area hOc4d [V3A]','Area hOc4la','Area hOc3d [V3d]','Area hOc4lp'};
% 'Thalamus',     {'Thal: Prefrontal','Thal: Motor','Thal: Parietal','Thal: Premotor','Thal: Somatosensory','Thal: Visual','Thal: Temporal'};
% 'Cerebellum',   {'Lobule IX (Hem)','Lobule VI (Verm)','Lobule VIIa crusI (Verm)','Lobule VIIa crusII (Verm)','Lobule VIIIb (Verm)','Lobule VIIb (Hem)','Lobule VIIIa (Verm)','Lobule X (Verm)','Lobule VIIIb (Hem)','Lobule VIIb (Verm)','Lobule VIIIa (Hem)','Lobule IX (Verm)','Lobule VIIa crusI (Hem)','Lobule I IV (Hem)','Lobule VIIa crusII (Hem)','Lobule X (Hem)','Lobule VI (Hem)','Lobule V (Hem)','Ventral Dentate Nucleus','Dorsal Dentate Nucleus','Interposed Nucleus','Fastigii Nucleus'};
% };

%% definition of anatomical areas: motor+sensory
AA = {
'FrontalPole',  {'Area Fp1','Area Fp2'};
'BasalForebrain',{'BF (Ch 1-3)','BF (Ch 4)'};
'OFC',          {'Area Fo1', 'Area Fo2','Area Fo3'};
'Broca',        {'Area 44','Area 45'};
'MotorSensory', {'Area 4a','Area 4p',  ...                  % M1 (v22c)
                 'Area 6','aufCS','auf CS',...              % M1 (v18 & FR hierarchical assignment)
                 'Area 1','Area 2','Area 3a','Area 3b'};    % S1 
'Insula',       {'Area Id1 ','Area Ig2 ','Area Ig1 '};
'Cingulum',     {'Area 25','Area s24','Area s32','Area 33'};
'Operculum',    {'Area OP1 [SII]','Area OP2 [PIVC]','Area OP3 [VS]','Area OP4 [PV]'};
'Auditory',     {'Area TE 1.0','Area TE 1.1','Area TE 1.2','Area TE 3',    'TE 3'};
'Amygdala',     {'Amygdala (CM)','Amygdala (SF)','Amygdala (AStr)','Amygdala (LB)'};
'Hippocampus',  {'DG (Hippocampus)','CA1 (Hippocampus)','CA2 (Hippocampus)','CA3 (Hippocampus)','HATA Region','Entorhinal Cortex','Subiculum'};
'IPL',          {'Area PFop (IPL)','Area PFm (IPL)','Area PF (IPL)','Area PFcm (IPL)','Area PGp (IPL)','Area PGa (IPL)','Area PFt (IPL)',     'IPC (PF)'};
'SPL/IPS',      {'Area 5M (SPL)','Area 5L (SPL)','Area 5Ci (SPL)','Area 7M (SPL)','Area 7A (SPL)','Area 7PC (SPL)','Area 7P (SPL)','SPL (7PC)',   'Area hIP1 (IPS)','Area hIP3 (IPS)','Area hIP2 (IPS)'};
'Visual',       {'Area FG1','Area FG2','Area  FG3','Area FG4','Area hOc1 [V1]','Area hOc2 [V2]','Area hOc4v [V4(v)]', 'Area hOc5 [V5/MT]','Area hOc3v [V3v]','Area hOc4d [V3A]','Area hOc4la','Area hOc3d [V3d]','Area hOc4lp'};
'Thalamus',     {'Thal: Prefrontal','Thal: Motor','Thal: Parietal','Thal: Premotor','Thal: Somatosensory','Thal: Visual','Thal: Temporal'};
'Cerebellum',   {'Lobule IX (Hem)','Lobule VI (Verm)','Lobule VIIa crusI (Verm)','Lobule VIIa crusII (Verm)','Lobule VIIIb (Verm)','Lobule VIIb (Hem)','Lobule VIIIa (Verm)','Lobule X (Verm)','Lobule VIIIb (Hem)','Lobule VIIb (Verm)','Lobule VIIIa (Hem)','Lobule IX (Verm)','Lobule VIIa crusI (Hem)','Lobule I IV (Hem)','Lobule VIIa crusII (Hem)','Lobule X (Hem)','Lobule VI (Hem)','Lobule V (Hem)','Ventral Dentate Nucleus','Dorsal Dentate Nucleus','Interposed Nucleus','Fastigii Nucleus'};
};

%% list of larger (defined) anatomical areas
if strcmp(mode, 'list_largeAA')
    aaName = cell(size(AA,1),1);
    for a = 1:size(AA,1)
       aaName{a} = AA{a,1};
    end
end

%% individual cytoarchitectonic (anatomical) areas
if strcmp(mode, 'list_individualAA')
    % probabilitic MAP (loading)
    global MAP;
    %MapName = '/home/hammer/Documents/MATLAB/spm12/toolbox/Anatomy/Anatomy_v22c_MPM.mat';        % !!! edit here the correct path !!! (as of Apr2017)
    MapName = 'C:\Users\Jiri\Documents\MATLAB\spm12\toolbox\Anatomy\Anatomy_v22c_MPM.mat';        % !!! edit here the correct path !!! (as of Apr2017)
    assert(exist(MapName,'file') == 2);
    se_getMap('anat',MapName);

    % individual cytoarchitectonic (anatomical) areas
    aaName = cell(size(MAP,2),1);
    for a = 1:size(MAP,2)
        aaName{a} = MAP(a).name;
    end
end

%% asignment to larger (defined) anatomical areas
if strcmp(mode, 'aa_one')
    aaName = 'n.a.';        % default
    for a = 1:size(AA,1)
        if ismember(ass_cytarch, AA{a,2})
            aaName = AA{a,1};
        end
    end
end

