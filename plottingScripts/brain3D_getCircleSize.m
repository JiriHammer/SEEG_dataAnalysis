function def_circle_size = brain3D_getCircleSize(chnls_MNI_VAL, chVals_lims, circleSizeLims)
% for example, assume circleSizeLims = [5, 30]

% (c) Jiri, Sep22

%% (1) linearly scale to chVals_lims from [-30 to 30]
def_circle_size = linTransform(chnls_MNI_VAL(:,4), chVals_lims, [-circleSizeLims(2), circleSizeLims(2)]);

%% (2) apply absolute value, so that small & high values -> large circles
def_circle_size = abs(def_circle_size);   

%% (3) clip close-to zero vals to min. circle size = 5                      
def_circle_size(def_circle_size<circleSizeLims(1)) = circleSizeLims(1);
