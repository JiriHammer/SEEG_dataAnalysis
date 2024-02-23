function label_AA = msSEI_get_label_AA(brainAtlas)
% returns labels of anatomic area atlas (e.g. Yeo7)

% (c) Jiri, Nov23

label_AA_yeo7 = {
    'Visual',           'VIS';
    'Dorsal Attention', 'DAN';
    'Somatomotor',      'SMT';
    'Ventral Attention','VAN';
    'Frontoparietal',   'FPN';
    'Limbic',           'LIM';
    'Default',          'DMN';
    };
label_AA_yeo17 = {
    'Visual peripheral',    'VIS P';
    'Visual central',       'VIS C';
    'Dorsal attention A',   'DAN A';
    'Dorsal attention B',   'DAN B';
    'Somatomotor A',        'SMT A';
    'Somatomotor B',        'SMT B';
    'Ventral attention',    'VAN';
    'Salience',             'SAL';
    'Control C',            'FPN C';
    'Control A',            'FPN A';
    'Control B',            'FPN B';
    'Limbic N9',            'LIM N9';
    'Limbic N10',           'LIM N10';
    'Default D (Auditory)', 'DMN D';
    'Default C',            'DMN C';
    'Default A',            'DMN A';
    'Default B',            'DMN B';
    };

assert(exist(['label_AA_' lower(brainAtlas)],'var') == 1);
label_AA = eval(['label_AA_' lower(brainAtlas)]);
