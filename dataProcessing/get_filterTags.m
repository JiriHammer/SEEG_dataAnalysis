function filterTags = get_filterTags(varName)
% return filter processing pipeline given the variable name
% typically designed for extraction of trajectory derivatives

% (c) Jiri, Apr17

%% notes
% - 'zeroMean' used until 19.07.2017 (large car game tuning analysis v5)

%% code
varName = lower(varName);
switch varName
    case 'xpos'                                 % 1-D task (car game)
        filterTags = {'raw'};
    case 'xvel'
        filterTags = {'derivative'};
    case 'xacc'
        filterTags = {'derivative','derivative'};      
    case 'abspos'
        filterTags = {'ampl'};        
    case 'absvel'
        filterTags = {'derivative','ampl'}; 
    case 'absacc'
        filterTags = {'derivative','derivative','ampl'};        
    case 'xvelPower'
        filterTags = {'derivative','n_power'}; 
    case 'xpos_dir'
        filterTags = {'raw'}; 
    case 'xvel_dir'
        filterTags = {'derivative'}; 
    case 'dddt_pos'
        filterTags = {'ampl','derivative','zeroOut_negative'};     % derivative of distance      
    case 'dddt_neg'
        filterTags = {'ampl','derivative','zeroOut_positive'};     % derivative of distance     
    case 'abs_dddt_pos'
        filterTags = {'ampl','derivative','zeroOut_negative','ampl'};     % derivative of distance      
    case 'abs_dddt_neg'
        filterTags = {'ampl','derivative','zeroOut_positive','ampl'};     % derivative of distance     
        
    case {'xpos_2d','ypos_2d'}                  % 2-D task (RTP, UBoat)
        filterTags = {'channelSelection'};
    case {'xvel_2d','yvel_2d'}
        filterTags = {'channelSelection','derivative'}; 
    case 'abspos_2d'
        filterTags = {'combine2complex','ampl'};         
    case 'absvel_2d'
        filterTags = {'derivative','combine2complex','ampl'}; 
    case 'pos_2d'
        filterTags = {'combine2complex'};                 
    case 'vel_2d'
        filterTags = {'derivative','combine2complex'};         
    case 'abs_dddt_pos_2d'
        filterTags = {'combine2complex','ampl','derivative','zeroOut_negative','ampl'};     % derivative of distance      
    case 'abs_dddt_neg_2d'
        filterTags = {'combine2complex','ampl','derivative','zeroOut_positive','ampl'};     % derivative of distance     
        
    otherwise
       error(['unknown variable name: ' varName]);
end
