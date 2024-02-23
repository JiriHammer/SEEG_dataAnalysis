function str_short = getAbbreviation(str_long)
% return abbreviation (str_short) for long string (str_long) if found in the list below, 
% else returns str_long

% (c) Jiri, Oct23

if strcmp(str_long, 'Default')
    str_short = 'DMN';
elseif strcmp(str_long, 'Dorsal Attention')
    str_short = 'DAN';
elseif strcmp(str_long, 'Salience')
    str_short = 'SN';    
else
    str_short = str_long;
end
