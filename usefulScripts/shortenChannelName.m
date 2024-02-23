function shortStr = shortenChannelName(longStr, maxLngth)
% shortens the anatomical assignemnt strings
% 'maxLngth' = maximum number of characters in 1 word

% (c) Jiri, Sep16

%% input
if nargin < 2
    maxLngth = 4;
end
tmp = longStr;

%% exchange Left (=L)
str2ex = 'Left';
k = strfind(tmp, str2ex);
if ~isempty(k)
    tmp(k:k+(length(str2ex))) = [];
    tmp = ['L ' tmp];
end

%% exchange Right (=R)
str2ex = 'Right';
k = strfind(tmp, str2ex);
if ~isempty(k)
    tmp(k:k+(length(str2ex))) = [];
    tmp = ['R ' tmp];
end

%% delete parenthesis at the end
i_parenthesis = strfind(tmp, '(');
if ~isempty(i_parenthesis)
    tmp(i_parenthesis:end) = [];    
end

%% truncate wrods longer than 4 characters to 3 chars
tmp = [' ' tmp ' '];
i_space = [strfind(tmp,' ')];
i_del = [];
for k = 2:size(i_space,2)
    if i_space(k) - i_space(k-1) > maxLngth
        i_del = cat(2, i_del, i_space(k-1)+maxLngth:i_space(k)-1);
    end
end
tmp(i_del) = [];

%% output
shortStr = tmp;
