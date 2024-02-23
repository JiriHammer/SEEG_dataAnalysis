function str_out = str_splitSortDelWM(str_in, delimiter)
% splits string to cell array based on delimiter '/'
% sorts alphabetically
% removes any occurances of 'wm'
% cats again to string

% (c) Jiri, Aug22

if nargin < 2
    delimiter = '/';    % default
end

if isempty(strfind(str_in, '/'))
    str_out = str_in;   % no delimiter found
    return;
end

% splits string to cell array based on delimiter '/'
str_split = split(str_in,delimiter);

% sorts alphabetically
str_sort = sort(str_split);

% remove 'wm'
str_out = [];
for s = 1:size(str_sort,1)
    if ~ismember(str_sort{s}, {'wm','WM','het','Het'})
        str_out = [str_out str_sort{s} delimiter];  % include delimiter
    end
end

% exclude final delimiter
str_out = str_out(1:end-1);
