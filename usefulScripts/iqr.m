function Q = iqr(X,d)
% calculates the inter-quartile range of X along the dimension d

% (c) Jiri, May17

if size(X,d) == 1
    warning(['Computing inter-quartile range along a singleton dimension d = ' num2str(d)]);
end

Q = prctile(X,75,d) - prctile(X,25,d);
