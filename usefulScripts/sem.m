function y = sem(X, d)
% computes standard error of the mean (SEM) of matrix 'X' along the dimension 'd'

% (c) Jiri, Mar17

y = nanstd(X, 0, d)./size(X,d)^0.5;
