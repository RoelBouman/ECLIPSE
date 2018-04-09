function [SCR, LDS, VAR]=pcafunction(X)
%%
%Function that perfoms PCA
% Input:
% - X data matrix
%% Output:
% - SCR scores
% - LDS loadings
% - VAR variance

% Written at Radboud University
%%

[U D V]=svd(X, 'econ');
SCR = U*D;
LDS = V;
VAR=diag(D.*D);
VAR = round(100*VAR/sum(VAR));
end