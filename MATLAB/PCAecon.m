%Function performs Principal Component Analysis ojn the imput data 
%Output:
%       - SCR(scores, which are the coordinates on the PC space).
%       - LDS(loading matrix, how much each variable contributes to a PC)
%       en VAR(variance explained per PC) 
%Input:
%data matrix

%% Written by R. Folcarelli
%%
function [SCR, LDS, VAR, D]=PCAecon(X)


[U,D,LDS]=svd(X,'econ');
SCR=U*D;
VAR=D.*D/sum(sum(D.*D));%fraction of variance.
reflect = diag(sign(sum(LDS.^3)));
SCR = SCR*reflect;
LDS = LDS*reflect;

VAR=diag(D.*D);
VAR = round(100*VAR/sum(VAR));

end 