function VAR = individual_VAR(S_auto, S_scores, LDS, pc_used)
%%
% Function used to calculate the variance per individual
%
% Input:
% S_auto    The struct containing the autoscaled data
% S_scores  The struct containing the scores data
% LDS       The Loadings
% pc_used
%
% Output:
% VAR       A m samples by n PCs double containing the variance explained
%
%
% Written by G.H. Tinnevelt and M. Kokla on 22-5-2015 at Radboud University
%%

N_PCs = length(pc_used);
VAR = zeros(length(S_auto), N_PCs);
% E = zeros(N_PCs, 1);
% Scores_data = vertcat(S_scores.Data);
for l2 = 1:N_PCs
    for l1 = 1:length(S_auto)
        VAR(l1,l2) = 1 - sum(var(S_auto(l1).Data - S_scores(l1).Data(:,l2)*LDS(:,pc_used(l2))'))/sum(var(S_auto(l1).Data));
    end
end

VAR = round(VAR*100,1);
end
