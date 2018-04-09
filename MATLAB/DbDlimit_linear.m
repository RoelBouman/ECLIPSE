% Function that uses a linear threshold for gating cells based on the KDE 
% The linear threshold just set to 0 (then negative will be the responders, positive the healthy)

% INPUT: 
% - S_scores_PC_used
%%%%%%%%%%%%
% OUTPUT:
% - inoroutDbD: struct with 0 for cells with positive density(so in the control) and 1 for cells with negative density(so in the response)
% - perc_DbD: percentage of cells with negative resulting density, shown
%   per individuals
% - count_DbD: count of cells with negative resulting density, shown
%   per individuals
% 
%% Dependencies:
% ComputeWeightedDensities2D
% ProcessDensities
% LocalDensity
% DensDiffThreshold

%% Written by R. Folcarelli (function based on the script written by R. Bouman)

%% 
function [inoroutDbD, perc_DbD, count_DbD]=DbDlimit_linear(S_scores_PC_used)

threshold = 0;

inoroutDbD=S_scores_PC_used;
perc_DbD=S_scores_PC_used;
count_DbD=S_scores_PC_used;



[DensityControl, DensityResponse, X, Y, ~, ~, ~, ~ ] = ComputeWeightedDensities2D(S_scores_PC_used);

[DensityControl, DensityResponse, DensityDiff, ~, ~] = ProcessDensities(DensityControl, DensityResponse);



for l1 = 1:length(S_scores_PC_used)
    
    densities = LocalDensity(DensityDiff,S_scores_PC_used(l1).Data,X(1,:),Y(:,1));
    
    inoroutDbD(l1).Data = ~(DensDiffThreshold(densities,threshold)); % 1 for cells with negative density(so in the patient situation)
    
    perc_DbD(l1).Data=((sum(inoroutDbD(l1).Data))/(length(inoroutDbD(l1).Data)))*100;
    
    count_DbD(l1).Data=sum(inoroutDbD(l1).Data);
end


end
