%% Remove based on Q statistic/SPE


% Function that calculates the Sum Prediction Error for the control model with the
% specified number of PCs. 95% is chosen as limit for the SPE estimated on
% control individuals
% INPUT:
%  - S_scaled: struct contaning the data preprocessed
%  - S_scores: struct contaning the scores per individuals
%  - controlmodel: SCA model estimated on the control individuals
%  - pc_used: combinations of two PCs used for building the control model
%%%%%%%%%
% OUTPUT:
%  - inoroutSPE: struct contaning the data preprocessed
%  - perc_aboveSPE: percentage of cells lying outside the SPElimit, shown
%  per individual
%  - count_SPE: counting of cells lying outside the SPElimit, shown per
%  individual
%  - SPElimit95: value of the SPE calculated on the control individuals

% Written by R. Folcarelli at Radboud University



function [inoroutSPE, perc_aboveSPE, count_SPE, SPElimit95]=SPElimitfunc(S_scaled, S_scores, controlmodel, pc_used)

S_scores_PC_used=S_scores;
for l1=1:length(S_scores_PC_used)
    S_scores_PC_used(l1).Data=S_scores(l1).Data(:,pc_used);
end

inoroutSPE=S_scores_PC_used;
perc_aboveSPE=S_scores_PC_used;
count_SPE=S_scores_PC_used;

%Calculate error matrix
S_error=S_scores_PC_used;
for l1=1:length(S_scores_PC_used)
    S_error(l1).Data=S_scaled(l1).Data-S_scores_PC_used(l1).Data*controlmodel.loadings(:,pc_used)';
end

%Evaluation of the residual error
S_SPE=S_error;
for l1 = 1:length(S_error)
    S_SPE(l1).Data = sum((S_error(l1).Data).^2,2);
end

SPE_control=vertcat(S_SPE(vertcat(S_scaled.Labels)==0)); %Calculation the 95% (or 99%) confidence interval upper limit mantaining the multiset structure
SPElimit95ind=SPE_control;
SPElimit99ind=SPE_control;

for i=1:length(SPE_control)
    SPElimit95ind(i).Data=(var(SPE_control(i).Data)/(2*mean(SPE_control(i).Data)))*chi2inv(.95,((2*mean(SPE_control(i).Data))^2)/var(SPE_control(i).Data));
    SPElimit99ind(i).Data=(var(SPE_control(i).Data)/(2*mean(SPE_control(i).Data)))*chi2inv(.99,((2*mean(SPE_control(i).Data))^2)/var(SPE_control(i).Data));
end

SPElimit95=mean(vertcat(SPElimit95ind.Data));
SPElimit99=mean(vertcat(SPElimit99ind.Data));

%Construct logical based on SPE limit.
for l1 = 1:length(S_SPE)
    if S_SPE(l1).Labels==1
        inoroutSPE(l1).Data=logical((S_SPE(l1).Data > SPElimit95));% cells above the limit are 1
    elseif S_SPE(l1).Labels==0
        inoroutSPE(l1).Data=logical(S_SPE(l1).Data > SPElimit99); % cells above the limit are 1
    end
end


for l1 = 1:length(perc_aboveSPE)
    perc_aboveSPE(l1).Data=(sum(inoroutSPE(l1).Data))/(size(inoroutSPE(l1).Data,1))*100;
end


for l1 = 1:length(count_SPE)
    count_SPE(l1).Data=sum(inoroutSPE(l1).Data);
end


end
