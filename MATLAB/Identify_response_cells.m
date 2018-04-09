% Function to identify the abnormal/responding cells based on the SPE limit and the Difference between Densities limit
% % INPUT: 
% - S_original: struct with original data
% - inoroutSPE: struct with 0 and 1 for cells below and above the limit,
% respectevely 
% - inoroutDbD: struct with 0 and 1 for cells with positive and negative
% density from the Difference between Densities
%%%%%%%%%%%%
% OUTPUT:
% - S_in_original: struct containing original data, after the removal of 'normal' cells 
% - percentage_eclipse: percentage of cells retained in the ECLIPSE
% analysis

%% Written by R. Folcarelli at Radboud University

function [S_in_original, percentage_eclipse]=Identify_response_cells(S_original, inoroutSPE, inoroutDbD)

inorout=S_original;

for l1 = 1:length(S_original)
        inorout(l1).Data=(inoroutDbD(l1).Data+inoroutSPE(l1).Data);
        inorout(l1).Data = logical(inorout(l1).Data);
end

S_in_original = S_original;
for i = 1:length(S_in_original)
    S_in_original(i).Data = S_in_original(i).Data(inorout(i).Data,:);
end

% calculate the percentage of cells retained in the ECLIPSE model
percentage_eclipse=S_original;
for i = 1:length(S_in_original)
    percentage_eclipse(i).Data = (size(S_in_original(i).Data,1)/size(S_original(i).Data,1))*100;
    percentage_eclipse(i).count=(size(S_in_original(i).Data,1));
end
end


