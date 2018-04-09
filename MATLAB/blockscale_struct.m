function [X] = blockscale_struct(S)
%%
%
% Function collects data from struct and scales each measurements according to the
% number of cells. 
%
% Input:
% struct S with i being the number of measurements and with
% the following fields:
% Data      The preprocessed data with m x n with m being the number of cells
% measured in that measurement and n being the number of surface proteins
% measured
% Labels    The labels with a label for each i
% ID        The ID labels with ID for each i
% Train     A logical vector with length i, with 1 belonging to the
% trainset and 0 belonging to the testset
%
%
% Output:
% The training data X which is scaled with the number of cells of each
% measurement
%
% Written by G.H. Tinnevelt on 12-October-2015 at Radboud University Nijmegen
% Edited by G.H. Tinnevelt to improve speed at Radboud University 24 May
% 2017
%%

if isfield(S(1), 'train')
    trainset = logical(vertcat(S.train));
else
    trainset = true(length(S),1);
end

S = S(trainset); 

for l1 = 1:length(S)
    S(l1).Data = S(l1).Data/sqrt(size(S(l1).Data,1));
end

X = vertcat(S.Data); 
end
