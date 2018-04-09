function [CH] = con_hul2(X, p)
%%
% Function used to predict the confidence hull based on a histogram. 
% Input:
% X is a 3-D matrix of m by n histogram by o samples
% p is the percentage of consistent pixels
%
% The function first sorts the histogram from high to low and cumaltively
% sums the vector t. The points used to create the confidence hull are 
% values of vector t below p.
%
% Output:
% CH is a matrix with the same size as the histogram with all the point of 
% the confidence hull
% CH_edge is the edge of the confidence hull
%
% Written by G.H. Tinnevelt 19 Septembre 2014
%%

if p >= 100 || p <= 0
    disp('Invalid percentage selected')
    return
end

M = zeros(size(X));
X_new = X;

for l1 = 1:size(X,3)
    x = X_new(:,:,l1);
    [t, idx] = sort(x(:), 'descend');
    t = cumsum(t);
    [aidx, bidx] = ind2sub(size(x), idx);
    for l2 = 1:find(t <= p/100, 1, 'last')
        M(aidx(l2), bidx(l2),l1) = 1;
    end
end
M_tot = sum(M,3)/size(X,3);
CH = zeros(size(X,1),size(X,2));
CH(M_tot >= 1-(p/100)) = 1;
% CH_edge = edge(CH);
% 
% [ind(:,2),ind(:,1)]=find(CH_edge==1) ; % Change X and Y
% ind(:,1)=edges(1,ind(:,1));
% ind(:,2)=edges(2,ind(:,2));
% CIpoints = ind;   

end    