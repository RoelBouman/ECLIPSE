function [Pr] = PRESSplot_control(S_scaled)
%%
% Function used to plot the PRESS
%
% Input
% S_scaled    Preprocess S struct
%
% Output
% Pr        1 x n vector with the PRESS with n being the number of PCs
%
%
% Written by G.H. Tinnevelt at Radboud University on 12 february 2016
%%

Labels = vertcat(S_scaled.Labels);
[X] = blockscale_struct(S_scaled(Labels == 0));

%% Press
[Pr]=pca_cv(X,[],0,round(size(X,1)/10)); 
figure
subplot(1,2,1)
plot(Pr)
title('PRESS', 'fontsize', 22)
xlabel('PCs', 'fontsize', 18)
ylabel('PRESS (a.u.)', 'fontsize', 18)
set(gca, 'Xtick', 1:size(X,2)-1, 'fontsize', 14)
%% Variance
[~, ~, Varp] = pcafunction(X); 
subplot(1,2,2)
plot(Varp)
title('Variance explained', 'fontsize', 22)
xlabel('PCs', 'fontsize', 18)
ylabel('Variance (%)', 'fontsize', 18)
set(gca, 'Xtick', 1:size(X,2), 'fontsize', 14)
axis([1 size(X,2) 0 100])
end