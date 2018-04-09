function [mu] = mean_of_means(S)
%%

%%
mu = zeros(1,size(S(1).Data,2)); 
for l1 = 1:length(S)
    mu = mu + mean(S(l1).Data);
end
mu = mu/length(S);
