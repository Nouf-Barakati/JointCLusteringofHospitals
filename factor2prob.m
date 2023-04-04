function [prob, clust] = factor2prob(H)
%% find the probability to belong to one of the clusters
% by Nouf ALbarakati 
s = sum(H,2);
H =(diag(s.^(-1)))*H;
[prob, clust] = max(H,[],2);
end
