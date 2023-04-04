
function [DomNetsFactors_all,DomNetsProbs_all,DomNetsClustValues] = ...
    clusteringNetworks(g,y,k,all_net_sim,HOSP_NETs_DIDs, AllHospID)
%% Cluster each domain networks using Matrix factorization
% g  # of domain_based Networks / # of nodes at the super network
% y  # Number of years
% k  # of clusters

% initializing Hospital networks (bottom) Clustering results 
%Symmetric Nonnegative Matrix factorization of all hospital networks over 4
%years
DomNetsFactors = cell(g,y);
% converted the matrix factorization into probabilities
DomNetsProbs = cell(g,y);
% the cluster and relative probibility for each hospital in every hospital
% network over 4 years
DomNetsClustValues = cell(g,y);

%% Step I : Cluster every  network in every year

for i = 1:g
    for j = 1:y
        A = all_net_sim{i,j};
        H = symnmf_newton(sparse(A), k);
        DomNetsFactors{i,j} = H;
        Dh = sum(H,2);
        Dh = diag(Dh.^(-1));
        H = Dh*H;
        DomNetsProbs{i,j} = H;
        [prob, H_clus] = max(H,[],2);
        DomNetsClustValues{i,j} = [prob, H_clus];
    end
end

% Incorporating all nodes in the hospital networks. add 0 for hospitals
% that are not presented in the network
DomNetsFactors_all = part2all(DomNetsFactors, HOSP_NETs_DIDs, AllHospID);
DomNetsProbs_all = part2all(DomNetsProbs, HOSP_NETs_DIDs, AllHospID);

