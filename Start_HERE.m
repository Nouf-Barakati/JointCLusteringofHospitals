clear 
% By Nouf ALbarakati Nouf@temple.edu
% Loading data
% HospNets: Hosp similarity networks 
% DisIDs 
% HospIDs

%% Initializations
g = 145; % # of Hospital_based Networks / # of nodes at the super network
y = 4;   % # of years
k = 3;   % # of clusters

%% Step I : Cluster every Hospital network in every year

% clustering Hospital Networks Seperatly using clusteringNetworks

% input: g: number of networks, y: number of years, k: number of clusters
% HospNets: a cell of all disease-specific hospital networks over 4 years {g X y}    
% HospNetID: the IDs of hospital included in the g X y networks
% HospID: All Hospital IDs included in the study

% Output: 
% DomNetsFactors_all: Symmetric Nonnegative Matrix factorization of all g
% hospital networks over y years
% DomNetsProbs_all: converted matrix factorization into probabilities of all g
% hospital networks over y years
% DomNetsClustValues: the cluster and relative probibility for each hospital in all hospital
% network over all year

[DomNetsFactors_all,DomNetsProbs_all,DomNetsClustValues] = ...
    clusteringNetworks(g,y,k,HospNets,HospNetIDs, HospID);

%% Step II : claculate the similarity among clustering structures of all Domain Networks
% Calculating the network Similarity for each cluster
% consider only nodes that are present in both networks 

% First: % Using Euclidean distance 
[DomNetsSim_euclidean,DomNetsSim_euclidean_sparse] = ...
    graph_similarities_Euclidean(DomNetsProbs_all,g,y,k,0.25);

% Second: % Using Cosine Similarity
[DomNetsSim_cosine,DomNetsSim_cosine_sparse] = ...
    graph_similarities_Cosine(DomNetsProbs_all,g,y,k,0.5);

%% Step III : Apply NoNClus algorithm to jointly cluster Hospital networks.
% joint clustering of hospital guided by
% model 1: Disease symptom similarity
% model 2: graph match similarity measured using euclidean distance 
% model 3: graph match similarity measured using cosine distance 

% NoNClus algorithm is used for the joint clustering
% by: Jingchao Ni (jingchao.ni@case.edu)
% Refernce : https://github.com/nijingchao/NoNClus

y=4;
g = 145;
k = 3; 
a =1;
t_u = 3*ones(g,1);
t_v = 3*ones(k,1);
MaxIter = 1000;
epsilon = 1e-6;

% MainNet: gXg : Disease Symptom similarity network 
MainNet1 = sparse(CCS_Sim_MeSH_red);
MainNet2 = cell(y,1);
MainNet3 = cell(y,1);

for i=1:y
    % MainNet: gXg : Disease similarity network extracted using Euclidean
    MainNet2(i,1) = {sparse(DomNetsSim_euclidean_sparse{i,1})};
    % MainNet: gXg : Disease similarity network extracted using Cosine
    MainNet3(i,1) = {sparse(DomNetsSim_cosine_sparse{i,1})};
end
% initializing NoN Clustering results 

NoN_KLD = cell(g,4);
NoN_KLD_labels = cell(g,4);
HospLabelAll_KLD = cell(g,4);
KLD_HospVal = cell(g,4);

% Model 1: using Disease Symptom Similarity Network (Mesh_Sim)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Running Phase I to share H data among all Domain networks
H = symnmf_newton(MainNet1, k);
Dh = sum(H,2);
Dh = diag(Dh.^(-1));
H = Dh*H;
[H_Vals_MeSH, H_idx_MeSH] = max(H,[],2);

%Running Phase II to share H data among all Domain networks
for i=1:4
    Us = NoNClus_phase2(DomKLDNets(:,i), DomKLDIDs(:,i), a, k, t_u, t_v, MaxIter, epsilon, H);
    NoN_KLD(:,i) = Us;
   [KLD_HospVal(:,i), NoN_KLD_labels(:,i),HospLabelAll_KLD(:,i)] = NoNClusLabelAllHosp(Us,DomKLDIDs(:,i), HospID); 
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Model 2: using Disease Similarity Network extracted using Euclidean
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NoN_KLD_Euclidean = cell(g,4);
NoN_KLD_labels_Euclidean = cell(g,4);
HospLabelAll_KLD_Euclidean = cell(g,4);
KLD_HospVal_Euclidean = cell(g,4);
H_Vals_Euclidean = cell(4,1);
H_idx_Euclidean = cell(4,1);

for i=1:4
    %Running Phase I to share H data among all Domain networks
    H = symnmf_newton(MainNet2{i,1}, k);
    Dh = sum(H,2);
    Dh = diag(Dh.^(-1));
    H = Dh*H;
    [Vals, H_idx] = max(H,[],2);
    H_Vals_Euclidean(i) = {Vals};
    H_idx_Euclidean(i) = {H_idx};

    %Running Phase II to share H data among all Domain networks
    Us = NoNClus_phase2(DomKLDNets(:,i), DomKLDIDs(:,i), a, k, t_u, t_v, MaxIter, epsilon, H);
    NoN_KLD_Euclidean(:,i) = Us;
   [KLD_HospVal_Euclidean(:,i), NoN_KLD_labels_Euclidean(:,i),HospLabelAll_KLD_Euclidean(:,i)] = ...
       NoNClusLabelAllHosp(Us,DomKLDIDs(:,i), HospID); 
end


% Model 3: using Disease Similarity Network extracted using Cosine
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NoN_KLD_cosine = cell(g,4);
NoN_KLD_labels_cosine = cell(g,4);
HospLabelAll_KLD_cosine = cell(g,4);
KLD_HospVal_cosine = cell(g,4);
H_Vals_cosine = cell(4,1);
H_idx_cosine = cell(4,1);

for i=1:4
    %Running Phase I to share H data among all Domain networks
    H = symnmf_newton(MainNet3{i,1}, k);
    Dh = sum(H,2);
    Dh = diag(Dh.^(-1));
    H = Dh*H;
    [Vals, H_idx] = max(H,[],2);
    H_Vals_cosine(i) = {Vals};
    H_idx_cosine(i) = {H_idx};

    %Running Phase II to share H data among all Domain networks
    Us = NoNClus_phase2(DomKLDNets(:,i), DomKLDIDs(:,i), a, k, t_u, t_v, MaxIter, epsilon, H);
    NoN_KLD_cosine(:,i) = Us;
   [KLD_HospVal_cosine(:,i), NoN_KLD_labels_cosine(:,i),HospLabelAll_KLD_cosine(:,i)] = ...
       NoNClusLabelAllHosp(Us,DomKLDIDs(:,i), HospID); 
end

