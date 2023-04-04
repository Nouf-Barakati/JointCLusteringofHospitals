function [DomNetsSim,DomNetsSimSparse] = ...
    graph_similarities_Cosine(DomNetsProbs_all,g,y,k,s)
%% Step II : claculate the similarity among clustering structures of all Domain Networks

% Calculating the network Similarity for each cluster
DomNetsSim = cell(y,1);
%DomNetsClustSim = cell(y,1);

for j = 1:y
    sim_matrix = zeros(g,g);
    %sim_matrix_clust = cell(k,1);
    %for l= 1:k
        %sim_matrix_clust{l}(:,:) = zeros(g,g);
    %end
    
    for i = 1:g  % first graph
        for ii = i:g  % second graph
            % finding the nodes present in both graphs
            %x = sum(DomNetsProbs_all{i,j},2) + sum(DomNetsProbs_all{ii,j},2) ==2 ;
            net_perm = perms(1:k);
     
            a = DomNetsProbs_all{i,j};%(x,:);
            cos_sim=0;
            for p = 1:length(net_perm)
                b = DomNetsProbs_all{ii,j}(:,net_perm(p,:));
                cos_sim(p,1) = dot(a,b)/(vecnorm(a).*vecnorm(b));
            end
            
            %cos_sim = sum(a.*b)./(sqrt(sum(a.^2)).*sqrt(sum(b.^2)));
            sim_matrix(i,ii) = max(cos_sim);% *not_comm_nodes);
            % Try all permutation of the second graph - get the maximun 
            
            %for l = 1:k
             %   sim_matrix_clust{l}(i,ii) = cos_sim(l);%*not_comm_nodes;  
            %end   
        end
    end
    sim_matrix = sim_matrix+transpose(triu(sim_matrix,1)); % copy upper triangle to bottom 
    DomNetsSim(j) = {sim_matrix};
    %DomNetsClustSim(j) = {sim_matrix_clust};
end

% sparse 

DomNetsSimSparse = cell(y,1);
%DomNetsClustSimSparse = cell(y,1);
for i=1:4
    DomNetsSimSparse(i) = {DomNetsSim{i}.*(DomNetsSim{i}>s)};
    %for l = 1:k
        %DomNetsClustSimSparse{i}{l}(:,:) = DomNetsClustSim{i}{l}(:,:).*(DomNetsClustSim{i}{l}(:,:)>s);
    %end
end
