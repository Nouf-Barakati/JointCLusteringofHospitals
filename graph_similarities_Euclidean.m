function [DomNetsSim,DomNetsSimSparse] = ...
    graph_similarities_Euclidean(DomNetsProbs_all,g,y,k,s)
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
            %x = sum(DomNetsProbs_all{i,j},2) + sum(DomNetsProbs_all{ii,j},2) >= 1.9 ;
            % finding the nodes present in first graph
            %n1 = ((sum(DomNetsProbs_all{i,j},2) >=.9) & (sum(DomNetsProbs_all{ii,j},2) ==0));
            % finding the nodes present in second graph
            %n2 = ((sum(DomNetsProbs_all{i,j},2) ==0) & (sum(DomNetsProbs_all{ii,j},2) >=.9));
            % finding the nodes present in neither of the graphs
            %nn = sum(DomNetsProbs_all{i,j},2) + sum(DomNetsProbs_all{ii,j},2) ==0 ;
            
            %diff = sqrt(sum((DomNetsProbs_all{i,j}(x,:)-DomNetsProbs_all{ii,j}(x,:)).^2))/sum(x);%% ## model 1 ##
            net_perm = perms(1:k);
            a = DomNetsProbs_all{i,j};
            
            for p = 1:length(net_perm)
                b = DomNetsProbs_all{ii,j}(:,net_perm(p,:));
                euc_sim(p,:) =  sqrt(sum((a-b).^2));%% ## model 1 ##
            end
            
            
            sim_matrix(i,ii) = min(sum(euc_sim,2)); % *not_comm_nodes);
            
            %for l = 1:k
                %sim_matrix_clust{l}(i,ii) = diff(l);
            %end   
        end
    end
    sim_matrix = sim_matrix+transpose(triu(sim_matrix,1)); % copy upper triangle to bottom 
    sim_matrix = norm2Sim(sim_matrix);
    %for l= 1:k
        %sim_matrix_clust{l}(:,:) = norm2Sim(sim_matrix_clust{l}(:,:));
    %end
    
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
