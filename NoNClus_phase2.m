function Us = NoNClus_phase2(DomNets, DomIDs, a, k, t_u, t_v, MaxIter, epsilon,H)

g = length(DomNets);

[Vals, H_idx] = max(H,[],2);

%% Phase II clustering

% Initialization of the domain node IDs in the hidden factor matrices

V_IDs = cell(k,1);

for i = 1:g
    
    V_IDs{H_idx(i)} = union(V_IDs{H_idx(i)}, DomIDs{i});
    
end

n_v = cellfun(@length,V_IDs);

% Generate mapping matrices O and D

O = cell(g,k);
D = cell(g,k);

for i = 1:g
    
    for j = 1:k
        
        [proj, IS, IV] = intersect(DomIDs{i},V_IDs{j});
        O{i,j} = sparse(IS,IV,ones(length(proj),1),length(DomIDs{i}),length(V_IDs{j}));
        D{i,j} = sparse(IS,IS,ones(length(proj),1),length(DomIDs{i}),length(DomIDs{i}));
        
    end
    
end

% Simultaneous clustering of the domain-specific networks
disp('Phase II : ');
tic;
Us = DomClus(DomNets, H, D, O, n_v, a, t_u, t_v, MaxIter, epsilon);
toc;
end