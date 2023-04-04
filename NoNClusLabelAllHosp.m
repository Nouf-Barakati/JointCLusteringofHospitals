function [HospVal, HospLabel,HospLabelAll] = NoNClusLabelAllHosp(Us,DomID, HospID)

%%% Input
%
% Us: a set of cluster indicator matrices

%% Initialization
h = length(HospID);
g = length(Us);
HospLabel = cell(g,1);
HospVal = cell(g,1);
HospLabelAll = cell(g,1);

%% labeling clusters

for i = 1:g    
    matchAll = zeros(h,1);
    U_i = Us{i};
    Dh = sum(U_i,2);
    Dh = diag(Dh.^(-1));
    U_i = Dh*U_i;
    [U_i_Val, U_i_idx] = max(U_i,[],2);
    HospLabel(i) = {U_i_idx};   
    HospVal(i) = {U_i_Val};  
    [~,bi] =ismember(DomID{i},HospID);
    matchAll(bi,1) =U_i_idx;
    HospLabelAll(i) ={matchAll};
end


end