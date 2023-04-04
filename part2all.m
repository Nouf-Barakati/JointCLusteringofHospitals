function all = part2all(part_data, par_IDs, All_IDs)
% a function to create networks data with all hopital ids with 0's on the
% factors if the hospital is not included in that network
% By Nouf Albarakati

h = 152;
g = 145;
y = 4;
k = 3;

all = cell(g,y);

for i = 1:g
    for j = 1:y
        p = par_IDs{i,j};
        pp = part_data{i,j}; % 
        [~,w] = size(pp);
        aa = zeros(h,w);
        [tf,~]=ismember(All_IDs,p);
        aa(tf,:) = pp;
        all{i,j} = aa;
    end
end

end