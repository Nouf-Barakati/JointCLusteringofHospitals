function [normNet] = norm2Sim(net)
% min_max normalization then transforming to similarity (S = 1-d)
% input: network befor min max normalization) -Distance Matrix
% output: normalized Net
% By Nouf Alabarakati
NetMin = min(min(net));
NetMax = max(max(net));
normNet = 1- ((net-NetMin)/(NetMax-NetMin));       

end