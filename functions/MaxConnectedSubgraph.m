function [xmax, W2max, W3max, T3max] = MaxConnectedSubgraph(x, c2, c3, W2, W3, T3)
% 
% return the largest connected component based on the hpergraph Laplacian c2*W2 + c3*W3


    % keep the largest connected graph
    G = graph(c2*W2 + c3*W3); % build combined adjacency matrix


    comp = conncomp(G);% find connected components
    [~,val] = max(histc(comp,unique(comp))); % find the largest component
    nodes = 1 : numnodes(G); 
    
    idx = comp==val;
    
    W2max = W2(idx,idx);
    W3max = W3(idx,idx);
    T3max = T3(idx,idx, idx);
    xmax = x(idx);
end