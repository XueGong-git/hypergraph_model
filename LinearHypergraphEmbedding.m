% linear_hypergraph_embedding  Embed nodes in hypergraph to a linear
% position
% structure
%
% INPUTS
% 
% - W2   2nd order adjacency matrix
% - W3   3rd order adjacency tensor
% - c2   Coefficient for 2nd order Laplacian
% - c3   Coefficient for 3rd order Laplacian

% OUTPUTS
% - x       Node embedding


%c2 = 1/2;
%c3 = 1/3;

function [x_est, V, lambda] = LinearHypergraphEmbedding(W2, W3, c2, c3, norm)




% calculate degree matrix D2 and D3
d2 = sum(W2, 2); D2 = diag(d2);
d3 = sum(W3, 2); D3 = diag(d3);

% calculate graph Laplacian
L2 = D2 - W2;
L3 = D3 - W3;
if norm == "true"
    d2_inv = 1./d2; d2_inv(isinf(d2_inv)|isnan(d2_inv)) = 0;
    d3_inv = 1./d3; d3_inv(isinf(d3_inv)|isnan(d3_inv)) = 0;
    L2 = diag(d2_inv)*L2;
    L3 = diag(d3_inv)*L3;
end
L = c2*L2 + c3*L3;


% solve eigenvectors
[V, lambda] = eigs(L,size(L,1),'smallestabs'); % all eigenvectors ranging from smallest eigenvalue to largest eigenvalue


% return embedding
x_est = V(:,2);

end
