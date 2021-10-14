% linear_hypergraph_embedding  Embed nodes in hypergraph to a linear
% position
% structure
%
% INPUTS
% 
% - W2   2nd order adjacency matrix
% - W3   3rd order adjacency tensor
% - gamma   Decay parameter for the linear hypergraph model
% - c2   Coefficient for 2nd order Laplacian
% - c3   Coefficient for 3rd order Laplacian

% OUTPUTS
% - x       Node embedding


%c2 = 1/2;
%c3 = 1/3;

function [x_est] = LinearHypergraphEmbedding(W2, W3, gamma, c2, c3)




% calculate degree matrix D2 and D3
d2 = sum(W2, 2); D2 = diag(d2);
d3 = sum(W3, 2); D3 = diag(d3);

% calculate graph Laplacian
L2 = D2 - W2;
L3 = D3 - W3;
L = c2*L2 + c3*L3;


% solve eigenvectors
[V, lambda] = eigs(L,size(L,1),'smallestabs'); % all eigenvectors ranging from smallest eigenvalue to largest eigenvalue

% plot eigenvalues
scatter(1:size(L,1), diag(lambda));

% return embedding
x_est = V(:,2);

end
