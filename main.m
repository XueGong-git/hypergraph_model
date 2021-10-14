% generate_linear_hypergraph  Generate unweighted hypergraph that shows a linear
% structure
%
% INPUTS
% 
% - n   Number of nodes
% - x   Node embedding
% - gamma   Decay parameter for the linear hypergraph model
%
% OUTPUTS
% - (A2, A3)       Hypergraph -- with a matrix + a tensor

clear


n = 25;
K = 5; % number of clusters 
a = 1; % noise;
m = n/K; % number of node per cluster
I = zeros(n);
f = zeros(n);

c2 = 1/2; % weight of simple edge
c3 = 1/3; % weight of simple edge


for input = 1:3
    for gamma = 5
        
        switch input
    
        case 1 
            x = linspace(1,n,n);
            data_shape = "uniform";
        case 2 
            x = sort(repmat(linspace(n/K,n,K),1,m)+(2*a*rand(1,n)-a)); % trophic levels
            data_shape = "cluster";
        case 3
            x = linspace(0,0,n); %overlapping x
            data_shape = "overlap";
        end

        
        [W2, W3] = GenerateLinearHypergraph(x, gamma, c2, c3, data_shape);
        
        %shuffle input adjacency matrix
        idx_rand = randperm(size(W2,1));% shuffle the nodes
        W2 = W2(idx_rand,idx_rand);
        W3 = W3(idx_rand,idx_rand);
        
        [x_est] = LinearHypergraphEmbedding(W2, W3, gamma, c2, c3);
        
                
        
        
        % plot estimated embedding
        s = scatter(x(idx_rand), x_est, 200, 'MarkerFaceColor','black','MarkerEdgeColor','none');
        alpha(s,0.3) % transparent color
        xlabel('x','FontSize', 13);
        ylabel('x*','FontSize', 13);
        set(gca,'fontsize',30);
        ax = gca;
        exportgraphics(ax,strcat('plots/linear_hygraph_embedding_', data_shape,'_gamma=', num2str(round(gamma,2)),'.eps'),'Resolution',300) 

    end
end
