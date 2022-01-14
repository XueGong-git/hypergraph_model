% LinearHypergraphEmbedding  Embed nodes in hypergraph to a linear
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
x_est = V(:,2);

%{
eigenvalues = diag(lambda);
scatter(1:min(25,size(L,2)), eigenvalues(1:min(25,size(L,2))), 100, 'MarkerFaceColor','black');
xlabel('i','FontSize', 13);
ylabel('Eigenvalue \lambda_i','FontSize', 13);
set(gca,'fontsize',30);
ax = gca;
exportgraphics(ax,strcat('plots/linear_L_eigenvalues_c3=',num2str(round(c3,1)),'.eps'),'Resolution',300) 


% plot top eigenvectors
t = tiledlayout(4,1);
ax1 = nexttile;
plot(V(:,1), 'Color', 'black');
ylabel('v_0')
ax2 = nexttile;
plot(V(:,2), 'Color', 'black');
ylabel('v_1')
ax3 = nexttile;
plot(V(:,3), 'Color', 'black');
ylabel('v_2')
ax4 = nexttile;
plot(V(:,4), 'Color', 'black');
ylabel('v_3')
% Link the axes
linkaxes([ax1,ax2,ax3,ax4],'x');
linkaxes([ax1,ax2,ax3,ax4],'y');
% Add shared title and axis labels
title(t,'Top Eigenvectors')
% Move plots closer together
xticklabels(ax1,{})
t.TileSpacing = 'compact';
axis([0 size(W2,2) -1 1])
exportgraphics(t,strcat('plots/linear_eigenvectors_c3=',num2str(round(c3,1)),'.eps'),'Resolution',300) 

figure
%}

% return embedding

end
