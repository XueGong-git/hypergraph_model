tic
clear


c2 = 1/10; % weight of simple edge
c3_array = [1, 3, 5, 7]/30;
gamma_array = 0:1:15; % gamma for likelihood plot
data_type = "highschool";
K = 9;
rand_linear = [];
rand_periodic = [];


if data_type == "highschool"
    [W2, W3, T3, E2, E3] = LoadHighSchool();
    label = readtable('highschool_label', 'ReadVariableNames', false);
elseif data_type == "senate_bill"
    [W2, W3, T3, E2, E3] = LoadSenateBill();

end

n = size(W2,2);
%shuffle input adjacency matrix
idx_rand = randperm(n);% shuffle the nodes
[~, idx_reverse] = sort(idx_rand);
W2 = W2(idx_rand,idx_rand);
W3 = W3(idx_rand,idx_rand);

for c3 = c3_array % weight of triangles


    %estimate embedding using linear spectral clustering
    [x_est_linear] = LinearHypergraphEmbedding(W2, W3, c2, c3, "false");
    [x_est_periodic] = PeriodicHypergraphEmbedding(W2, W3, c2, c3, "false");

    %normalize the estimated embedding to the same range
    %x_est_linear = x_est_linear*norm(x_est_periodic,2)/norm(x_est_linear,2);        

    %reverse to the input order
    x_est_linear = x_est_linear(idx_reverse);
    x_est_periodic = x_est_periodic(idx_reverse);
    W2 = W2(idx_reverse, idx_reverse);
    W3 = W3(idx_reverse, idx_reverse);

    %reorder nodes according to the embedding
    [~,idx_linear] = sort(x_est_linear);
    W2_reorder_linear = W2(idx_linear,idx_linear);
    W3_reorder_linear = W3(idx_linear,idx_linear);
    
    %reorder nodes according to the embedding
    [~,idx_periodic] = sort(x_est_periodic);
    W2_reorder_periodic = W2(idx_periodic,idx_periodic);
    W3_reorder_periodic = W3(idx_periodic,idx_periodic);
    
    %compare likelihood
    lnP_linear = []; lnP0_linear = [];
    lnP_periodic = []; lnP0_periodic = [];

    for test_gamma = gamma_array
        [lnP_linear(end+1), lnP0_linear(end+1)] = CalculateModelLikelihood(x_est_linear, W2, T3, c2, c3, test_gamma, "linear");
        [lnP_periodic(end+1), lnP0_periodic(end+1)] = CalculateModelLikelihood(x_est_periodic, W2, T3, c2, c3, test_gamma, "periodic");

    end
    
    %maximum likelihood of gamma
        [~, max_linear_idx] = max(lnP_linear);
        gamma_max_linear = gamma_array(max_linear_idx); 
        [~, max_periodic_idx] = max(lnP_periodic);
        gamma_max_periodic = gamma_array(max_periodic_idx); 

        
        % plot likelihood
        plt = plot(gamma_array, lnP_linear, 'b', 'LineWidth',1.5);
        hold on;
        plot(gamma_array, lnP_periodic, '-r', 'LineWidth',1.5);
        plt = plot(gamma_array, lnP_linear, 'b', 'LineWidth',1.5);
        %plot(gamma_array, lnP0_linear, '--b', 'LineWidth',1.5);
        %plot(gamma_array, lnP0_periodic, '--r', 'LineWidth',1.5);
        plot(gamma_max_linear, lnP_linear(max_linear_idx), 'ok', 'MarkerSize',10, 'LineWidth',2);
        plot(gamma_max_periodic, lnP_periodic(max_periodic_idx), 'or', 'MarkerSize',10, 'LineWidth',2);
        legend({'Linear','Periodic','MLE'},'FontSize', 20,'Location','southeast');
        xlabel('\gamma','FontSize', 13);
        ylabel('Log-likelihood','FontSize', 13);
        set(gca,'fontsize',30);
        set(gca,'XLim',[0 max(gamma_array)])
        plt.LineWidth = 2;
        ax = gca;
        exportgraphics(ax,strcat('plots/model_comparison_', data_type,'.eps'),'Resolution',300) 
        hold off;

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Generate plots
    % 2-d scatter plot for edges
    figure
    s = scatter(x_est_linear(E2(:,1)),x_est_linear(E2(:,2)),'MarkerFaceColor','black','MarkerEdgeColor','none');
    alpha(s,0.3) % transparent color
    set(gca,'fontsize',30);
    ax = gca;
    exportgraphics(ax,strcat('plots/periodic_',data_type,'_edges.eps'),'Resolution',300) 

    % 3-d scatter plot for triangle
    s = scatter3(E3(:,1),E3(:,2),E3(:,3),10,'MarkerFaceColor','black','MarkerEdgeColor','none');
    xlim([0 n])
    xlim([0 n])
    zlim([0 n])
    alpha(s,0.2) % transparent color
    set(gca,'fontsize',30);
    ax = gca;
    exportgraphics(ax,strcat('plots/',data_type,'_triangles.eps'),'Resolution',300) 

    % plot reordered W2
    imagesc(W2_reorder_linear,[0,1]); %plot color map of original matrix
    colormap(flipud(gray(2)));
    set(gca,'FontSize',30) ;
    ax = gca;% Requires R2020a or later
    exportgraphics(ax,strcat('plots/',data_type,'_W2_reorder_linear_c3=',num2str(round(c3,1)),'.eps'),'Resolution',300) 

    % plot reordered W2
    imagesc(W2_reorder_periodic,[0,1]); %plot color map of original matrix
    colormap(flipud(gray(2)));
    set(gca,'FontSize',30) ;
    ax = gca;% Requires R2020a or later
    exportgraphics(ax,strcat('plots/',data_type,'_W2_reorder_periodic_c3=',num2str(round(c3,1)),'.eps'),'Resolution',300) 



    % plot reordered W3
    imagesc(W3_reorder_linear); % plot color map of original matrix
    colormap(flipud(gray(256)));colorbar
    set(gca,'FontSize',30) ;
    set(gca,'ColorScale','log')
    ax = gca;% Requires R2020a or later
    exportgraphics(ax,strcat('plots/',data_type,'_W3_reorder_linear_c3=',num2str(round(c3,1)),'.eps'),'Resolution',300) 

    % plot reordered W3
    imagesc(W3_reorder_periodic); % plot color map of original matrix
    colormap(flipud(gray(256)));colorbar
    set(gca,'FontSize',30) ;
    set(gca,'ColorScale','log')
    ax = gca;% Requires R2020a or later
    exportgraphics(ax,strcat('plots/',data_type,'_W3_reorder_periodic_c3=',num2str(round(c3,1)),'.eps'),'Resolution',300) 

    

    cluster_input = table2array(label);
    cluster_est_linear = kmeans(x_est_linear, K);
    cluster_est_periodic = kmeans(x_est_periodic, K);
    rand_linear(end+1) = CalculateRandIndex(cluster_input, cluster_est_linear);
    rand_periodic(end+1) = CalculateRandIndex(cluster_input, cluster_est_periodic);



    

  %{
  % plot eigenvalues of L
    [V, lambda] = eigs(L,n,'smallestabs');
    eigenvalues = diag(lambda);
    scatter(1:40, eigenvalues(1:40), 100, 'MarkerFaceColor','black');
    xlabel('i','FontSize', 13);
    ylabel('Eigenvalue \lambda_i','FontSize', 13);
    set(gca,'fontsize',30);
    ax = gca;
    exportgraphics(ax,strcat('plots/periodic_',data_type,'_eigenvalues_c3=',num2str(round(c3,1)),'.eps'),'Resolution',300) 


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
    axis([0 n -1 1])
    exportgraphics(t,strcat('plots/',data_type,'_eigenvectors_c3=',num2str(round(c3,1)),'.eps'),'Resolution',300) 

    figure

        
    % plot original W2 and W3
    imagesc(W2,[0,1]); %plot color map of original matrix
    colormap(flipud(gray(2)));
    set(gca,'FontSize',30) ;
    ax = gca;% Requires R2020a or later
    exportgraphics(ax,strcat('plots/periodic_',data_type,'_W2_input.eps'),'Resolution',300) 


    imagesc(W3); % plot color map of original matrix
    colormap(flipud(gray(256)));colorbar
    set(gca,'FontSize',30) ;
    set(gca,'ColorScale','log')
    ax = gca;% Requires R2020a or later
    exportgraphics(ax,strcat('plots/periodic_',data_type,'_W3_input_.eps'),'Resolution',300) 

    % plot L2 eigenvalues
    % eigenvalues of L2 and L3
    d2 = sum(W2, 2); D2 = diag(d2);
    d3 = sum(W3, 2); D3 = diag(d3);
    d2_inv = 1./d2; d2_inv(isinf(d2_inv)|isnan(d2_inv)) = 0;
    d3_inv = 1./d3; d3_inv(isinf(d3_inv)|isnan(d3_inv)) = 0;
    L2 = D2 - W2;
    L3 = D3 - W3;
    %normalize Laplacian
    L2 = diag(d2_inv)*L2;
    L3 = diag(d3_inv)*L3;
    [V2, lambda2] = eigs(L2,n,'smallestabs'); % all eigenvectors ranging from smallest eigenvalue to largest eigenvalue
    [V3, lambda3] = eigs(L3,n,'smallestabs'); % all eigenvectors ranging from smallest eigenvalue to largest eigenvalue


    eigenvalues2 = diag(lambda2);
    scatter(1:40, eigenvalues2(1:40), 100, 'MarkerFaceColor','black');
    xlabel('i','FontSize', 13);
    ylabel('Eigenvalue \lambda_i','FontSize', 13);
    set(gca,'fontsize',30);
    ax = gca;
    exportgraphics(ax,strcat('plots/periodic_',data_type,'_L2_eigenvalues.eps'),'Resolution',300) 

    % plot L3 eigenvalues
    eigenvalues3 = diag(lambda3);
    scatter(1:40, eigenvalues3(1:40), 100, 'MarkerFaceColor','black');
    xlabel('i','FontSize', 13);
    ylabel('Eigenvalue \lambda_i','FontSize', 13);
    set(gca,'fontsize',30);
    ax = gca;
    exportgraphics(ax,strcat('plots/periodic_',data_type,'_L3_eigenvalues.eps'),'Resolution',300) 

    % plot  eigenvectors of L2
    t = tiledlayout(4,1);
    ax1 = nexttile;
    plot(V2(:,1), 'Color', 'black');
    ylabel('v_0')
    ax2 = nexttile;
    plot(V2(:,2), 'Color', 'black');
    ylabel('v_1')
    ax3 = nexttile;
    plot(V2(:,3), 'Color', 'black');
    ylabel('v_2')
    ax4 = nexttile;
    plot(V2(:,4), 'Color', 'black');
    ylabel('v_3')
    % Link the axes
    linkaxes([ax1,ax2,ax3,ax4],'x');
    linkaxes([ax1,ax2,ax3,ax4],'y');
    % Move plots closer together
    xticklabels(ax1,{})
    t.TileSpacing = 'compact';
    axis([0 n -1 1])
    exportgraphics(t,strcat('plots/periodic_',data_type,'_L2_eigenvectors.eps'),'Resolution',300) 


    % plot  eigenvectors of L2
    t = tiledlayout(4,1);
    ax1 = nexttile;
    plot(V3(:,1), 'Color', 'black');
    ylabel('v_0')
    ax2 = nexttile;
    plot(V3(:,2), 'Color', 'black');
    ylabel('v_1')
    ax3 = nexttile;
    plot(V3(:,3), 'Color', 'black');
    ylabel('v_2')
    ax4 = nexttile;
    plot(V3(:,4), 'Color', 'black');
    ylabel('v_3')
    % Link the axes
    linkaxes([ax1,ax2,ax3,ax4],'x');
    linkaxes([ax1,ax2,ax3,ax4],'y');
    % Move plots closer together
    xticklabels(ax1,{})
    t.TileSpacing = 'compact';
    axis([0 n -1 1])
    exportgraphics(t,strcat('plots/periodic_',data_type,'_L3_eigenvectors.eps'),'Resolution',300) 


%}
    
end

plt = plot(c3_array, rand_linear, 'b', 'LineWidth',1.5);
hold on;
plot(c3_array, rand_periodic, '-r', 'LineWidth',1.5);
legend({'Linear','Periodic'},'FontSize', 20,'Location','southeast');
xlabel('c_3','FontSize', 13);
ylabel('Rand Index','FontSize', 13);
set(gca,'fontsize',30);
set(gca,'YLim',[0.5 1.1])
plt.LineWidth = 2;
ax = gca;
exportgraphics(ax,strcat('plots/rand_linear_', data_type,'.eps'),'Resolution',300) 
hold off;

toc
load handel
sound(y,Fs)