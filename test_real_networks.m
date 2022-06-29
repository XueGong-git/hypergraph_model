clear all 
clc
close all

% folder that contains algorithm scripts
addpath 'functions';
set(groot,'defaultFigureVisible','off') % 'on' to turn back on.  
tic


c2 = 0; % weight of simple edge
c3_array = 0.3;
gamma_array = 0:0.2:2; % gamma for likelihood plot
data_type = "highschool";
K = 9;
rand_linear = [];
rand_periodic = [];
max_lnP_linear = [];
max_lnP_periodic = [];
% likelihood
lnP_linear = zeros(length(c3_array),length(gamma_array));
lnP_periodic = zeros(length(c3_array),length(gamma_array));

filename = strcat(data_type,'.mat');

if isfile(filename)
    if data_type == "highschool"
        load(filename,'W2', 'W3', 'T3', 'E2', 'E3', 'label')
    elseif data_type == "senate_bill"
        load(filename,'W2', 'W3', 'T3', 'E2', 'E3')
    end

else
    if data_type == "highschool"
        [W2, W3, T3, E2, E3] = LoadHighSchool();
        label = readtable('raw_data/highschool_label', 'ReadVariableNames', false);
        save('highschool','W2', 'W3', 'T3', 'E2', 'E3','label');
    elseif data_type == "senate_bill"
        [W2, W3, T3, E2, E3] = LoadSenateBill();
        save('senate_bill','W2', 'W3', 'T3', 'E2', 'E3');
    end
end


n_nodes = size(W2,2);
n_edge = sum(W2,'all')/2;
n_triangle = sum(T3, 'all')/6;
edge_density = 2*n_edge/(n_nodes*(n_nodes-1));
triangle_density = 6*n_triangle/(n_nodes*(n_nodes-1)*(n_nodes-2));

n = size(W2,2);
%shuffle input adjacency matrix
idx_rand_in = randperm(n);% shuffle the nodes
[~, idx_reverse] = sort(idx_rand_in);
W2_in = W2(idx_rand_in,idx_rand_in);
W3_in = W3(idx_rand_in,idx_rand_in);
T3_in = T3(idx_rand_in, idx_rand_in, idx_rand_in);

for ii = 1:length(c3_array)
    c3 = c3_array(ii);

    %check degree distribution
    degree_eff = sum(c2*W2_in + c3*W3_in, 2);
    norm_eta = c2*n_edge + c3*n_triangle;
    f=figure;
    histogram(degree_eff);
    saveas(f,strcat('plots/highschool_degree_eff_c3=',num2str(round(c3,2)),'.eps'));

    %trim  nodes with highest degrees
    keepIDs = (degree_eff< max(degree_eff)*1.5) & (degree_eff> min(degree_eff)*0.4);
    W2 = W2_in(keepIDs, keepIDs);
    T3 = T3_in(keepIDs, keepIDs, keepIDs);
    W3 = sum(T3, 3);
    idx_rand = idx_rand_in(keepIDs);
    [~, idx_reverse] = sort(idx_rand);

    %check degree distribution again
    degree_eff = sum(c2*W2 + c3*W3, 2);
    f=figure;
    histogram(degree_eff);
    saveas(f,strcat('plots/highschool_degree_eff_trim_c3=',num2str(round(c3,2)),'_c2_',num2str(round(c2,2)),'.eps'));

    
    %estimate embedding using linear spectral clustering
    [x_est_linear] = LinearHypergraphEmbedding(W2, W3, c2, c3, "false", 3);
    [x_est_periodic] = PeriodicHypergraphEmbedding(W2, W3, c2, c3, "false");
                
    
    %calculate eta
    [~, eta_linear_est] = CalculateModelLikelihood(x_est_linear, W2, T3, c2, c3, 2, "linear");
    [~, eta_periodic_est] = CalculateModelLikelihood(x_est_periodic, W2, T3, c2, c3, 2, "linear");

    %scale estimation
    x_est_linear = x_est_linear*sqrt(norm_eta/eta_linear_est);  
    x_est_periodic = x_est_periodic*sqrt(norm_eta/eta_periodic_est);  

    [~, eta_linear_scaled] = CalculateModelLikelihood(x_est_linear, W2, T3, c2, c3, 2, "linear");
    [~, eta_periodic_scaled] = CalculateModelLikelihood(x_est_linear, W2, T3, c2, c3, 2, "linear");


    %reorder nodes according to the embedding
    [~,idx_linear] = sort(x_est_linear);
    W2_reorder_linear = W2(idx_linear,idx_linear);
    W3_reorder_linear = W3(idx_linear,idx_linear);
    
    %reorder nodes according to the embedding
    [~,idx_periodic] = sort(x_est_periodic);
    W2_reorder_periodic = W2(idx_periodic,idx_periodic);
    W3_reorder_periodic = W3(idx_periodic,idx_periodic);

    for jj = 1:length(gamma_array)
        gamma = gamma_array(jj);
        [lnP_linear(ii,jj), ~] = CalculateModelLikelihood(x_est_linear, W2, T3, c2, c3, gamma, "linear");
        [lnP_periodic(ii,jj), ~] = CalculateModelLikelihood(x_est_periodic, W2, T3, c2, c3, gamma, "periodic");

    end
    
    %maximum likelihood of gamma
    [max_lnP_linear(end+1), max_linear_idx] = max(lnP_linear(ii,:));
    gamma_max_linear = gamma_array(max_linear_idx); 
    [max_lnP_periodic(end+1), max_periodic_idx] = max(lnP_periodic(ii,:));
    gamma_max_periodic = gamma_array(max_periodic_idx); 

    % plot likelihood
    plt = plot(gamma_array, lnP_linear(ii,:), 'b', 'LineWidth',1.5);
    hold on;
    plot(gamma_array, lnP_periodic(ii,:), '-r', 'LineWidth',1.5);
    plot(gamma_max_linear, lnP_linear(ii, max_linear_idx), 'ok', 'MarkerSize',10, 'LineWidth',2);
    plot(gamma_max_periodic, lnP_periodic(ii, max_periodic_idx), 'or', 'MarkerSize',10, 'LineWidth',2);
    legend({'Linear','Periodic','MLE'},'FontSize', 20,'Location','best');
    xlabel('\gamma','FontSize', 13);
    ylabel('Log-likelihood','FontSize', 13);
    set(gca,'fontsize',30);
    set(gca,'XLim',[0 max(gamma_array)])
    plt.LineWidth = 2;
    ax = gca;
    exportgraphics(ax,strcat('plots/model_comparison_', data_type,'_c2_',num2str(round(c2,2)),'_c3_',num2str(round(c3,2)),'.eps'),'Resolution',300) 
    hold off;

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Generate plots
    %{
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
%}
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

    %{
    
    cluster_input = table2array(label);
    cluster_input = cluster_input(keepIDs);
    cluster_est_linear = kmeans(x_est_linear, K);
    cluster_est_periodic = kmeans(x_est_periodic, K);
    rand_linear(end+1) = CalculateRandIndex(cluster_input, cluster_est_linear);
    rand_periodic(end+1) = CalculateRandIndex(cluster_input, cluster_est_periodic);

    %}
    

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

%{
%plot rand index
plt = plot(c3_array, rand_linear, 'b', 'LineWidth',1.5);
hold on;
plot(c3_array, rand_periodic, '-r', 'LineWidth',1.5);
legend({'Linear','Periodic'},'FontSize', 20,'Location','southeast');
xlabel('c_3','FontSize', 13);
ylabel('Rand Index','FontSize', 13);
set(gca,'fontsize',30);
set(gca,'YLim',[0.5 1.1])
ax = gca;
exportgraphics(ax,strcat('plots/rand_linear_', data_type,'.eps'),'Resolution',300) 
hold off;
%}

plt = plot(c3_array, max_lnP_linear, 'b', 'LineWidth',1.5);
hold on;
plot(c3_array, max_lnP_periodic, '-r', 'LineWidth',1.5);
legend({'Linear','Periodic'},'FontSize', 20,'Location','southeast');
xlabel('c_3','FontSize', 13);
ylabel('Max log likelihood','FontSize', 13);
set(gca,'fontsize',30);
ax = gca;
exportgraphics(ax,strcat('plots/max_lnP', data_type,'_c2_',num2str(round(c2,2)),'.eps'),'Resolution',300) 
hold off;

% plot heatmap of linear log-likelihood
h = heatmap(round(lnP_linear, 2),'Colormap',parula);
h.Title = strcat('LnP Periodic' );
h.XLabel = '\gamma';
h.YLabel = 'c3_{inv}';
h.XData = round(gamma_array,2);
h.YData = round(c3_array,2);
ax = gca;% Requires R2020a or later
caxis([-2000000, 0]);
exportgraphics(ax,strcat('plots/lnP_linear_', data_type,'_c2_',num2str(round(c2,2)),'.eps'),'Resolution',300) 

% plot heatmap of periodic log-likelihood
h = heatmap(round(lnP_periodic, 2),'Colormap',parula);
h.Title = strcat('LnP Periodic' );
h.XLabel = '\gamma';
h.YLabel = 'c3_{inv}';
h.XData = round(gamma_array,2);
h.YData = round(c3_array,2);
ax = gca;% Requires R2020a or later
caxis([-2000000, 0]);
exportgraphics(ax,strcat('plots/lnP_periodic_', data_type,'_c2_',num2str(round(c2,2)),'.eps'),'Resolution',300) 


toc
%load handel
%sound(y,Fs)