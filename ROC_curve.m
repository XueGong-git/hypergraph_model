clear all 
clc
close all

% folder that contains algorithm scripts
addpath 'functions';
set(groot,'defaultFigureVisible','off') % 'on' to turn back on.  
tic

K = 9; %number of clusters


c2 = 1; % weight of simple edge
c3_array = 0:0.1:1;
gamma_array = 0.5; % gamma for likelihood plot
data_name = "synthetic_lin";
train_ratio = 0.8;
train_file = strcat('processed_data/', data_name,'_train.mat');
test_file = strcat('processed_data/', data_name,'_test.mat');
m = 20; % number of nodes per cluster
gamma_inp = 18; % decay parameter for synthetic networks
c2_inp = 1;
c3_inp = 1/3;
%%%%%%%%% Parameters only applicable to synthetic graphs %%%%%%%%%
if data_name == "synthetic_lin" || data_name == "synthetic_per"
    n = K*m; % total number of nodes
end
%%%%%%%%%%%%%%%%%%%%%% 

rand_linear = [];
rand_periodic = [];
max_lnP_linear = [];
max_lnP_periodic = [];
AUC_edge_lin = [];
AUC_edge_per = [];
AUC_triangle_lin = [];
AUC_triangle_per = [];
% likelihood
lnP_linear = zeros(length(c3_array),length(gamma_array));
lnP_periodic = zeros(length(c3_array),length(gamma_array));



if ismember (data_name, ["contact-high-school", "contact-primary-school", "coauth-DBLP", "email-Eu", "email-Enron" ] )
    if  ~isfile(train_file) || ~isfile(test_file)
        [train_list, test_list, n] = SplitData(data_name,train_ratio);
        [W2, W3, T3, E2, E3] = LoadSimplex(train_list, n);
        save(strcat('processed_data/', data_name, '_train'),'W2', 'W3', 'T3', 'E2', 'E3');
        [W2, W3, T3, E2, E3] = LoadSimplex(train_list, n);
        save(strcat('processed_data/', data_name, '_test'),'W2', 'W3', 'T3', 'E2', 'E3');
    end
    train = load(train_file,'W2', 'W3', 'T3', 'E2', 'E3');
    test  = load(test_file,'W2', 'W3', 'T3', 'E2', 'E3');
    n = size(test.W2,1);
elseif data_name == "senate_bill"
    if isfile(filename)
        load(filename,'W2', 'W3', 'T3', 'E2', 'E3', 'label')
    else
        [W2, W3, T3, E2, E3] = LoadSenateBill();
        save('processed_data/senate_bill','W2', 'W3', 'T3', 'E2', 'E3');
    end
elseif data_name == "synthetic_lin"
    %if isfile(filename)
     %   load(filename,'W2', 'W3', 'T3', 'label')
    %else
        a = 0.01; % noise parameter
        x = sort(repmat(linspace(0,2, K),1,m)+(2*a*rand(1,n)-a)); % input node embedding
        [E2, E3, W2, W3, T3] = GenerateLinearHypergraph(x, gamma_inp, c2_inp, c3_inp);   
        label = repmat(linspace(1,K,K),1,m);
        split_index = randperm(n, round(train_ratio*n));
        save('processed_data/synthetic_lin','W2', 'W3', 'T3', 'label');
    %end
elseif data_name == "synthetic_per"
    a = 0.01*pi; % noise parameter
    x = sort(repmat(linspace(-pi,pi,K),1,m)+(2*a*rand(1,n)-a)); % angles from -pi to pi
    [W2, W3, T3] = GeneratePeriodicHypergraph(x, gamma_inp, c2_inp, c3_inp);
    label = repmat(linspace(1,K,K),1,m);

end

%extract largest connected component
%[x, W2, W3, T3] = MaxConnectedSubgraph(x, c2, c3, W2, W3, T3);    


%plot the input adjacency matrix
imagesc(train.W2,[0,1]); %plot color map of original matrix
colormap(flipud(gray(2)));
set(gca,'FontSize',30) ;
ax = gca;% Requires R2020a or later
exportgraphics(ax,strcat('plots/',data_name,'_gamma=', num2str(gamma_inp),'_input_W2_c2_',num2str(round(c2_inp,2)),'.eps'),'Resolution',300) 


%plot the input triangle adjacency matrix
imagesc(train.W3); %plot color map of original matrix
colormap(flipud(gray(256)));
set(gca,'FontSize',30) ;
ax = gca;% Requires R2020a or later
exportgraphics(ax,strcat('plots/',data_name,'_gamma=', num2str(gamma_inp),'_input_W3_c3_',num2str(round(c3_inp,2)),'.eps'),'Resolution',300) 


n_nodes = size(train.W2,2);
n_edge = sum(train.W2,'all')/2;
n_triangle = sum(train.T3, 'all')/6;
edge_density = 2*n_edge/(n_nodes*(n_nodes-1));
triangle_density = 6*n_triangle/(n_nodes*(n_nodes-1)*(n_nodes-2));

%shuffle input adjacency matrix
idx_rand_in = randperm(n_nodes);% shuffle the nodes
[~, idx_reverse] = sort(idx_rand_in);
W2_in = train.W2;
W3_in = train.W3;
T3_in = train.T3;

for ii = 1:length(c3_array)
    c3 = c3_array(ii);
    norm_eta = c2*n_edge + c3*n_triangle;

    
    %{
    %check degree distribution
    degree_eff = sum(c2*W2_in + c3*W3_in, 2);
    f=figure;
    histogram(degree_eff);
    saveas(f,strcat('plots/highschool_degree_eff_c3=',num2str(round(c3,2)),'.eps'));

    if data_name == "highschool"

        %trim  nodes with highest degrees
        keepIDs = (degree_eff< max(degree_eff)*0.5) & (degree_eff> min(degree_eff)*4);
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

    end
    %}
    
    %estimate embedding using linear spectral clustering
    [x_est_linear] = LinearHypergraphEmbedding(W2_in, W3_in, c2, c3, "false");
    [x_est_periodic] = PeriodicHypergraphEmbedding(W2_in, W3_in, c2, c3, "false");
                
    
    %calculate eta
    [~, eta_linear_est, ~] = CalculateModelLikelihood(x_est_linear, W2_in, T3_in, c2, c3, 2, "linear");
    [~, eta_periodic_est, ~] = CalculateModelLikelihood(x_est_periodic, W2_in, T3_in, c2, c3, 2, "linear");

    %scale estimation
    x_est_linear = x_est_linear*sqrt(norm_eta/eta_linear_est);  
    x_est_periodic = x_est_periodic*sqrt(norm_eta/eta_periodic_est);  

    [~, eta_linear_scaled, ~] = CalculateModelLikelihood(x_est_linear, W2_in, T3_in, c2, c3, 2, "linear");
    [~, eta_periodic_scaled, ~] = CalculateModelLikelihood(x_est_linear, W2_in, T3_in, c2, c3, 2, "linear");


    %reorder nodes according to the embedding
    [~,idx_linear] = sort(x_est_linear);
    W2_reorder_linear = W2_in(idx_linear,idx_linear);
    W3_reorder_linear = W3_in(idx_linear,idx_linear);
    
    %reorder nodes according to the embedding
    [~,idx_periodic] = sort(x_est_periodic);
    W2_reorder_periodic = W2_in(idx_periodic,idx_periodic);
    W3_reorder_periodic = W3_in(idx_periodic,idx_periodic);

    for jj = 1:length(gamma_array)
        gamma = gamma_array(jj);
        [lnP_linear(ii,jj), ~] = CalculateModelLikelihood(x_est_linear, W2_in, T3_in, c2, c3, gamma, "linear");
        [lnP_periodic(ii,jj), ~] = CalculateModelLikelihood(x_est_periodic, W2_in, T3_in, c2, c3, gamma, "periodic");
    end
    
    %maximum likelihood of gamma
    [max_lnP_linear(end+1), max_linear_idx] = max(lnP_linear(ii,:));
    gamma_max_linear = gamma_array(max_linear_idx); 
    
    [~, ~, P_edge_lin, P_triangles_lin] = CalculateModelLikelihood(x_est_linear, W2_in, T3_in, c2, c3, gamma_max_linear, "linear");
    
    %predict edge
    [X,Y,T,AUC_edge_lin(end+1)] = perfcurve(test.W2(:),P_edge_lin(:),1,'XCrit','reca','YCrit','prec');
    plot(X,Y)
    xlabel('Recall') 
    ylabel('Precision')
    title('ROC for edge prediction (linear)')
    set(gca,'fontsize',30);
    ax = gca;
    exportgraphics(ax,strcat('plots/ROC_edge_lin_', data_name,'_c2_',num2str(round(c2,2)),'_c3_',num2str(round(c3,2)),'.eps'),'Resolution',300) 

    %triangle
    [X,Y,T,AUC_triangle_lin(end+1)] = perfcurve(test.T3(:),P_triangles_lin(:),1,'XCrit','reca','YCrit','prec');
    plot(X,Y)
    xlabel('Recall') 
    ylabel('Precision')
    title('ROC for triangle prediction (linear)')
    set(gca,'fontsize',30);
    ax = gca;
    exportgraphics(ax,strcat('plots/ROC_triangle_lin_', data_name,'_c2_',num2str(round(c2,2)),'_c3_',num2str(round(c3,2)),'.eps'),'Resolution',300) 

    % periodic embedding
    [max_lnP_periodic(end+1), max_periodic_idx] = max(lnP_periodic(ii,:));
    gamma_max_periodic = gamma_array(max_periodic_idx); 
    [~, ~, P_edge_per, P_triangles_per] = CalculateModelLikelihood(x_est_periodic, W2_in, T3_in, c2, c3, gamma_max_periodic, "periodic");
    
    %predict edge
    [X,Y,T,AUC_edge_per(end+1)] = perfcurve(test.W2(:),P_edge_per(:),1,'XCrit','reca','YCrit','prec');
    plot(X,Y)
    xlabel('Recall') 
    ylabel('Precision')
    title('ROC for edge prediction (periodic)')
    set(gca,'fontsize',30);
    ax = gca;
    exportgraphics(ax,strcat('plots/ROC_edge_per_', data_name,'_c2_',num2str(round(c2,2)),'_c3_',num2str(round(c3,2)),'.eps'),'Resolution',300) 

    
    %triangle
    [X,Y,T,AUC_triangle_per(end+1)] = perfcurve(test.T3(:),P_triangles_per(:),1,'XCrit','reca','YCrit','prec');
    plot(X,Y)
    xlabel('Recall') 
    ylabel('Precision')
    title('ROC for triangle prediction (periodic)')
    set(gca,'fontsize',30);
    ax = gca;
    exportgraphics(ax,strcat('plots/ROC_triangle_per_', data_name,'_c2_',num2str(round(c2,2)),'_c3_',num2str(round(c3,2)),'.eps'),'Resolution',300) 

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
    exportgraphics(ax,strcat('plots/model_comparison_', data_name,'_c2_',num2str(round(c2,2)),'_c3_',num2str(round(c3,2)),'.eps'),'Resolution',300) 
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
    exportgraphics(ax,strcat('plots/periodic_',data_name,'_edges.eps'),'Resolution',300) 

    % 3-d scatter plot for triangle
    s = scatter3(E3(:,1),E3(:,2),E3(:,3),10,'MarkerFaceColor','black','MarkerEdgeColor','none');
    xlim([0 n])
    xlim([0 n])
    zlim([0 n])
    alpha(s,0.2) % transparent color
    set(gca,'fontsize',30);
    ax = gca;
    exportgraphics(ax,strcat('plots/',data_name,'_triangles.eps'),'Resolution',300) 
%}
    
    % plot reordered W2
    imagesc(W2_reorder_linear,[0,1]); %plot color map of original matrix
    colormap(flipud(gray(2)));
    set(gca,'FontSize',30) ;
    ax = gca;% Requires R2020a or later
    exportgraphics(ax,strcat('plots/',data_name,'_W2_reorder_linear_c3=',num2str(round(c3,1)),'.eps'),'Resolution',300) 

    % plot reordered W2
    imagesc(W2_reorder_periodic,[0,1]); %plot color map of original matrix
    colormap(flipud(gray(2)));
    set(gca,'FontSize',30) ;
    ax = gca;% Requires R2020a or later
    exportgraphics(ax,strcat('plots/',data_name,'_W2_reorder_periodic_c3=',num2str(round(c3,1)),'.eps'),'Resolution',300) 

    % plot reordered W3
    imagesc(W3_reorder_linear); % plot color map of original matrix
    colormap(flipud(gray(256)));colorbar
    set(gca,'FontSize',30) ;
    set(gca,'ColorScale','log')
    ax = gca;% Requires R2020a or later
    exportgraphics(ax,strcat('plots/',data_name,'_W3_reorder_linear_c3=',num2str(round(c3,1)),'.eps'),'Resolution',300) 

    % plot reordered W3
    imagesc(W3_reorder_periodic); % plot color map of original matrix
    colormap(flipud(gray(256)));colorbar
    set(gca,'FontSize',30) ;
    set(gca,'ColorScale','log')
    ax = gca;% Requires R2020a or later
    exportgraphics(ax,strcat('plots/',data_name,'_W3_reorder_periodic_c3=',num2str(round(c3,1)),'.eps'),'Resolution',300) 

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
    exportgraphics(ax,strcat('plots/periodic_',data_name,'_eigenvalues_c3=',num2str(round(c3,1)),'.eps'),'Resolution',300) 


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
    exportgraphics(t,strcat('plots/',data_name,'_eigenvectors_c3=',num2str(round(c3,1)),'.eps'),'Resolution',300) 

    figure

        
    % plot original W2 and W3
    imagesc(W2,[0,1]); %plot color map of original matrix
    colormap(flipud(gray(2)));
    set(gca,'FontSize',30) ;
    ax = gca;% Requires R2020a or later
    exportgraphics(ax,strcat('plots/periodic_',data_name,'_W2_input.eps'),'Resolution',300) 


    imagesc(W3); % plot color map of original matrix
    colormap(flipud(gray(256)));colorbar
    set(gca,'FontSize',30) ;
    set(gca,'ColorScale','log')
    ax = gca;% Requires R2020a or later
    exportgraphics(ax,strcat('plots/periodic_',data_name,'_W3_input_.eps'),'Resolution',300) 

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
    exportgraphics(ax,strcat('plots/periodic_',data_name,'_L2_eigenvalues.eps'),'Resolution',300) 

    % plot L3 eigenvalues
    eigenvalues3 = diag(lambda3);
    scatter(1:40, eigenvalues3(1:40), 100, 'MarkerFaceColor','black');
    xlabel('i','FontSize', 13);
    ylabel('Eigenvalue \lambda_i','FontSize', 13);
    set(gca,'fontsize',30);
    ax = gca;
    exportgraphics(ax,strcat('plots/periodic_',data_name,'_L3_eigenvalues.eps'),'Resolution',300) 

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
    exportgraphics(t,strcat('plots/periodic_',data_name,'_L2_eigenvectors.eps'),'Resolution',300) 


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
    exportgraphics(t,strcat('plots/periodic_',data_name,'_L3_eigenvectors.eps'),'Resolution',300) 


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
exportgraphics(ax,strcat('plots/rand_linear_', data_name,'.eps'),'Resolution',300) 
hold off;
%}
plt = plot(c3_array, max_lnP_linear, 'b', 'LineWidth',1.5);
hold on;
plot(c3_array, max_lnP_periodic, '-r', 'LineWidth',1.5);
legend({'Linear','Periodic'},'FontSize', 20,'Location','best');
xlabel('c_3','FontSize', 13);
ylabel('Max log likelihood','FontSize', 13);
set(gca,'fontsize',30);
ax = gca;
exportgraphics(ax,strcat('plots/max_lnP', data_name,'_c2_',num2str(round(c2,2)),'.eps'),'Resolution',300) 
hold off;

% AUC for edge prediction
plt = plot(c3_array, AUC_edge_lin, 'b', 'LineWidth',1.5);
hold on;
plot(c3_array, AUC_edge_per, '-r', 'LineWidth',1.5);
legend({'Linear','Periodic'},'FontSize', 20,'Location','best');
xlabel('c_3','FontSize', 13);
ylabel('AUC for edge prediction','FontSize', 13);
set(gca,'fontsize',30);
ax = gca;
exportgraphics(ax,strcat('plots/AUC_edge_', data_name,'_c2_',num2str(round(c2,2)),'_gamma_',num2str(round(gamma_inp,2)),'.eps'),'Resolution',300) 
hold off;

% AUC for triange prediction
plt = plot(c3_array, AUC_triangle_lin, 'b', 'LineWidth',1.5);
hold on;
plot(c3_array, AUC_triangle_per, '-r', 'LineWidth',1.5);
legend({'Linear','Periodic'},'FontSize', 20,'Location','best');
xlabel('c_3','FontSize', 13);
ylabel('AUC for triangle prediction','FontSize', 13);
set(gca,'fontsize',30);
ax = gca;
exportgraphics(ax,strcat('plots/AUC_triangle_', data_name,'_c2_',num2str(round(c2,2)),'_gamma_',num2str(round(gamma_inp,2)),'.eps'),'Resolution',300) 
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
exportgraphics(ax,strcat('plots/lnP_linear_', data_name,'_c2_',num2str(round(c2,2)),'.eps'),'Resolution',300) 

% plot heatmap of periodic log-likelihood
h = heatmap(round(lnP_periodic, 2),'Colormap',parula);
h.Title = strcat('LnP Periodic' );
h.XLabel = '\gamma';
h.YLabel = 'c3_{inv}';
h.XData = round(gamma_array,2);
h.YData = round(c3_array,2);
ax = gca;% Requires R2020a or later
caxis([-2000000, 0]);
exportgraphics(ax,strcat('plots/lnP_periodic_', data_name,'_c2_',num2str(round(c2,2)),'.eps'),'Resolution',300) 


[X,Y,T,AUC_triangle_arith] = perfcurve(test.T3(:),CalculateArithMean(test.W3),1,'XCrit','reca','YCrit','prec');
[X,Y,T,AUC_triangle_harm] = perfcurve(test.T3(:),CalculateHarmMean(test.W3),1,'XCrit','reca','YCrit','prec');
[X,Y,T,AUC_triangle_geo] = perfcurve(test.T3(:),CalculateGeoMean(test.W3),1,'XCrit','reca','YCrit','prec');
[X,Y,T,AUC_triangle_rand] = perfcurve(test.T3(:), rand(n*n*n,1),1,'XCrit','reca','YCrit','prec');

[X,Y,T,AUC_edge_rand] = perfcurve(test.W2(:), rand(n*n,1),1,'XCrit','reca','YCrit','prec');
[X,Y,T,AUC_edge_mean] = perfcurve(test.W2(:), test.W2(:),1,'XCrit','reca','YCrit','prec');


toc
%load handel
%sound(y,Fs)