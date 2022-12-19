clear all 
clc
close all

% folder that contains algorithm scripts
addpath 'functions';
set(groot,'defaultFigureVisible','off') % 'on' to turn figures back on.  
tic


c2 = 1; % weight of simple edge
c3_array = 0:0.1:2;
gamma_array = 0:100:2000; % gamma for likelihood plot
data_type = "highschool";
%data_type = "primaryschool";
%data_type = "senate-committees";
rand_linear = [];
rand_periodic = [];
max_lnP_linear = [];
max_lnP_periodic = [];
% likelihood
lnP_linear = zeros(length(c3_array),length(gamma_array));
lnP_periodic = zeros(length(c3_array),length(gamma_array));

filename = strcat(data_type,'.mat');

if data_type == "highschool"
    trim = "false";
    K = 9;
    n_eig = 3;
    if isfile(filename)
        load(filename,'W2', 'W3', 'T3', 'E2', 'E3', 'label')
    else
        [W2, W3, T3, E2, E3] = LoadHighSchool();
        label = readtable('raw_data/node-labels-contact-high-school.txt', 'ReadVariableNames', false);
        save('highschool','W2', 'W3', 'T3', 'E2', 'E3','label');
    end
elseif data_type == "primaryschool"
    trim = "false";
    K = 11;
    n_eig = 4;
    if isfile(filename)
        load(filename,'W2', 'W3', 'T3', 'E2', 'E3', 'label')
    else
        [W2, W3, T3, E2, E3] = LoadPrimarySchool();
        label = readtable('raw_data/primaryschool_label', 'ReadVariableNames', false);
        save('primaryschool','W2', 'W3', 'T3', 'E2', 'E3','label');
    end
elseif data_type == "senate_bill"
    trim = "true";
    K = 2;
    n_eig = 1;
    if isfile(filename)
        load(filename,'W2', 'W3', 'T3', 'E2', 'E3')
    else [W2, W3, T3, E2, E3] = LoadSenateBill();
        label = readtable('raw_data/senate-bills/node-labels-senate-bills.txt', 'ReadVariableNames', false);
        save('senate_bill','W2', 'W3', 'T3', 'E2', 'E3','label');
    end
elseif data_type == "senate-committees"
    K = 2;
    n_eig = 1;
    if isfile(filename)
        load(filename,'W2', 'W3', 'T3', 'E2', 'E3')
    else [W2, W3, T3, E2, E3] = LoadSenateCommittees();
        label = readtable('raw_data/senate-committees/node-labels-senate-committees.txt', 'ReadVariableNames', false);
        save('house_bill','W2', 'W3', 'T3', 'E2', 'E3','label');
    end
end


n_nodes = size(W2,2);
n_edge = sum(W2,'all')/2;
n_triangle = sum(T3, 'all')/6;
edge_density = 2*n_edge/(n_nodes*(n_nodes-1));
triangle_density = 6*n_triangle/(n_nodes*(n_nodes-1)*(n_nodes-2));

n = size(W2,2);
cluster_input = table2array(label);
%shuffle input adjacency matrix
idx_rand_in = randperm(n);% shuffle the nodes
%[~, idx_reverse] = sort(idx_rand_in);
W2 = W2(idx_rand_in,idx_rand_in);
T3 = T3(idx_rand_in, idx_rand_in, idx_rand_in);
W3 = sum(T3, 3);
cluster_input = cluster_input(idx_rand_in);

for ii = 1:length(c3_array)
    c3 = c3_array(ii);

    %check degree distribution
    degree_eff = sum(c2*W2 + c3*W3, 2);
    norm_eta = c2*n_edge + c3*n_triangle;
    f=figure;
    histogram(degree_eff);
    saveas(f,strcat('plots/highschool_degree_eff_c3=',num2str(round(c3,2)),'.eps'));

    %trim  nodes with highest degrees
    if trim == "true"
        keepIDs = (degree_eff > quantile(degree_eff,0.02));
        W2 = W2(keepIDs, keepIDs);
        T3 = T3(keepIDs, keepIDs, keepIDs);
        W3 = sum(T3, 3);
        cluster_input = cluster_input(keepIDs);
        %idx_rand = idx_rand_in(keepIDs);
        %[~, idx_reverse] = sort(idx_rand);
        %label = label(keepIDs,1);

        %check degree distribution again
        degree_eff = sum(c2*W2 + c3*W3, 2);
        f=figure;
        histogram(degree_eff);
        saveas(f,strcat('plots/highschool_degree_eff_trim_c3=',num2str(round(c3,2)),'_c2_',num2str(round(c2,2)),'.eps'));
    end
    %estimate embedding using linear spectral clustering
    [x_est_linear] = LinearHypergraphEmbedding(W2, W3, c2, c3, "false", n_eig); %number of eigen vectors
    [x_est_periodic] = PeriodicHypergraphEmbedding(W2, W3, c2, c3, "false");
                
    
    %calculate eta
    [~, eta_linear_est] = CalculateModelLikelihood(x_est_linear, W2, T3, c2, c3, 2, "linear");
    [~, eta_periodic_est] = CalculateModelLikelihood(x_est_periodic, W2, T3, c2, c3, 2, "periodic");

    %scale estimation
    x_est_linear = x_est_linear*sqrt(eta_periodic_est/eta_linear_est);  
    %x_est_periodic = x_est_periodic*sqrt(norm_eta/eta_periodic_est); 

    [~, eta_linear_scaled] = CalculateModelLikelihood(x_est_linear, W2, T3, c2, c3, 2, "linear");
    %[~, eta_periodic_scaled] = CalculateModelLikelihood(x_est_periodic, W2, T3, c2, c3, 2, "periodic");


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
    plt = plot(gamma_array, lnP_linear(ii,:), '-*b', 'LineWidth',1.5);
    hold on;
    plot(gamma_array, lnP_periodic(ii,:), '-or', 'LineWidth',1.5);
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
    
    
    %cluster_input = cluster_input(keepIDs);
    cluster_est_linear = kmeans(x_est_linear, K);
    cluster_est_periodic = kmeans(x_est_periodic, K);
    rand_linear(end+1) = CalculateRandIndex(cluster_input, cluster_est_linear, 'adjusted');
    rand_periodic(end+1) = CalculateRandIndex(cluster_input, cluster_est_periodic, 'adjusted');

    
end


save(strcat('test_',data_type,'.mat'), 'rand_linear', ...
    'rand_periodic', 'x_est_linear', 'x_est_periodic', 'max_lnP_linear', 'max_lnP_periodic', 'edge_density', 'triangle_density')

%load(strcat('test_',data_type,'.mat'), 'rand_linear', ...
 %   'rand_periodic', 'x_est_linear', 'x_est_periodic', 'max_lnP_linear', 'max_lnP_periodic', 'edge_density', 'triangle_density')



%plot rand index
plt = plot(c3_array, rand_linear, '-*b', 'LineWidth',1.5);
hold on;
plot(c3_array, rand_periodic, '--or', 'LineWidth',1.5);
legend({'Linear','Periodic'},'FontSize', 20,'Location','southeast');
xlabel('c_3*','FontSize', 13);
ylabel('ARI','FontSize', 13);
set(gca,'fontsize',30);
set(gca,'YLim')
ax = gca;
exportgraphics(ax,strcat('plots/rand_linear_', data_type,'.eps'),'Resolution',300) 
hold off;


plt = plot(c3_array, max_lnP_linear, '-*b', 'LineWidth',1.5);
hold on;
plot(c3_array, max_lnP_periodic, '--or', 'LineWidth',1.5);
[global_max_lnP_lin, max_lin_idx] = max(max_lnP_linear);
[global_max_lnP_per, max_per_idx] = max(max_lnP_periodic);
plot(c3_array(max_lin_idx), global_max_lnP_lin, 'pr', 'MarkerSize',15, 'MarkerFaceColor','red');
plot(c3_array(max_per_idx), global_max_lnP_per, 'pr', 'MarkerSize',15, 'MarkerFaceColor','red');
legend({'Linear','Periodic'},'FontSize', 20,'Location','southeast');
xlabel('c_3*','FontSize', 13);
ylabel('Maximum LnP','FontSize', 13);
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
 


save(strcat('results/',data_type,'_clustering.mat'), 'rand_linear', ...
    'rand_periodic',  'max_lnP_linear', 'max_lnP_periodic', 'edge_density', 'triangle_density', ...
      'n','K','c2', 'c3_array','gamma_array', 'trim', 'n_eig')


toc
load chirp.mat
sound(y,Fs)