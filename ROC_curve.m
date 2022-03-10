clear all 
clc
close all

% folder that contains algorithm scripts
addpath 'functions';
set(groot,'defaultFigureVisible','off') % 'on' to turn back on.  
tic

K = 9; %number of clusters


c2 = 1; % weight of simple edge
c3_array = 0:0.3:1;
gamma_array = 0.5; % gamma for likelihood plot
data_name = "contact-primary-school";
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
    %if  ~isfile(train_file) || ~isfile(test_file)
        [train_list, test_list, n] = SplitData(data_name,train_ratio);
        [W2, W3, T3, E2, E3] = LoadSimplex(train_list, n);
        save(strcat('processed_data/', data_name, '_train'),'W2', 'W3', 'T3', 'E2', 'E3');
        [W2, W3, T3, E2, E3] = LoadSimplex(test_list, n);
        save(strcat('processed_data/', data_name, '_test'),'W2', 'W3', 'T3', 'E2', 'E3');
    %end
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
%idx_rand_in = randperm(n_nodes);% shuffle the nodes
%[~, idx_reverse] = sort(idx_rand_in);
W2_in = train.W2;
W3_in = train.W3;
T3_in = train.T3;

for ii = 1:length(c3_array)
    c3 = c3_array(ii);
    norm_eta = c2*n_edge + c3*n_triangle;

    
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

    
    %find indices of open triangles
    trg_open_idx = find(1-train.T3(:));
    
    %triangle
    [X,Y,T,AUC_triangle_lin(end+1)] = perfcurve(test.T3(trg_open_idx),P_triangles_lin(trg_open_idx),1,'XCrit','reca','YCrit','prec');
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
    [X,Y,T,AUC_triangle_per(end+1)] = perfcurve(test.T3(trg_open_idx),P_triangles_per(trg_open_idx),1,'XCrit','reca','YCrit','prec');
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

end

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

score_arith = CalculateArithMean(test.W3);
score_harm = CalculateHarmMean(test.W3);
score_geom = CalculateGeoMean(test.W3);


[X,Y,T,AUC_triangle_arith] = perfcurve(test.T3(trg_open_idx),score_arith(trg_open_idx),1,'XCrit','reca','YCrit','prec');
[X,Y,T,AUC_triangle_harm] = perfcurve(test.T3(trg_open_idx),score_harm(trg_open_idx),1,'XCrit','reca','YCrit','prec');
[X,Y,T,AUC_triangle_geo] = perfcurve(test.T3(trg_open_idx),score_geom(trg_open_idx),1,'XCrit','reca','YCrit','prec');
[X,Y,T,AUC_triangle_rand] = perfcurve(test.T3(trg_open_idx), rand(size(trg_open_idx)),1,'XCrit','reca','YCrit','prec');

[X,Y,T,AUC_edge_rand] = perfcurve(test.W2(:), rand(n*n,1),1,'XCrit','reca','YCrit','prec');
[X,Y,T,AUC_edge_mean] = perfcurve(test.W2(:), test.W2(:),1,'XCrit','reca','YCrit','prec');


toc
%load handel
%sound(y,Fs)