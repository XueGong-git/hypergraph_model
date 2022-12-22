% Predict triangles in real data
%
% INPUTS
% - W3 Projected matrix which counts common triangles 
% 
%
% OUTPUTS
% - S_triangle A vector of scores for  triangles
%


function [AUC_triangle_lin, AUC_triangle_geo_shuffle] = PredictTriangles(data_name, train_ratio, c3_array, gamma_array, trim, n_eig)


c2 = 1; % weight of diadic edge
train_file = strcat('processed_data/', data_name,'_train_', num2str(train_ratio),'.mat');
test_file = strcat('processed_data/', data_name,'_test_',num2str(1-train_ratio),'.mat');

if data_name == "synthetic_lin" || data_name == "synthetic_per" % synthetic data set
    K = 4; %number of clusters
    m = 60; % number of nodes per cluster
    n_nodes = K*m; % total number of nodes
end

% initialize variables
max_lnP_linear = []; % find max lnP for linear embedding
max_lnP_periodic = [];  % find max lnP for linear embedding
AUC_triangle_lin = []; % area under curve with linear model
AUC_triangle_per = []; % area under curve for periodic model
lnP_linear = zeros(length(c3_array),length(gamma_array)); % lnP for linear embedding
lnP_periodic = zeros(length(c3_array),length(gamma_array)); % lnP for linear embedding
gamma_max_linear = [];
gamma_max_periodic = [];

% load data
if ismember (data_name, ["contact-high-school", "contact-primary-school", "email-Eu", "email-Enron" ] )
    if  ~isfile(train_file) || ~isfile(test_file)
        [train_list, test_list, n_nodes] = SplitData(data_name,train_ratio); % split data into train and test by time
        [W2, W3, T3, E2, E3] = LoadSimplex(train_list, n_nodes); % save list of simplices into matrix
        save(train_file,'W2', 'W3', 'T3', 'E2', 'E3');
        [W2, W3, T3, E2, E3] = LoadSimplex(test_list, n_nodes); % save list of simplices into matrix
        save(test_file,'W2', 'W3', 'T3', 'E2', 'E3');
    end
    train = load(train_file,'W2', 'W3', 'T3', 'E2', 'E3'); % load training data
    test  = load(test_file,'W2', 'W3', 'T3', 'E2', 'E3'); % load testing data
    
elseif data_name == "synthetic_lin" % generates synthetic linear graph

    a = 0.01; % noise parameter
    x = sort(repmat(linspace(0,2,K),1,m)+(2*a*rand(1,n_nodes)-a)); % input node embedding
    [E2, W2, W3, T3] = GenerateLinearHypergraph(x, gamma_inp, c2_inp, c3_inp, 1, data_name);
    [train_list, test_list, n_nodes] = SplitData(data_name,train_ratio); %split data into train and test by time
    [W2, W3, T3, ~] = LoadSimplex(train_list, n_nodes); %save list of simplices into matrix
    save(train_file,'W2', 'W3', 'T3', 'E2');
    [W2, W3, T3, ~] = LoadSimplex(test_list, n_nodes);
    save(test_file,'W2', 'W3', 'T3', 'E2');

elseif data_name == "synthetic_per" % generates synthetic periodic data
    
    a = 0.01*pi; % noise parameter
    x = sort(repmat(linspace(-pi,pi,K),1,m)+(2*a*rand(1,n_nodes)-a)); % angles from -pi to pi
    [train_list, test_list, n_nodes] = SplitData(data_name,train_ratio); %split data into train and test by time
    [W2, W3, T3, ~] = LoadSimplex(train_list, n_nodes); %save list of simplices into matrix
    save(train_file,'W2', 'W3', 'T3', 'E2', 'E3');
    [W2, W3, T3, ~] = LoadSimplex(test_list, n_nodes);
    save(test_file,'W2', 'W3', 'T3', 'E2', 'E3');
    
end

W2_in = train.W2;
W3_in = train.W3;
T3_in = train.T3;

if trim == "true" % trim off nodes with high or low degree
    % check degree distribution
    degree_eff = sum(W2_in + W3_in, 2);

    % trim  nodes with highest degrees
    keepIDs = (degree_eff > quantile(degree_eff,0.02));
    %keepIDs = (degree_eff > quantile(degree_eff,0.02))&(degree_eff< quantile(degree_eff,0.98));
    W2_in = W2_in(keepIDs, keepIDs);
    T3_in = T3_in(keepIDs, keepIDs, keepIDs);
    W3_in = sum(T3_in, 3);

   
    % trim test data accordingly
    test.W2 = test.W2(keepIDs, keepIDs);
    test.T3 = test.T3(keepIDs, keepIDs, keepIDs);
    test.W3 = sum(test.T3, 3);
    
end

% extract largest connected component
[idx, W2_in, W3_in, T3_in] = MaxConnectedSubgraph(1, 1, W2_in, W3_in, T3_in);    

% get only the nodes in the connected component in the training data
test.W2 = test.W2(idx,idx);
test.W3 = test.W3(idx,idx);
test.T3 = test.T3(idx,idx, idx); 

% check number of nodes and edge/triangle density
n_nodes = size(train.W2,2);
n_edge = sum(train.W2,'all')/2;
n_triangle = sum(train.T3, 'all')/6;
edge_density = 2*n_edge/(n_nodes*(n_nodes-1));
triangle_density = 6*n_triangle/(n_nodes*(n_nodes-1)*(n_nodes-2));


%plot the input adjacency matrix
imagesc(train.W2,[0,1]); %plot color map of original matrix
colormap(flipud(gray(2)));
set(gca,'FontSize',30) ;
ax = gca;% Requires R2020a or later
exportgraphics(ax,strcat('plots/',data_name,'_input_W2.eps'),'Resolution',300) 

%plot the input triangle adjacency matrix
imagesc(train.W3); %plot color map of original matrix
colormap(flipud(gray(256)));
set(gca,'FontSize',30) ;
ax = gca;% Requires R2020a or later
exportgraphics(ax,strcat('plots/',data_name,'_input_W3.eps'),'Resolution',300) 


for ii = 1:length(c3_array)
    c3 = c3_array(ii);
    
    %estimate embedding using linear spectral clustering 
    [x_est_linear] = LinearHypergraphEmbedding(W2_in, W3_in, c2, c3, "false", n_eig); 
    [x_est_periodic] = PeriodicHypergraphEmbedding(W2_in, W3_in, c2, c3, "false");
    
    %calculate eta
    [~, eta_linear_est, ~] = CalculateModelLikelihood(x_est_linear, W2_in, T3_in, c2, c3, 2, "linear");
    [~, eta_periodic_est, ~] = CalculateModelLikelihood(x_est_periodic, W2_in, T3_in, c2, c3, 2, "linear");

    %scale estimation
    x_est_linear = x_est_linear*sqrt(eta_periodic_est/eta_linear_est);  
    
    %x_est_periodic = x_est_periodic*sqrt(norm_eta/eta_periodic_est);  
    [~, eta_linear_scaled, ~] = CalculateModelLikelihood(x_est_linear, W2_in, T3_in, c2, c3, 2, "linear");
    %[~, eta_periodic_scaled, ~] = CalculateModelLikelihood(x_est_linear, W2_in, T3_in, c2, c3, 2, "linear");


    %reorder nodes according to the embedding
    [~,idx_linear] = sort(x_est_linear);
    W2_reorder_linear = W2_in(idx_linear,idx_linear);
    W3_reorder_linear = W3_in(idx_linear,idx_linear);
    
    %reorder nodes according to the embedding
    [~,idx_periodic] = sort(x_est_periodic);
    W2_reorder_periodic = W2_in(idx_periodic,idx_periodic);
    W3_reorder_periodic = W3_in(idx_periodic,idx_periodic);
    
    
    for jj = 1:length(gamma_array) % return likelihood for each gamma test points
        gamma = gamma_array(jj);
        [lnP_linear(ii,jj), ~] = CalculateModelLikelihood(x_est_linear, W2_in, T3_in, c2, c3, gamma, "linear");
        [lnP_periodic(ii,jj), ~] = CalculateModelLikelihood(x_est_periodic, W2_in, T3_in, c2, c3, gamma, "periodic");
    end
    
    % maximum likelihood of gamma
    [max_lnP_linear(end+1), max_linear_idx] = max(lnP_linear(ii,:));
    gamma_max_linear(end+1) = gamma_array(max_linear_idx); 
    [max_lnP_periodic(end+1), max_periodic_idx] = max(lnP_periodic(ii,:));
    gamma_max_periodic(end+1) = gamma_array(max_periodic_idx); 
    
   
    % plot likelihood
    plt = plot(gamma_array, lnP_linear(ii,:), 'b', 'LineWidth',1.5);
    hold on;
    plot(gamma_array, lnP_periodic(ii,:), '-r', 'LineWidth',1.5);
    plot(gamma_max_linear(end), lnP_linear(ii, max_linear_idx), 'ok', 'MarkerSize',10, 'LineWidth',2);
    plot(gamma_max_periodic(end), lnP_periodic(ii, max_periodic_idx), 'or', 'MarkerSize',10, 'LineWidth',2);
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

% find c3 and gamma that gives maximum likelihood
% linear embedding
[~, max_c3_lin_idx] = max(max_lnP_linear);
max_c3_lin = c3_array(max_c3_lin_idx); 
max_gamma_lin = gamma_max_linear(max_c3_lin_idx);

% periodic embedding
[~, max_c3_per_idx] = max(max_lnP_periodic);
max_c3_per = c3_array(max_c3_per_idx); 
max_gamma_per = gamma_max_periodic(max_c3_per_idx);


% plot max likelihood on training data 
plt = plot(c3_array, max_lnP_linear, 'b', 'LineWidth',1.5);
hold on;
plot(c3_array, max_lnP_periodic, '-r', 'LineWidth',1.5);
plot(max_c3_lin, max(max_lnP_linear), 'or', 'MarkerSize',10, 'LineWidth',2);
plot(max_c3_per, max(max_lnP_periodic), 'or', 'MarkerSize',10, 'LineWidth',2);
legend({'Linear','Periodic'},'FontSize', 20,'Location','best');
xlabel('c_3','FontSize', 13);
ylabel('Max log likelihood','FontSize', 13);
set(gca,'fontsize',30);
ax = gca;
exportgraphics(ax,strcat('plots/max_lnP_', data_name,'_c2_',num2str(round(c2,2)),'.eps'),'Resolution',300) 
hold off;

% test AUC with c3 that gives rise to highest likelihood

% calculate likelihood of edges and triangles with, actually the gamma
% chosen shouldn't affect the AUC
[~, ~, P_edge_lin, P_triangles_lin] = CalculateModelLikelihood(x_est_linear, W2_in, T3_in, c2, max_c3_lin, max_gamma_lin, "linear");

%find indices of open triangles
trg_open_idx = find(1-T3_in(:));

%triangle
[X_lin,Y_lin,T_lin,AUC_triangle_lin(end+1)] = perfcurve(test.T3(trg_open_idx),P_triangles_lin(trg_open_idx),1,'XCrit','reca','YCrit','prec');
plot(X_lin,Y_lin)
xlabel('Recall') 
ylabel('Precision')
title('ROC for triangle prediction (linear)')
set(gca,'fontsize',30);
ax = gca;
exportgraphics(ax,strcat('plots/ROC_triangle_lin_', data_name,'_c2_',num2str(round(c2,2)),'_c3_',num2str(round(c3,2)),'.eps'),'Resolution',300) 

% periodic embedding
[~, ~, P_edge_per, P_triangles_per] = CalculateModelLikelihood(x_est_periodic, W2_in, T3_in, c2, max_c3_per, max_gamma_per, "periodic");

%triangle
[X_per,Y_per,T_per,AUC_triangle_per(end+1)] = perfcurve(test.T3(trg_open_idx),P_triangles_per(trg_open_idx),1,'XCrit','reca','YCrit','prec');
plot(X_per,Y_per)
xlabel('Recall') 
ylabel('Precision')
title('ROC for triangle prediction (periodic)')
set(gca,'fontsize',30);
ax = gca;
exportgraphics(ax,strcat('plots/ROC_triangle_per_', data_name,'_c2_',num2str(round(c2,2)),'_c3_',num2str(round(c3,2)),'.eps'),'Resolution',300) 

%%%%%% ROC for mean scores %%%%%%%%%
score_arith = CalculateArithMean(W3_in);
score_harm = CalculateHarmMean(W3_in);
score_geom = CalculateGeoMean(W3_in);

%%%%%%%%% arithmetic mean %%%%%%%%%%%
[X_arith,Y_arith,T_arith,AUC_triangle_arith] = perfcurve(test.T3(trg_open_idx),score_arith(trg_open_idx),1,'XCrit','reca','YCrit','prec','TVals', linspace(max(score_arith(trg_open_idx), [], 'all'),1/3));
arith_open = score_arith(trg_open_idx);
arith_sort = sort(arith_open); [~,~,min_nonzero] = find(arith_sort,1); % find the smallest non-zero
arith_zero_idx = find(arith_open==0); % index with zero score
arith_open(arith_zero_idx) = rand(length(arith_zero_idx),1)*min_nonzero; %random number from 0 to smallest non-zero

%plot new ROC curve
[X_harm,Y_harm,T_harm,AUC_triangle_arith_shuffle] = perfcurve(test.T3(trg_open_idx),arith_open,1,'XCrit','reca','YCrit','prec');

%%%%%%%%%% harmonic mean %%%%%%%%%%%%
[X_harm,Y_harm,T_harm,AUC_triangle_harm] = perfcurve(test.T3(trg_open_idx),score_harm(trg_open_idx),1,'XCrit','reca','YCrit','prec','TVals', linspace(max(score_harm(trg_open_idx), [], 'all'),1));
% assign random score from 0 to 1 to triplets with score_harm = 0
harm_open = score_harm(trg_open_idx);
harm_sort = sort(harm_open); [~,~,min_nonzero] = find(harm_sort,1); % find the smallest non-zero
harm_zero_idx = find(harm_open==0); % index with zero score
harm_open(harm_zero_idx) = rand(length(harm_zero_idx),1)*min_nonzero; %random number from 0 to smallest non-zero
%plot new ROC curve
[X_harm,Y_harm,T_harm,AUC_triangle_harm_shuffle] = perfcurve(test.T3(trg_open_idx),harm_open,1,'XCrit','reca','YCrit','prec');


%%%%%%%% geometric mean %%%%%%%%%
[X_geo,Y_geo,T_geo, AUC_triangle_geo] = perfcurve(test.T3(trg_open_idx),score_geom(trg_open_idx),1,'XCrit','reca','YCrit','prec','TVals', linspace(max(score_geom(trg_open_idx), [], 'all'),1));
% assign random score from 0 to 1 to triplets with score_harm = 0
geom_open = score_geom(trg_open_idx);
geo_sort = sort(geom_open); [~,~,min_nonzero] = find(geo_sort,1); % find the smallest non-zero
geo_zero_idx = find(geom_open==0); % index with zero score
geom_open(geo_zero_idx) = rand(length(geo_zero_idx),1)*min_nonzero; %random number from 0 to smallest non-zero
%plot new ROC curve
[X_geo,Y_geo,T_geo, AUC_triangle_geo_shuffle] = perfcurve(test.T3(trg_open_idx),geom_open,1,'XCrit','reca','YCrit','prec');


%%%%%% random score %%%%%%%%%%
[X_rand,Y_rand,T_rand, AUC_triangle_rand] = perfcurve(test.T3(trg_open_idx), rand(size(trg_open_idx)),1,'XCrit','reca','YCrit','prec');

% among test data,  how many of triangles have zero mean geometric mean
triangle_idx = find(test.T3(:)==1);
zero_geom_idx = find(score_geom == 0);
[val,pos]=intersect(triangle_idx,zero_geom_idx);
length(pos);
length(pos)/length(triangle_idx);

% scatter plot of AUC
plot(X_lin,Y_lin,'Marker','s')
hold on
plot(X_per,Y_per,'Marker','s')
plot(X_geo,Y_geo,'Marker','s')
plot(X_harm,Y_harm,'Marker','s')
hold off
xlabel('Recall')
ylabel('Precision')
title('ROC curve')
legend('Linear model', 'Geometric Mean', 'Harmonic Mean')


% plot area under curve
area(X_lin,Y_lin, 'FaceAlpha', 0.8);
hold on
area(X_geo,Y_geo, 'FaceAlpha', 0.8);
area(X_harm,Y_harm, 'FaceAlpha', 0.8);
ax = gca;
exportgraphics(ax,strcat('plots/ROC_comparison_', data_name,'_c2_',num2str(round(c2,2)),'_c3_',num2str(round(c3,2)),'.eps'),'Resolution',300) 

% plot area under curve
area(X_geo,Y_geo, 'FaceAlpha', 0.8);
hold on
area(X_lin,Y_lin, 'FaceAlpha', 0.8);
legend('Geometric score', 'Linear model')
xlabel('Recall') 
ylabel('Precision')
ax = gca;
set(gca,'fontsize',20)
exportgraphics(ax,strcat('plots/ROC_comparison_', data_name,'_train_',num2str(round(train_ratio,2)),'.eps'),'Resolution',300) 
hold off
save(strcat("results/AUC_parameters_",data_name, "_train=", num2str(train_ratio), '.mat'), 'c3_array', 'gamma_array', 'gamma_max_periodic', 'max_c3_lin', 'max_gamma_lin', ...
    'max_lnP_linear', 'max_lnP_periodic','gamma_max_linear', 'max_c3_per', 'max_gamma_per');

T = table(AUC_triangle_rand, AUC_triangle_lin, AUC_triangle_per, AUC_triangle_arith, AUC_triangle_arith_shuffle, AUC_triangle_geo, AUC_triangle_geo_shuffle, ...
    AUC_triangle_harm,AUC_triangle_harm_shuffle);
writetable(T,strcat("results/AUC_",data_name, "_train=", num2str(train_ratio), '.dat'));

end