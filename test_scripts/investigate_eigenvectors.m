clear all 
clc
close all

% folder that contains algorithm scripts
addpath 'functions';
set(groot,'defaultFigureVisible','off') % 'on' to turn back on.  
tic

%parameter for synthetic graph
K = 9; %number of clusters
gamma_inp = 18; % decay parameter for synthetic networks
c2_inp = 1;
c3_inp = 1;

c2 = 1; % weight of diadic edge
c3_array = 0:0.2:1; % weight of triadic edge
gamma_array = 0:2:10; % gamma for likelihood plot
data_name = "contact-high-school";
train_ratio = 0.8;
train_file = strcat('processed_data/', data_name,'_train.mat');
test_file = strcat('processed_data/', data_name,'_test.mat');
m = 10; % number of nodes per cluster

if data_name == "synthetic_lin" || data_name == "synthetic_per"
    n_nodes = K*m; % total number of nodes
end

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


% load data
if ismember (data_name, ["contact-high-school", "contact-primary-school", "coauth-DBLP", "email-Eu", "email-Enron" ] )
    if  ~isfile(train_file) || ~isfile(test_file)
        [train_list, test_list, n_nodes] = SplitData(data_name,train_ratio); %split data into train and test by time
        [W2, W3, T3, E2, E3] = LoadSimplex(train_list, n_nodes); %save list of simplices into matrix
        save(strcat('processed_data/', data_name, '_train'),'W2', 'W3', 'T3', 'E2', 'E3');
        [W2, W3, T3, E2, E3] = LoadSimplex(test_list, n_nodes);
        save(strcat('processed_data/', data_name, '_test'),'W2', 'W3', 'T3', 'E2', 'E3');
    end
    train = load(train_file,'W2', 'W3', 'T3', 'E2', 'E3');
    test  = load(test_file,'W2', 'W3', 'T3', 'E2', 'E3');
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
        x = sort(repmat(linspace(0,2, K),1,m)+(2*a*rand(1,n_nodes)-a)); % input node embedding
        [E2, E3, W2, W3, T3] = GenerateLinearHypergraph(x, gamma_inp, c2_inp, c3_inp);   
        save(strcat('processed_data/', data_name, '_train'),'W2', 'W3', 'T3', 'E2', 'E3');
        [E2, E3, W2, W3, T3] = GenerateLinearHypergraph(x, gamma_inp, c2_inp, c3_inp);   
        save(strcat('processed_data/', data_name, '_test'),'W2', 'W3', 'T3', 'E2', 'E3');
        train = load(train_file,'W2', 'W3', 'T3', 'E2', 'E3');
        test  = load(test_file,'W2', 'W3', 'T3', 'E2', 'E3');
    %end
elseif data_name == "synthetic_per"
    a = 0.01*pi; % noise parameter
    x = sort(repmat(linspace(-pi,pi,K),1,m)+(2*a*rand(1,n_nodes)-a)); % angles from -pi to pi
    [W2, W3, T3] = GeneratePeriodicHypergraph(x, gamma_inp, c2_inp, c3_inp);
    label = repmat(linspace(1,K,K),1,m);

end




% check number of nodes and edge/triangle density
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
              
    %plot linear embedding
    scatter(1:length(x_est_linear), sort(x_est_linear))
    ax = gca;
    exportgraphics(ax,strcat('plots/embedding_lin_', data_name,'_c3_',num2str(round(c3,2)),'.eps'),'Resolution',300) 

    
    %plot periodic embedding
    scatter(sin(sort(x_est_periodic)), cos(sort(x_est_periodic)))
    xlim([-1 1])
    ylim([-1 1])
    ax = gca;
    exportgraphics(ax,strcat('plots/embedding_per_', data_name,'_c3_',num2str(round(c3,2)),'.eps'),'Resolution',300) 

    %calculate eta
    [~, eta_linear_est, ~] = CalculateModelLikelihood(x_est_linear, W2_in, T3_in, c2, c3, 2, "linear");
    [~, eta_periodic_est, ~] = CalculateModelLikelihood(x_est_periodic, W2_in, T3_in, c2, c3, 2, "linear");

    %scale estimation
    x_est_linear = x_est_linear*sqrt(norm_eta/eta_linear_est);  
    x_est_periodic = x_est_periodic*sqrt(norm_eta/eta_periodic_est);  

    [~, eta_linear_scaled, ~] = CalculateModelLikelihood(x_est_linear, W2_in, T3_in, c2, c3, 2, "linear");
    [~, eta_periodic_scaled, ~] = CalculateModelLikelihood(x_est_linear, W2_in, T3_in, c2, c3, 2, "linear");
    

end


load handel
sound(y,Fs)