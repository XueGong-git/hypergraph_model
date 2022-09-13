clear all 
clc
close all

%%check how many triangles in test data are formed by triplets with geom mean
%%= 0
train_ratio = 0.8;

% folder that contains algorithm scripts
addpath 'functions';
set(groot,'defaultFigureVisible','off') % 'on' to turn back on.  
tic

%parameter for synthetic graph
K = 9; %number of clusters
gamma_inp = 18; % decay parameter for synthetic networks
c2_inp = 1;
c3_inp = 1/3;

c2 = 1; % weight of diadic edge
c3_array = 0.3; %0:0.1:1; % weight of triadic edge
gamma_array = 10:4:30; % gamma for likelihood plot
data_name = "contact-high-school";
train_file = strcat('processed_data/', data_name,'_train_', num2str(train_ratio),'.mat');
test_file = strcat('processed_data/', data_name,'_test_',num2str(1-train_ratio),'.mat');
m = 10; % number of nodes per cluster

if data_name == "synthetic_lin" || data_name == "synthetic_per"
    n_nodes = K*m; % total number of nodes
end


% likelihood
lnP_linear = zeros(length(c3_array),length(gamma_array));
lnP_periodic = zeros(length(c3_array),length(gamma_array));


% load data
if ismember (data_name, ["contact-high-school", "contact-primary-school", "coauth-DBLP", "email-Eu", "email-Enron" ] )
    if  ~isfile(train_file) || ~isfile(test_file)
        [train_list, test_list, n_nodes] = SplitData(data_name,train_ratio); %split data into train and test by time
        [W2, W3, T3, E2, E3] = LoadSimplex(train_list, n_nodes); %save list of simplices into matrix
        save(train_file,'W2', 'W3', 'T3', 'E2', 'E3');
        [W2, W3, T3, E2, E3] = LoadSimplex(test_list, n_nodes);
        save(test_file,'W2', 'W3', 'T3', 'E2', 'E3');
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



score_arith = CalculateArithMean(train.W3);
score_harm = CalculateHarmMean(train.W3);
score_geom = CalculateGeoMean(train.W3);



% how many of these triangles have zero mean geometric mean
triangle_idx = find(test.T3(:)==1);
zero_geom_idx = find(score_geom == 0);
[val,pos]=intersect(triangle_idx,zero_geom_idx);
length(pos)
length(pos)/length(triangle_idx)

load handel
sound(y,Fs)
