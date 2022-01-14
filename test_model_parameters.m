clear all 
clc
close all

% folder that contains algorithm scripts
addpath 'functions';
set(groot,'defaultFigureVisible','on') % 'on' to turn back on.  
tic


n = 50;
K = 5; % number of clusters 


c2 = 1; % weight of simple edge
gamma_array = 0:1:10; % gamma for likelihood plot
gamma_input = 0:1:10; % gamma for generating graph
c3_inp_array = 0:0.1:0.9; % weight of triangles for input graph
c3_inv_array = 0:0.1:0.9; % weight of triangles for inverse problem

for input_shape = "periodic"
    
for input = 2
    rand_linear = zeros(length(c3_inp_array), length(c3_inv_array),length(gamma_input));
    rand_periodic = zeros(size(rand_linear));
    max_lnP_linear = zeros(size(rand_linear));
    max_lnP_linear_scaled = zeros(size(rand_linear));
    max_lnP_periodic = zeros(size(rand_linear));
    triangle_density =  zeros(length(c3_inp_array),length(gamma_input));
    edge_density =  zeros(length(c3_inp_array),length(gamma_input));

    for ii = 1:length(c3_inp_array)
        c3_inp = c3_inp_array(ii);
        for kk = 1:length(gamma_input)
                gamma = gamma_input(kk);

                % generate inputs
                [x, W2, W3, T3, data_type] = GenerateHygraph(n, K, gamma, c2, c3_inp, input_shape, input);
                
                %calculate edge density
                n_nodes = size(W2,2); n_edge = sum(W2,'all')/2; n_triangle = sum(T3, 'all')/6;
                edge_density(ii, kk) = 2*n_edge/(n_nodes*(n_nodes-1));
                triangle_density(ii, kk) = 6*n_edge/(n_nodes*(n_nodes-1)*(n_nodes-2));
                
                %extract largest connected component
                [x, W2, W3, T3] = MaxConnectedSubgraph(x, c2, c3_inp, W2, W3, T3);    
                n_nodes = size(W2,2);
                
                %shuffle input adjacency matrix
                idx_rand = randperm(size(W2,1));% shuffle the nodes
                [~, idx_reverse] = sort(idx_rand); % index to undo the shuffle
                W2 = W2(idx_rand,idx_rand); W3 = W3(idx_rand,idx_rand);
                x = x(idx_rand);

            for jj = 1:length(c3_inv_array)
                c3_inv = c3_inv_array(jj);

                %estimate embedding using linear spectral clustering
                x_est_linear = LinearHypergraphEmbedding(W2, W3, c2, c3_inv, "false");
                x_est_periodic = PeriodicHypergraphEmbedding(W2, W3, c2, c3_inv, "false");
                
                %normalize the estimated embedding to the same range
                x_est_linear = x_est_linear*norm(x,2)/norm(x_est_linear,2);        
                
                %calculate eta
                [~, eta_linear] = CalculateModelLikelihood(x_est_linear, W2, T3, c2, c3_inv, 2, "linear");
                [~, eta_periodic] = CalculateModelLikelihood(x_est_periodic, W2, T3, c2, c3_inv, 2, "periodic");
            
                %scale linear to the same incoherence
                x_est_linear_scaled = x_est_linear*eta_periodic/eta_linear;  
                
                %compare likelihood
                lnP_linear = zeros(1, length(gamma_array));  lnP_linear_scaled = zeros(1, length(gamma_array)); 
                lnP_periodic = zeros(1, length(gamma_array)); 
                

                for ll = 1:length(gamma_array)
                    test_gamma = gamma_array(ll);
                    [lnP_linear(ll), ~] = CalculateModelLikelihood(x_est_linear, W2, T3, c2, c3_inv, test_gamma, "linear");
                    [lnP_linear_scaled(ll), ~] = CalculateModelLikelihood(x_est_linear_scaled, W2, T3, c2, c3_inv, test_gamma, "linear");
                    [lnP_periodic(ll), ~] = CalculateModelLikelihood(x_est_periodic, W2, T3, c2, c3_inv, test_gamma, "periodic");
                end
                
                %maximum likelihood of gamma
                max_lnP_linear(ii,jj,kk) = max(lnP_linear);
                max_lnP_linear_scaled(ii,jj,kk) = max(lnP_linear_scaled);
                max_lnP_periodic(ii,jj,kk) = max(lnP_periodic);
        
                if data_type == "cluster"
                    cluster_input = kmeans(transpose(x), K);
                    cluster_est_linear = kmeans(x_est_linear, K);
                    cluster_est_periodic = kmeans(x_est_periodic, K);
                    rand_linear(ii,jj,kk) = CalculateRandIndex(cluster_input, cluster_est_linear);
                    rand_periodic(ii,jj,kk) = CalculateRandIndex(cluster_input, cluster_est_periodic);
        
                end
                
            end
        end
    end


%%%%%%%%%%%%%%%% SAVE DATA %%%%%%%%%%%%%%%%%%%%%%%%
save(strcat('mat/',input_shape,'_',data_type, '.mat'), 'rand_linear', 'rand_periodic',  'max_lnP_linear', 'max_lnP_periodic', 'edge_density', 'triangle_density');
%%%%%%%%%%%%%%%% PLOTS %%%%%%%%%%%%%%%%%

select_gamma_in = 6;
gamma_plot = gamma_input(select_gamma_in);

% plot heatmap of linear rand

figure 
h = heatmap(round(rand_linear(:,:,select_gamma_in), 2));
h.Title = strcat('Linear Rand, \gamma_{in}=', num2str(gamma_plot) );
h.XLabel = 'c3_{inp}';
h.YLabel = 'c3_{inv}';
h.XData = round(c3_inp_array,2);
h.YData = round(c3_inv_array,2);
ax = gca;% Requires R2020a or later
exportgraphics(ax,strcat('plots/',input_shape,'_',data_type,'gamma=', num2str(gamma_plot),'_rand_lin.eps'),'Resolution',300) 

% plot heatmap of periodic rand
h = heatmap(round(rand_periodic(:,:,select_gamma_in), 2));
h.Title = strcat('Periodic Rand, \gamma_{in}=', num2str(gamma_plot) );
h.XLabel = 'c3_{inp}';
h.YLabel = 'c3_{inv}';
h.XData = round(c3_inp_array,2);
h.YData = round(c3_inv_array,2);
ax = gca;% Requires R2020a or later
exportgraphics(ax,strcat('plots/',input_shape,'_',data_type,'gamma=', num2str(gamma_plot),'_rand_per.eps'),'Resolution',300) 

% plot heatmap of linear - periodic rand
h = heatmap(round(rand_linear(:,:,select_gamma_in)-rand_periodic(:,:,select_gamma_in), 2),'Colormap',parula);
h.Title = strcat('Linear Rand - Periodic Rand, \gamma_{in}=', num2str(gamma_plot) );
h.XLabel = 'c3_{inp}';
h.YLabel = 'c3_{inv}';
h.XData = round(c3_inp_array,2);
h.YData = round(c3_inv_array,2);
ax = gca;% Requires R2020a or later
exportgraphics(ax,strcat('plots/',input_shape,'_',data_type,'gamma=', num2str(gamma_plot),'_rand_lin_vs_per.eps'),'Resolution',300) 


% plot heatmap of linear max log-likelihood

 
h = heatmap(round(max_lnP_linear(:,:,select_gamma_in), 2));
h.Title = strcat('Max LnP Linear, \gamma_{in}=', num2str(gamma_plot) );
h.XLabel = 'c3_{inp}';
h.YLabel = 'c3_{inv}';
h.XData = round(c3_inp_array,2);
h.YData = round(c3_inv_array,2);
ax = gca;% Requires R2020a or later
exportgraphics(ax,strcat('plots/',input_shape,'_',data_type,'gamma=', num2str(gamma_plot),'_maxlnP_lin.eps'),'Resolution',300) 

% plot heatmap of periodic max log-likelihood
h = heatmap(round(max_lnP_periodic(:,:,select_gamma_in), 2));
h.Title = strcat('Max LnP Periodic, \gamma_{in}=', num2str(gamma_plot) );
h.XLabel = 'c3_{inp}';
h.YLabel = 'c3_{inv}';
h.XData = round(c3_inp_array,2);
h.YData = round(c3_inv_array,2);
ax = gca;% Requires R2020a or later
exportgraphics(ax,strcat('plots/',input_shape,'_',data_type,'gamma=', num2str(gamma_plot),'_maxlnP_per.eps'),'Resolution',300) 

% plot heatmap of linear - periodic log-likelihood
h = heatmap(round(max_lnP_linear(:,:,select_gamma_in)-max_lnP_periodic(:,:,select_gamma_in), 2),'Colormap',parula);
h.Title = strcat('Max LnP Linear-Periodic, \gamma_{in}=', num2str(gamma_plot) );
h.XLabel = 'c3_{inp}';
h.YLabel = 'c3_{inv}';
h.XData = round(c3_inp_array,2);
h.YData = round(c3_inv_array,2);
ax = gca;% Requires R2020a or later
exportgraphics(ax,strcat('plots/',input_shape,'_',data_type,'gamma=', num2str(gamma_plot),'_maxlnP_lin_vs_per.eps'),'Resolution',300) 


end
end
toc
load handel
sound(y,Fs)