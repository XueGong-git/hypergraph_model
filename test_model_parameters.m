clear all 
clc
close all

% folder that contains algorithm scripts
addpath 'functions';
set(groot,'defaultFigureVisible','off') % 'on' to turn back on.  
tic


n = 250;
K = 5; % number of clusters 
ntrial = 40;

c2 = 1; % weight of simple edge
gamma_array = 1; % gamma for likelihood plot
gamma_input = 1; % gamma for generating graph
c3_inp_array = 0:0.1:1; % weight of triangles for input graph
c3_inv_array = 0:0.1:1; % weight of triangles for inverse problem

for input_shape = ["linear", "periodic"]
    
for input = 2
    rand_linear = zeros(length(c3_inp_array), length(c3_inv_array), ntrial);
    rand_periodic = zeros(size(rand_linear));
    max_lnP_linear = zeros(size(rand_linear));
    max_lnP_linear_scaled = zeros(size(rand_linear));
    max_lnP_periodic = zeros(size(rand_linear));
    correlation_lin = zeros(size(rand_linear)); %correlation between estimated and truth
    correlation_per = zeros(size(rand_linear)); %correlation between estimated and truth
    triangle_density =  zeros(length(c3_inp_array),ntrial);
    edge_density =  zeros(length(c3_inp_array),ntrial);
    for ii = 1:length(c3_inp_array)
        c3_inp = c3_inp_array(ii);
        disp('************************')
        disp(['c3_inp = ', num2str(c3_inp)])

        for kk = 1:ntrial
                gamma = gamma_input;

                disp(['Trial ', num2str(kk)])
                %disp(['Gamma ', num2str(gamma)])

                % generate inputs
                [x, W2, W3, T3, data_type, cluster_input] = GenerateHygraph(n, K, gamma, c2, c3_inp, input_shape);
                
                
                %calculate edge density
                n_nodes = size(W2,2); n_edge = sum(W2,'all')/2; n_triangle = sum(T3, 'all')/6;
                edge_density(ii, kk) = edge_density(ii, kk) + 2*n_edge/(n_nodes*(n_nodes-1)*ntrial);
                triangle_density(ii, kk) = triangle_density(ii, kk) + 6*n_triangle/(n_nodes*(n_nodes-1)*(n_nodes-2)*ntrial);
                %sprintf('Number of nodes: \n %f \n', n_nodes)

                
                %sort nodes
                [x_sort, idx_sort] = sort(x);
                %plot the adjacency matrix
                imagesc(W2(idx_sort,idx_sort),[0,1]); %plot color map of original matrix
                colormap(flipud(gray(2)));
                set(gca,'FontSize',30) ;
                ax = gca;% Requires R2020a or later
                exportgraphics(ax,strcat('plots/',input_shape,'_gamma=', num2str(gamma),'_input_W2.eps'),'Resolution',300) 


                %plot the adjacency matrix
                imagesc(W3(idx_sort, idx_sort)); %plot color map of original matrix
                colormap(flipud(gray(256)));
                set(gca,'FontSize',30) ;
                ax = gca;% Requires R2020a or later
                colorbar
                exportgraphics(ax,strcat('plots/',input_shape,'_gamma=', num2str(gamma),'_input_W3_c3_',num2str(c3_inp),'.eps'),'Resolution',300) 

                %}


            for jj = 1:length(c3_inv_array)
                c3_inv = c3_inv_array(jj);
                %disp(['c3_inv = ', num2str(c3_inv)])

                %estimate embedding using linear spectral clustering
                x_est_linear = LinearHypergraphEmbedding(W2, W3, c2, c3_inv, "false", 1);
                x_est_periodic = PeriodicHypergraphEmbedding(W2, W3, c2, c3_inv, "false");
                
                %normalize the estimated embedding to the same range
                %x_est_linear = x_est_linear*norm(x,2)/norm(x_est_linear,2);        
                
                %calculate eta
                [~, eta_linear_est] = CalculateModelLikelihood(x_est_linear, W2, T3, c2, c3_inv, 2, "linear");
                [~, eta_linear_true] = CalculateModelLikelihood(x', W2, T3, c2, c3_inv, 2, "linear");

                %scale estimation
                x_est_linear = x_est_linear*sqrt(eta_linear_true/eta_linear_est);  
                
                %calculate correlation
                corr_mat_lin = corrcoef(x_est_linear, x);
                corr_mat_per = corrcoef(x_est_periodic, x);
                correlation_lin(ii,jj,kk) = abs(corr_mat_lin(1,2));
                correlation_per(ii,jj,kk) = abs(corr_mat_per(1,2));



                %compare likelihood
                lnP_linear = zeros(1, length(gamma_array));  lnP_linear_scaled = zeros(1, length(gamma_array)); 
                lnP_periodic = zeros(1, length(gamma_array)); 
                

                %for ll = 1:length(gamma_array)
                %    test_gamma = gamma_array(ll);
                %    [lnP_linear(ll), ~] = CalculateModelLikelihood(x_est_linear, W2, T3, c2, c3_inv, test_gamma, "linear");
                   % [lnP_linear_scaled(ll), ~] = CalculateModelLikelihood(x_est_linear_scaled, W2, T3, c2, c3_inv, test_gamma, "linear");
                %    [lnP_periodic(ll), ~] = CalculateModelLikelihood(x_est_periodic, W2, T3, c2, c3_inv, test_gamma, "periodic");
                %end
                
                %maximum likelihood of gamma
                %max_lnP_linear(ii,jj,kk) = max_lnP_linear(ii,jj,kk) + max(lnP_linear)/ntrial;
                %max_lnP_linear_scaled(ii,jj,kk) = max_lnP_linear_scaled(ii,jj,kk) + max(lnP_linear_scaled)/ntrial;
                %max_lnP_periodic(ii,jj,kk) = max_lnP_periodic(ii,jj,kk) + max(lnP_periodic)/ntrial;

                max_lnP_linear(ii,jj,kk) = CalculateModelLikelihood(x_est_linear, W2, T3, c2, c3_inv, gamma, "linear");
                max_lnP_periodic(ii,jj,kk) = CalculateModelLikelihood(x_est_periodic, W2, T3, c2, c3_inv, gamma, "periodic");
        
                if data_type == "cluster"
                    %cluster_input = kmeans(transpose(x), K); %% USE TRUE CLUSTER?%%
                    cluster_est_linear = kmeans(x_est_linear, K);
                    cluster_est_periodic = kmeans(x_est_periodic, K);
                    rand_linear(ii,jj,kk) = CalculateRandIndex(cluster_input, cluster_est_linear, 'adjusted');
                    rand_periodic(ii,jj,kk) = CalculateRandIndex(cluster_input, cluster_est_periodic, 'adjusted');
        
                end
            end
            toc
        end
        toc
        end
    


%%%%%%%%%%%%%%%% SAVE DATA %%%%%%%%%%%%%%%%%%%%%%%%
save(strcat(input_shape,'_',data_type, '_gamma=', num2str(gamma_input,2),'.mat'), 'rand_linear', ...
    'rand_periodic',  'max_lnP_linear', 'max_lnP_periodic', 'edge_density', 'triangle_density', ...
    'correlation_lin', 'correlation_per', 'n','K','ntrial','c2','gamma_array', 'gamma_input' ,...
    'c3_inp_array', 'c3_inv_array')


end


end

toc
%beep on
%beep
load handel
sound(y,Fs)