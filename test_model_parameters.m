clear all 
clc
close all

% folder that contains algorithm scripts
addpath 'functions';
set(groot,'defaultFigureVisible','off') % 'on' to turn back on.  
tic


n = 100;
K = 5; % number of clusters 
ntrial = 40;

c2 = 1; % weight of simple edge
gamma_array = 0.5; % gamma for likelihood plot
gamma_input = 0.5; % gamma for generating graph
c3_inp_array = 0:0.1:0.4; % weight of triangles for input graph
c3_inv_array = 0:0.1:0.4; % weight of triangles for inverse problem

for input_shape = "linear"
    
for input = 2
    rand_linear = zeros(length(c3_inp_array), length(c3_inv_array),length(gamma_input));
    rand_periodic = zeros(size(rand_linear));
    max_lnP_linear = zeros(size(rand_linear));
    max_lnP_linear_scaled = zeros(size(rand_linear));
    max_lnP_periodic = zeros(size(rand_linear));
    correlation_lin = zeros(size(rand_linear)); %correlation between estimated and truth
    correlation_per = zeros(size(rand_linear)); %correlation between estimated and truth
    triangle_density =  zeros(length(c3_inp_array),length(gamma_input));
    edge_density =  zeros(length(c3_inp_array),length(gamma_input));
    for ii = 1:length(c3_inp_array)
        c3_inp = c3_inp_array(ii);
        for kk = 1:length(gamma_input)
                gamma = gamma_input(kk);
                
                %repeated trials
                for iterator = 1:ntrial

                % generate inputs
                [x, W2, W3, T3, data_type] = GenerateHygraph(n, K, gamma, c2, c3_inp, input_shape, input);

                %calculate edge density
                n_nodes = size(W2,2); n_edge = sum(W2,'all')/2; n_triangle = sum(T3, 'all')/6;
                edge_density(ii, kk) = edge_density(ii, kk) + 2*n_edge/(n_nodes*(n_nodes-1)*ntrial);
                triangle_density(ii, kk) = triangle_density(ii, kk) + 6*n_triangle/(n_nodes*(n_nodes-1)*(n_nodes-2)*ntrial);
                
                %extract largest connected component
                [x, W2, W3, T3] = MaxConnectedSubgraph(x, c2, c3_inp, W2, W3, T3);    
                n_nodes = size(W2,2);
                
                %plot the adjacency matrix
                imagesc(W2,[0,1]); %plot color map of original matrix
                colormap(flipud(gray(2)));
                set(gca,'FontSize',30) ;
                ax = gca;% Requires R2020a or later
                exportgraphics(ax,strcat('plots/',input_shape,'_',data_type,'_gamma=', num2str(gamma),'_input.eps'),'Resolution',300) 


                %plot the adjacency matrix
                imagesc(W3); %plot color map of original matrix
                colormap(flipud(gray(256)));
                set(gca,'FontSize',30) ;
                ax = gca;% Requires R2020a or later
                exportgraphics(ax,strcat('plots/',input_shape,'_',data_type,'_gamma=', num2str(gamma),'_input_W3_c3_',num2str(c3_inp),'.eps'),'Resolution',300) 


                %shuffle input adjacency matrix
                %idx_rand = randperm(size(W2,1));% shuffle the nodes
                %[~, idx_reverse] = sort(idx_rand); % index to undo the shuffle
                %W2 = W2(idx_rand,idx_rand); W3 = W3(idx_rand,idx_rand);
                %x = x(idx_rand);

            for jj = 1:length(c3_inv_array)
                c3_inv = c3_inv_array(jj);

                %estimate embedding using linear spectral clustering
                x_est_linear = LinearHypergraphEmbedding(W2, W3, c2, c3_inv, "false");
                x_est_periodic = PeriodicHypergraphEmbedding(W2, W3, c2, c3_inv, "false");
                
                %normalize the estimated embedding to the same range
                %x_est_linear = x_est_linear*norm(x,2)/norm(x_est_linear,2);        
                
                %calculate eta
                [~, eta_linear_est] = CalculateModelLikelihood(x_est_linear, W2, T3, c2, c3_inv, 2, "linear");
                [~, eta_linear_true] = CalculateModelLikelihood(x, W2, T3, c2, c3_inv, 2, "linear");

                %scale estimation
                x_est_linear = x_est_linear*sqrt(eta_linear_true/eta_linear_est);  

                %scale linear to the same incoherence
                corr_mat_lin = corrcoef(x_est_linear, x);
                corr_mat_per = corrcoef(x_est_periodic, x);
                correlation_lin(ii,jj,kk) = correlation_lin(ii,jj,kk) + abs(corr_mat_lin(1,2))/ntrial;
                correlation_per(ii,jj,kk) = correlation_per(ii,jj,kk) + abs(corr_mat_per(1,2))/ntrial;


                %reverse to the input order
                %x = x(idx_reverse);
                %x_est_linear = x_est_linear(idx_reverse);
                %x_est_periodic = x_est_periodic(idx_reverse);
                %W2 = W2(idx_reverse, idx_reverse);
                %W3 = W3(idx_reverse, idx_reverse);

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

                max_lnP_linear(ii,jj,kk) = max_lnP_linear(ii,jj,kk) + CalculateModelLikelihood(x_est_linear, W2, T3, c2, c3_inv, gamma, "linear")/ntrial;
                max_lnP_periodic(ii,jj,kk) = max_lnP_periodic(ii,jj,kk) + CalculateModelLikelihood(x_est_periodic, W2, T3, c2, c3_inv, gamma, "periodic")/ntrial;
        
                if data_type == "cluster"
                    cluster_input = kmeans(transpose(x), K);
                    cluster_est_linear = kmeans(x_est_linear, K);
                    cluster_est_periodic = kmeans(x_est_periodic, K);
                    rand_linear(ii,jj,kk) = rand_linear(ii,jj,kk) + CalculateRandIndex(cluster_input, cluster_est_linear)/ntrial;
                    rand_periodic(ii,jj,kk) = rand_periodic(ii,jj,kk) + CalculateRandIndex(cluster_input, cluster_est_periodic)/ntrial;
        
                end
            end
            end
        end
    end


%%%%%%%%%%%%%%%% SAVE DATA %%%%%%%%%%%%%%%%%%%%%%%%
save(strcat('mat/',input_shape,'_',data_type, '.mat'), 'rand_linear', 'rand_periodic',  'max_lnP_linear', 'max_lnP_periodic', 'edge_density', 'triangle_density');
%%%%%%%%%%%%%%%% PLOTS %%%%%%%%%%%%%%%%%

    for select_gamma_in = 1:length(gamma_input)
        gamma_plot = gamma_input(select_gamma_in);

        % plot heatmap of linear rand
        if data_type == "cluster"
        figure 
        h = heatmap(round(rand_linear(:,:,select_gamma_in), 2));
        %caxis([0, 1]);
        h.Title = strcat('Linear Rand, \gamma_{in}=', num2str(gamma_plot) );
        h.XLabel = 'c3_{inv}';
        h.YLabel = 'c3_{inp}';
        h.XData = round(c3_inv_array,2);
        h.YData = round(c3_inp_array,2);
        ax = gca;% Requires R2020a or later
        exportgraphics(ax,strcat('plots/',input_shape,'_',data_type,'_gamma=', num2str(gamma_plot),'_rand_lin.eps'),'Resolution',300) 
        
        % plot heatmap of periodic rand
        h = heatmap(round(rand_periodic(:,:,select_gamma_in), 2));
        %caxis([0, 1]);
        h.Title = strcat('Periodic Rand, \gamma_{in}=', num2str(gamma_plot) );
        h.XLabel = 'c3_{inv}';
        h.YLabel = 'c3_{inp}';
        h.XData = round(c3_inv_array,2);
        h.YData = round(c3_inp_array,2);
        ax = gca;% Requires R2020a or later
        exportgraphics(ax,strcat('plots/',input_shape,'_',data_type,'_gamma=', num2str(gamma_plot),'_rand_per.eps'),'Resolution',300) 
        
        % plot heatmap of linear - periodic rand
        h = heatmap(round(rand_linear(:,:,select_gamma_in)-rand_periodic(:,:,select_gamma_in), 2));
        %caxis([-0.5, 0.5]);
        h.Title = strcat('Linear Rand - Periodic Rand, \gamma_{in}=', num2str(gamma_plot) );
        h.XLabel = 'c3_{inv}';
        h.YLabel = 'c3_{inp}';
        h.XData = round(c3_inv_array,2);
        h.YData = round(c3_inp_array,2);
        ax = gca;% Requires R2020a or later
        exportgraphics(ax,strcat('plots/',input_shape,'_',data_type,'_gamma=', num2str(gamma_plot),'_rand_lin_vs_per.eps'),'Resolution',300) 
        
        end

        % plot heatmap of linear max log-likelihood
        h = heatmap(round(max_lnP_linear(:,:,select_gamma_in), 2),'Colormap',parula);
        h.CellLabelFormat = '%.2f';
        h.Title = strcat('Max LnP Linear, \gamma_{in}=', num2str(gamma_plot) );
        h.XLabel = 'c3_{inv}';
        h.YLabel = 'c3_{inp}';
        h.XData = round(c3_inv_array,2);
        h.YData = round(c3_inp_array,2);
        ax = gca;% Requires R2020a or later
        exportgraphics(ax,strcat('plots/',input_shape,'_',data_type,'_gamma=', num2str(gamma_plot),'_maxlnP_lin.eps'),'Resolution',300) 
        
        % plot heatmap of periodic max log-likelihood
        h = heatmap(round(max_lnP_periodic(:,:,select_gamma_in), 2),'Colormap',parula);
        h.Title = strcat('Max LnP Periodic, \gamma_{in}=', num2str(gamma_plot) );
        h.XLabel = 'c3_{inv}';
        h.YLabel = 'c3_{inp}';
        h.XData = round(c3_inv_array,2);
        h.YData = round(c3_inp_array,2);
        ax = gca;% Requires R2020a or later
        exportgraphics(ax,strcat('plots/',input_shape,'_',data_type,'_gamma=', num2str(gamma_plot),'_maxlnP_per.eps'),'Resolution',300) 
        
        % plot heatmap of linear - periodic log-likelihood
        h = heatmap(round(max_lnP_linear(:,:,select_gamma_in)-max_lnP_periodic(:,:,select_gamma_in), 2),'Colormap',parula);
        h.Title = strcat('Max LnP Linear-Periodic, \gamma_{in}=', num2str(gamma_plot) );
        h.XLabel = 'c3_{inv}';
        h.YLabel = 'c3_{inp}';
        h.XData = round(c3_inv_array,2);
        h.YData = round(c3_inp_array,2);
        ax = gca;% Requires R2020a or later
        exportgraphics(ax,strcat('plots/',input_shape,'_',data_type,'_gamma=', num2str(gamma_plot),'_maxlnP_lin_vs_per.eps'),'Resolution',300) 
      
        % plot heatmap of linear correlation
        h = heatmap(round(correlation_lin(:,:,select_gamma_in), 2),'Colormap',parula);
        h.Title = strcat('Correlation betweeen estimation and truth (linear) \gamma_{in}=', num2str(gamma_plot) );
        h.XLabel = 'c3_{inv}';
        h.YLabel = 'c3_{inp}';
        h.XData = round(c3_inv_array,2);
        h.YData = round(c3_inp_array,2);
        ax = gca;% Requires R2020a or later
        exportgraphics(ax,strcat('plots/',input_shape,'_',data_type,'_gamma=', num2str(gamma_plot),'_corr_lin.eps'),'Resolution',300) 
        
        % plot heatmap of periodic correlation
        h = heatmap(round(correlation_per(:,:,select_gamma_in), 2),'Colormap',parula);
        h.Title = strcat('Correlation betweeen estimation and truth (periodic) \gamma_{in}=', num2str(gamma_plot) );
        h.XLabel = 'c3_{inv}';
        h.YLabel = 'c3_{inp}';
        h.XData = round(c3_inv_array,2);
        h.YData = round(c3_inp_array,2);
        ax = gca;% Requires R2020a or later
        exportgraphics(ax,strcat('plots/',input_shape,'_',data_type,'_gamma=', num2str(gamma_plot),'_corr_per.eps'),'Resolution',300) 
     
        % plot heatmap of linear - periodic correlation
        h = heatmap(round(correlation_lin(:,:,select_gamma_in)-correlation_per(:,:,select_gamma_in), 2),'Colormap',parula);
        h.Title = strcat('Difference between correlations \gamma_{in}=', num2str(gamma_plot) );
        h.XLabel = 'c3_{inv}';
        h.YLabel = 'c3_{inp}';
        h.XData = round(c3_inv_array,2);
        h.YData = round(c3_inp_array,2);
        ax = gca;% Requires R2020a or later
        exportgraphics(ax,strcat('plots/',input_shape,'_',data_type,'_gamma=', num2str(gamma_plot),'_corr_lin_vs_per.eps'),'Resolution',300) 


        % dual axis
        plt = plot(c3_inp_array, triangle_density(:,select_gamma_in), '--k', 'LineWidth',1.5);
        hold on;
        plot(c3_inp_array, edge_density(:,select_gamma_in), '--b', 'LineWidth',1.5);
        legend({'Triangle density','Edge density'},'FontSize', 20,'Location','best');
        xlabel('c3_{inp}','FontSize', 13);
        ylabel('Density','FontSize', 13);
        set(gca,'fontsize',25);
        plt.LineWidth = 2;
        ax = gca;
        exportgraphics(ax,strcat('plots/',input_shape,'_',data_type,'_gamma=', num2str(gamma_plot),'_density.eps'),'Resolution',300) 
        hold off;
        
    end

end


end

toc
beep on
%load handel
%sound(y,Fs)