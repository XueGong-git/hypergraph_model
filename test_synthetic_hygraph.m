tic
clear
set(groot,'defaultFigureVisible','off') % 'on' to turn back on.  
addpath ('functions')

m = 50; %number of nodes per cluster
K = 5; % number of clusters 
n = m*K; % number of nodes
c2 = 1; % weight of simple edge
c3 = 1/3; % weight of triangles
gamma_array = 0:1:10; % gamma for likelihood plot
gamma_input = 0:1:10; % gamma for generating graph
ntrial = 40;


for input_shape = "linear" % "periodic" or "linear"
    
for data_type = "cluster"
    rand_linear = [];
    rand_periodic = [];
    triangle_density = [];
    edge_density = [];
    max_lnP_linear = [];
    max_lnP_linear_scaled = [];
    max_lnP_periodic = [];
    
    for gg = 1:length(gamma_input)
        gamma = gamma_input(gg);    
        
        disp('************************')
        disp(['gamma_in = ', num2str(gamma)])

        %%%% run multiple iterations
        for trial_id = 1:ntrial
            disp(['Trial ', num2str(trial_id)])
     
            % generate inputs
            [x, W2, W3, T3, data_type, cluster_input] = GenerateHygraph(n, K, gamma, c2, c3, input_shape);

            %calculate edge density
            n_nodes(trial_id, gg) = size(W2,2);
            n_edge(trial_id, gg) = sum(W2,'all')/2;
            n_triangle(trial_id, gg) = sum(T3, 'all')/6;

            W2_old = W2;
            W3_old = W3;

            figure
    


            %estimate embedding using linear spectral clustering
            [x_est_linear] = LinearHypergraphEmbedding(W2, W3, c2, c3, "false", 1);
            [x_est_periodic] = PeriodicHypergraphEmbedding(W2, W3, c2, c3, "false");

            %normalize the estimated embedding to the same range
            x_est_linear = x_est_linear*norm(x,2)/norm(x_est_linear,2);        

            %reorder nodes according to the embedding
            [~,idx_periodic] = sort(x_est_periodic);
            W2_reorder_periodic = W2(idx_periodic,idx_periodic);
            W3_reorder_periodic = W3(idx_periodic,idx_periodic);


            %%%%%%%%%%% plots %%%%%%%
 
%{            
            %%%%%% input adjacency Dyadic matrix w2 %%%%%%
            cla(gca,'reset')
            imagesc(W2,[0,1]); %plot color map of original matrix
            colormap(flipud(gray(2)));
            set(gca,'FontSize',30) ;
            ax = gca;% Requires R2020a or later
            exportgraphics(ax,strcat('plots/',input_shape,'_',data_type,'_W2_inp.eps'),'Resolution',300) 

            %%%%%% input adjacency Dyadic matrix w3 %%%%%%
            cla(gca,'reset')
            imagesc(W3); %plot color map of original matrix
            colormap(flipud(gray(256)));
            colorbar
            set(gca,'FontSize',30) ;
            ax = gca;% Requires R2020a or later
            exportgraphics(ax,strcat('plots/',input_shape,'_',data_type,'_W3_inp.eps'),'Resolution',300) 
%}
            
            %{

            % plot reordered W2
            imagesc(W2_reorder_linear,[0,1]); %plot color map of original matrix
            colormap(flipud(gray(2)));
            set(gca,'FontSize',30) ;
            ax = gca;% Requires R2020a or later
            exportgraphics(ax,strcat('plots/',input_shape,'_',data_type,'_W2_reorder_linear.eps'),'Resolution',300) 

            % plot reordered W2
            imagesc(W2_reorder_periodic,[0,1]); %plot color map of original matrix
            colormap(flipud(gray(2)));
            set(gca,'FontSize',30) ;
            ax = gca;% Requires R2020a or later
            exportgraphics(ax,strcat('plots/',input_shape,'_',data_type,'_W2_reorder_periodic.eps'),'Resolution',300) 
    %}

            %compare likelihood
            lnP_linear = [];  
            lnP_periodic = []; 
            lnP_linear_scaled = []; 

            %calculate eta
            [~, eta_linear] = CalculateModelLikelihood(x_est_linear, W2, T3, c2, c3, 2, "linear");
            [~, eta_periodic] = CalculateModelLikelihood(x_est_periodic, W2, T3, c2, c3, 2, "periodic");


            %scale linear to the same incoherence
            x_est_linear_scaled = x_est_linear*eta_periodic/eta_linear;  
            for test_gamma = gamma_array
                [lnP_linear(end+1), ~] = CalculateModelLikelihood(x_est_linear, W2, T3, c2, c3, test_gamma, "linear");
                [lnP_linear_scaled(end+1), ~] = CalculateModelLikelihood(x_est_linear_scaled, W2, T3, c2, c3, test_gamma, "linear");
                [lnP_periodic(end+1), ~] = CalculateModelLikelihood(x_est_periodic, W2, T3, c2, c3, test_gamma, "periodic");
            end

            %maximum likelihood
            [max_lnP_linear(trial_id, gg), max_linear_idx(gg, trial_id)] = max(lnP_linear);
            [max_lnP_linear_scaled(trial_id, gg), max_linear_scaled_idx(gg, trial_id)] = max(lnP_linear_scaled);
            [max_lnP_periodic(trial_id, gg), max_periodic_idx(gg, trial_id)] = max(lnP_periodic);




            if data_type == "cluster"
                cluster_est_linear = kmeans(x_est_linear, K);
                cluster_est_periodic = kmeans(x_est_periodic, K);
                rand_linear(trial_id, gg) = CalculateRandIndex(cluster_input, cluster_est_linear, 'adjusted');
                rand_periodic(trial_id, gg) = CalculateRandIndex(cluster_input, cluster_est_periodic, 'adjusted');

            end

            %{
            % plot likelihood
            plt = plot(gamma_array, lnP_linear, 'b', 'LineWidth',1.5);
            hold on;
            plot(gamma_array, lnP_linear_scaled, 'b--', 'LineWidth',1.5);
            plot(gamma_array, lnP_periodic, 'r', 'LineWidth',1.5);
            plot(gamma_max_linear, lnP_linear(max_linear_idx), 'ob', 'MarkerSize',10, 'LineWidth',2);
            plot(gamma_array, -1*eta_linear*gamma_array, ':b', 'LineWidth',1);
            plot(gamma_array, -1*eta_periodic*gamma_array, ':r', 'LineWidth',1);
            plot(gamma_max_periodic, lnP_periodic(max_periodic_idx), 'or', 'MarkerSize',10, 'LineWidth',2);
            plot(gamma_max_linear_scaled, lnP_linear_scaled(max_linear_scaled_idx), 'ob', 'MarkerSize',10, 'LineWidth',2);
            xline(gamma,'-',{'True \gamma'},'fontsize',20)
            legend({'Linear','Linear scaled','Periodic','MLE', 'Slope = -Periodic Incoherence', 'Slope = -Linear Incoherence'},'FontSize', 20,'Location','best');
            xlabel('\gamma','FontSize', 13);
            ylabel('Log-likelihood','FontSize', 13);
            set(gca,'fontsize',30);
            set(gca,'XLim',[0 max(gamma_array)])
            plt.LineWidth = 2;
            ax = gca;
            exportgraphics(ax,strcat('plots/',input_shape,'_',data_type, '_model_comparison_gamma_', num2str(round(gamma,2)),'.eps'),'Resolution',300) 
            hold off;
        %}
        toc
        end
    toc
    end
    
end
end

save(strcat(input_shape,'_',data_type, '_c3=', num2str(c3,2) ,'.mat'), ...
    'c2','c3', 'rand_linear', 'rand_periodic',  'max_lnP_linear', 'max_linear_idx', ...
'max_lnP_linear_scaled', 'max_linear_scaled_idx', 'max_lnP_periodic', 'max_periodic_idx', ...
'n_edge', 'n_nodes', 'n_triangle', 'gamma_array', 'gamma_input');








toc
beep on
beep

%load handel
%sound(y,Fs)