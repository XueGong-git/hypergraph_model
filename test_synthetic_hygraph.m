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


for input_shape = "periodic" % "periodic" or "linear"
    
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
    %{
            % plot L2 eigenvalues
            % eigenvalues of L2 and L3
            d2 = sum(W2, 2); D2 = diag(d2);
            d3 = sum(W3, 2); D3 = diag(d3);
            [V2, lambda2] = eigs(D2-W2,n,'smallestabs'); % all eigenvectors ranging from smallest eigenvalue to largest eigenvalue
            [V3, lambda3] = eigs(D3-W3,n,'smallestabs'); % all eigenvectors ranging from smallest eigenvalue to largest eigenvalue


            eigenvalues2 = diag(lambda2);
            scatter(1:min(n, 30), eigenvalues2(1:min(n, 30)), 100, 'MarkerFaceColor','black');
            xlabel('i','FontSize', 13);
            ylabel('Eigenvalue \lambda_i','FontSize', 13);
            set(gca,'fontsize',30);
            ax = gca;
            exportgraphics(ax,strcat('plots/',data_type,'_L2_eigenvalues.eps'),'Resolution',300) 

            % plot L3 eigenvalues
            eigenvalues3 = diag(lambda3);
            scatter(1:min(n, 30), eigenvalues3(1:min(n, 30)), 100, 'MarkerFaceColor','black');
            xlabel('i','FontSize', 13);
            ylabel('Eigenvalue \lambda_i','FontSize', 13);
            set(gca,'fontsize',30);
            ax = gca;
            exportgraphics(ax,strcat('plots/',data_type,'_L3_eigenvalues.eps'),'Resolution',300) 

            if size(V2,2)>=4
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
                exportgraphics(t,strcat('plots/',data_type,'_L2_eigenvectors.eps'),'Resolution',300) 
            end

            if size(V3,2)>=4

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
                exportgraphics(t,strcat('plots/',data_type,'_L3_eigenvectors.eps'),'Resolution',300) 
            end

    %}


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

            
            
            %{
           %plot estimated embedding

            s = scatter(x, x_est_linear, 200, 'MarkerFaceColor','black','MarkerEdgeColor','none');
            alpha(s,0.3) % transparent color
            xlabel('x','FontSize', 13);
            ylabel('x*','FontSize', 13);
            set(gca,'fontsize',30);
            ax = gca;
            exportgraphics(ax,strcat('plots/',input_shape,'_',data_type,'_linear_embedding_gamma=', num2str(round(gamma,2)),'.eps'),'Resolution',300) 


            s = scatter(x, x_est_periodic, 200, 'MarkerFaceColor','black','MarkerEdgeColor','none');
            alpha(s,0.3) % transparent color
            xlabel('x','FontSize', 13);
            ylabel('x*','FontSize', 13);
            set(gca,'fontsize',30);
            ax = gca;
            exportgraphics(ax,strcat('plots/',input_shape,'_',data_type,'_periodic_embedding_gamma=', num2str(round(gamma,2)),'.eps'),'Resolution',300) 




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
    'rand_linear', 'rand_periodic',  'max_lnP_linear', 'max_linear_idx', ...
'max_lnP_linear_scaled', 'max_linear_scaled_idx', 'max_lnP_periodic', 'max_periodic_idx', ...
'n_edge', 'n_nodes', 'n_triangle');








toc
beep on
beep

%load handel
%sound(y,Fs)