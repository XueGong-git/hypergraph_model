tic
clear


n = 100;
K = 5; % number of clusters 
m = n/K; % number of node per cluster
a = 0.2*m; % noise;


c2 = 1/2; % weight of simple edge
c3 = 1/3; % weight of simple edge



for input = 2
    for gamma = 0.025
        
        switch input
    
        case 1 
            x = linspace(1,2*pi,n);
            data_type = "uniform";
        case 2 
            x = sort(repmat(linspace(1,2*pi,K),1,m)+(2*a*rand(1,n)-a)); % trophic levels
            data_type = "cluster";
        case 3
            x = linspace(0,0,n); %overlapping x
            data_type = "overlap";
        end

        
        [W2, W3, T3] = GenerateLinearHypergraph(x, gamma, c2, c3, data_type);
        
        data_type = "highschool";

        figure
        % plot original W2 and W3
        imagesc(W2,[0,1]); %plot color map of original matrix
        colormap(flipud(gray(2)));
        set(gca,'FontSize',30) ;
        ax = gca;% Requires R2020a or later
        exportgraphics(ax,strcat('plots/',data_type,'_W2_input.eps'),'Resolution',300) 


        imagesc(W3); % plot color map of original matrix
        colormap(flipud(gray(256)));colorbar
        set(gca,'FontSize',30) ;
        set(gca,'ColorScale','log')
        ax = gca;% Requires R2020a or later
        exportgraphics(ax,strcat('plots/',data_type,'_W3_input_.eps'),'Resolution',300) 


        % plot L2 eigenvalues
        % eigenvalues of L2 and L3
        d2 = sum(W2, 2); D2 = diag(d2);
        d3 = sum(W3, 2); D3 = diag(d3);
        [V2, lambda2] = eigs(D2-W2,n,'smallestabs'); % all eigenvectors ranging from smallest eigenvalue to largest eigenvalue
        [V3, lambda3] = eigs(D3-W3,n,'smallestabs'); % all eigenvectors ranging from smallest eigenvalue to largest eigenvalue


        eigenvalues2 = diag(lambda2);
        scatter(1:40, eigenvalues2(1:40), 100, 'MarkerFaceColor','black');
        xlabel('i','FontSize', 13);
        ylabel('Eigenvalue \lambda_i','FontSize', 13);
        set(gca,'fontsize',30);
        ax = gca;
        exportgraphics(ax,strcat('plots/',data_type,'_L2_eigenvalues.eps'),'Resolution',300) 

        % plot L3 eigenvalues
        eigenvalues3 = diag(lambda3);
        scatter(1:40, eigenvalues3(1:40), 100, 'MarkerFaceColor','black');
        xlabel('i','FontSize', 13);
        ylabel('Eigenvalue \lambda_i','FontSize', 13);
        set(gca,'fontsize',30);
        ax = gca;
        exportgraphics(ax,strcat('plots/',data_type,'_L3_eigenvalues.eps'),'Resolution',300) 

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


        
        %shuffle input adjacency matrix
        idx_rand = randperm(size(W2,1));% shuffle the nodes
        W2 = W2(idx_rand,idx_rand);
        W3 = W3(idx_rand,idx_rand);
        
        %estimate embedding using linear spectral clustering
        [x_est_linear] = LinearHypergraphEmbedding(W2, W3, c2, c3, "false");
        [x_est_periodic] = PeriodicHypergraphEmbedding(W2, W3, c2, c3, "false");

                
        
        
        % plot estimated embedding
        figure
        s = scatter(x(idx_rand), x_est_linear, 200, 'MarkerFaceColor','black','MarkerEdgeColor','none');
        alpha(s,0.3) % transparent color
        xlabel('x','FontSize', 13);
        ylabel('x*','FontSize', 13);
        set(gca,'fontsize',30);
        ax = gca;
        exportgraphics(ax,strcat('plots/linear_hygraph_embedding_', data_type,'_gamma=', num2str(round(gamma,2)),'.eps'),'Resolution',300) 

        figure
        s = scatter(x(idx_rand), x_est_periodic, 200, 'MarkerFaceColor','black','MarkerEdgeColor','none');
        alpha(s,0.3) % transparent color
        xlabel('x','FontSize', 13);
        ylabel('x*','FontSize', 13);
        set(gca,'fontsize',30);
        ax = gca;
        exportgraphics(ax,strcat('plots/periodic_hygraph_embedding_', data_type,'_gamma=', num2str(round(gamma,2)),'.eps'),'Resolution',300) 

        %compare likelihood
        lnP_linear = [];
        lnP_periodic = [];
        gamma_array = [0:0.005:0.025]
        for test_gamma = gamma_array
            lnP_linear(end+1) = CalculateModelLikelihood(x, W2, T3, c3, test_gamma, "linear");
            lnP_periodic(end+1) = CalculateModelLikelihood(x, W2, T3, c3, test_gamma, "periodic");
        end
        
        %maximum likelihood of gamma
        [~, max_linear_idx] = max(lnP_linear);
        gamma_max_linear = test_gamma(max_linear_idx); 
        [~, max_periodic_idx] = max(lnP_periodic);
        gamma_max_periodic = test_gamma(max_periodic_idx); 

        
        % plot likelihood
        plt = plot(gamma_array, lnP_linear, 'LineWidth',1.5);
        hold on;
        plot(gamma_array, lnP_periodic, '--r', 'LineWidth',1.5);
        plot(gamma_max_linear, lnP_linear(max_linear_idx), 'ok', 'MarkerSize',10, 'LineWidth',2);
        plot(gamma_max_periodic, lnP_periodic(max_periodic_idx), 'or', 'MarkerSize',10, 'LineWidth',2);
        legend({'Linear','Periodic', 'MLE'},'FontSize', 20,'Location','southwest');
        xlabel('\gamma','FontSize', 13);
        ylabel('Log-likelihood','FontSize', 13);
        set(gca,'fontsize',30);
        plt.LineWidth = 2;
        ax = gca;
        exportgraphics(ax,strcat('plots/model_comparison_input_', num2str(input),'.eps'),'Resolution',300) 
        hold off;

    end
end
toc
