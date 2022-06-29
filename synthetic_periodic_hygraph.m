tic
clear


n = 100;
K = 5; % number of clusters 
m = n/K; % number of node per cluster
a = 0.05*pi; % noise;


c2 = 1/2; % weight of simple edge
c3 = 1/3; % weight of simple edge
gamma_input = 4; % gamma for generating graph
gamma_array = 0:0.25:8; % gamma for likelihood plot
rand_linear = [];
rand_periodic = [];



for input = 
    for gamma = gamma_input
        
        switch input
    
        case 1 
            x = linspace(0,2*pi,n);%uniform from 0 to 2pi
            data_type = "uniform";
        case 2 
            x = sort(repmat(linspace(-pi,pi,K),1,m)+(2*a*rand(1,n)-a)); % angles from -pi to pi
            data_type = "cluster";
        case 3
            x = linspace(0,0,n); %overlapping x
            data_type = "overlap";
        end

        
        [W2, W3, T3] = GeneratePeriodicHypergraph(x, gamma, c2, c3, data_type);
        
        figure
        % plot original W2 and W3
        imagesc(W2,[0,1]); %plot color map of original matrix
        colormap(flipud(gray(2)));
        set(gca,'FontSize',30) ;
        ax = gca;% Requires R2020a or later
        exportgraphics(ax,strcat('plots/periodic_',data_type,'_W2_input.eps'),'Resolution',300) 


        imagesc(W3); % plot color map of original matrix
        colormap(flipud(gray(256)));colorbar
        set(gca,'FontSize',30) ;
        set(gca,'ColorScale','log')
        ax = gca;% Requires R2020a or later
        exportgraphics(ax,strcat('plots/periodic_',data_type,'_W3_input_.eps'),'Resolution',300) 

        
        %shuffle input adjacency matrix
        idx_rand = randperm(size(W2,1));% shuffle the nodes
        [~, idx_reverse] = sort(idx_rand);
        W2 = W2(idx_rand,idx_rand);
        W3 = W3(idx_rand,idx_rand);

        %estimate embedding using linear spectral clustering
        [x_est_linear] = LinearHypergraphEmbedding(W2, W3, c2, c3, "false");
        [x_est_periodic] = PeriodicHypergraphEmbedding(W2, W3, c2, c3, "false");

        %normalize the estimated embedding to the same range
        x_est_linear = x_est_linear*norm(x,2)/norm(x_est_linear,2);        
        
        %reverse to the input order
        x_est_linear = x_est_linear(idx_reverse);
        x_est_periodic = x_est_periodic(idx_reverse);
        W2 = W2(idx_reverse, idx_reverse);
        W3 = W3(idx_reverse, idx_reverse);
        
        % plot estimated embedding
        figure
        s = scatter(x, x_est_linear, 200, 'MarkerFaceColor','black','MarkerEdgeColor','none');
        alpha(s,0.3) % transparent color
        xlabel('x','FontSize', 13);
        ylabel('x*','FontSize', 13);
        set(gca,'fontsize',30);
        ax = gca;
        exportgraphics(ax,strcat('plots/linear_hygraph_embedding_', data_type,'_gamma=', num2str(round(gamma,2)),'.eps'),'Resolution',300) 

        figure
        s = scatter(x, x_est_periodic, 200, 'MarkerFaceColor','black','MarkerEdgeColor','none');
        alpha(s,0.3) % transparent color
        xlabel('x','FontSize', 13);
        ylabel('x*','FontSize', 13);
        set(gca,'fontsize',30);
        ax = gca;
        exportgraphics(ax,strcat('plots/periodic_hygraph_embedding_', data_type,'_gamma=', num2str(round(gamma,2)),'.eps'),'Resolution',300) 

        %compare likelihood
        lnP_linear = []; lnP0_linear = [];lnP_linear_test = [];
        lnP_periodic = []; lnP0_periodic = [];
        
         
        
        for test_gamma = gamma_array
            [lnP_linear(end+1), lnP0_linear(end+1)] = CalculateModelLikelihood(x_est_linear, W2, T3, c2, c3, test_gamma, "linear");
            [lnP_periodic(end+1), lnP0_periodic(end+1)] = CalculateModelLikelihood(x_est_periodic, W2, T3, c2, c3, test_gamma, "periodic");
            [lnP_linear_test(end+1), ~] = CalculateModelLikelihood(x_est_linear*1.5, W2, T3, c2, c3, test_gamma, "linear");

        end
        
        %maximum likelihood of gamma
        [~, max_linear_idx] = max(lnP_linear);
        gamma_max_linear = gamma_array(max_linear_idx); 
        [~, max_periodic_idx] = max(lnP_periodic);
        gamma_max_periodic = gamma_array(max_periodic_idx); 

        
        % plot likelihood
        plt = plot(gamma_array, lnP_linear, 'b', 'LineWidth',1.5);
        hold on;
        plot(gamma_array, lnP_periodic, '-r', 'LineWidth',1.5);
        %plot(gamma_array, lnP_linear_test, '-k', 'LineWidth',1.5);
        %plot(gamma_array, lnP0_linear, '--b', 'LineWidth',1.5);
        %plot(gamma_array, lnP0_periodic, '--r', 'LineWidth',1.5);
        plot(gamma_max_linear, lnP_linear(max_linear_idx), 'ok', 'MarkerSize',10, 'LineWidth',2);
        plot(gamma_max_periodic, lnP_periodic(max_periodic_idx), 'or', 'MarkerSize',10, 'LineWidth',2);
        xline(gamma,'-',{'True \gamma'},'fontsize',20)
        legend({'Linear','Periodic','MLE'},'FontSize', 20,'Location','southwest');
        xlabel('\gamma','FontSize', 13);
        ylabel('Log-likelihood','FontSize', 13);
        set(gca,'fontsize',30);
        plt.LineWidth = 2;
        ax = gca;
        exportgraphics(ax,strcat('plots/model_comparison_periodic_', data_type,'_gamma_', num2str(round(gamma,2)),'.eps'),'Resolution',300) 
        hold off;
        
        
         %perform K-mean and compare with input
        if data_type == "cluster"
            cluster_input = kmeans(transpose(x), K);
            cluster_est_linear = kmeans(x_est_linear, K);
            cluster_est_periodic = kmeans(x_est_periodic, K);
            rand_linear(end+1) = CalculateRandIndex(cluster_input, cluster_est_linear);
            rand_periodic(end+1) = CalculateRandIndex(cluster_input, cluster_est_periodic);


        end

    end
    
    
if data_type == "cluster"

    % plot rand_index
    plt = plot(gamma_input, rand_linear, 'b', 'LineWidth',1.5);
    hold on;
    plot(gamma_input, rand_periodic, '-r', 'LineWidth',1.5);
    legend({'Linear','Periodic'},'FontSize', 20,'Location','southeast');
    xlabel('\gamma','FontSize', 13);
    ylabel('Rand Index','FontSize', 13);
    set(gca,'fontsize',30);
    plt.LineWidth = 2;
    ax = gca;
    exportgraphics(ax,strcat('plots/rand_periodic_', data_type,'.eps'),'Resolution',300) 
    hold off;
end
    
end
toc
sound(sin(1:3000));
