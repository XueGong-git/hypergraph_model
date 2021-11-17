tic
clear
set(groot,'defaultFigureVisible','off') % 'on' to turn back on.  


n = 100;
K = 5; % number of clusters 
m = n/K; % number of nodes per cluster
a = 0.05;
%a = 0.1*2*pi/K; % noise;


c2 = 1; % weight of simple edge
c3 = 1/3; % weight of triangles
gamma_input = [1 3 5 10]; % gamma for generating graph
%gamma_array = gamma_input;
gamma_array = 0:0.25:4; % gamma for likelihood plot

rand_linear = [];
rand_periodic = [];


for input = 1
    for gamma = gamma_input
        
        switch input
    
        case 1 
            x = linspace(0,2,n);%uniform from 0 to 2pi
            data_type = "uniform";
        case 2 
            x = sort(repmat(linspace(0,2, K),1,m)+(2*a*rand(1,n)-a));
            %x = sort(repmat(linspace(-pi,pi,K),1,m)+(2*a*rand(1,n)-a)); % angles from -pi to pi
            data_type = "cluster";
        case 3
            x = linspace(0,0,n); %overlapping x
            data_type = "overlap";
        end

        
        [W2, W3, T3] = GenerateLinearHypergraph(x, gamma, c2, c3, data_type);
        

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
        
        %shuffle input adjacency matrix
        idx_rand = randperm(size(W2,1));% shuffle the nodes
        W2 = W2(idx_rand,idx_rand);
        W3 = W3(idx_rand,idx_rand);
        [~, idx_reverse] = sort(idx_rand);
        
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
        
       %plot estimated embedding
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
        lnP_linear = []; lnP0_linear = [];
        lnP_periodic = []; lnP0_periodic = [];

        for test_gamma = gamma_array
            [lnP_linear(end+1), lnP0_linear(end+1)] = CalculateModelLikelihood(x_est_linear, W2, T3, c2, c3, test_gamma, "linear");
            [lnP_periodic(end+1), lnP0_periodic(end+1)] = CalculateModelLikelihood(x_est_periodic, W2, T3, c2, c3, test_gamma, "periodic");
            
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
        %plot(gamma_array, lnP0_linear, '--b', 'LineWidth',1.5);
        %plot(gamma_array, lnP0_periodic, '--r', 'LineWidth',1.5);
        plot(gamma_max_linear, lnP_linear(max_linear_idx), 'ok', 'MarkerSize',10, 'LineWidth',2);
        plot(gamma_max_periodic, lnP_periodic(max_periodic_idx), 'or', 'MarkerSize',10, 'LineWidth',2);
        xline(gamma,'-',{'True \gamma'},'fontsize',20)
        legend({'Linear','Periodic','MLE'},'FontSize', 20,'Location','southeast');
        xlabel('\gamma','FontSize', 13);
        ylabel('Log-likelihood','FontSize', 13);
        set(gca,'fontsize',30);
        set(gca,'XLim',[0 max(gamma_array)])
        plt.LineWidth = 2;
        ax = gca;
        exportgraphics(ax,strcat('plots/model_comparison_linear_', data_type,'_gamma_', num2str(round(gamma,2)),'.eps'),'Resolution',300) 
        hold off;

        if data_type == "cluster"
            
            cluster_input = kmeans(transpose(x), K);
            cluster_est_linear = kmeans(x_est_linear, K);
            cluster_est_periodic = kmeans(x_est_periodic, K);
            rand_linear(end+1) = CalculateRandIndex(cluster_input, cluster_est_linear);
            rand_periodic(end+1) = CalculateRandIndex(cluster_input, cluster_est_periodic);

        end
    
        
    end
    
             %perform K-mean and compare with input


end

% plot rand_index

if data_type == "cluster"

    plt = plot(gamma_input, rand_linear, 'b', 'LineWidth',1.5);
    hold on;
    plot(gamma_input, rand_periodic, '-r', 'LineWidth',1.5);
    legend({'Linear','Periodic'},'FontSize', 20,'Location','southeast');
    xlabel('\gamma','FontSize', 13);
    ylabel('Rand Index','FontSize', 13);
    set(gca,'fontsize',30);
    set(gca,'YLim',[0.5 1.1])
    plt.LineWidth = 2;
    ax = gca;
    exportgraphics(ax,strcat('plots/rand_linear_', data_type,'.eps'),'Resolution',300) 
    hold off;
end

toc
