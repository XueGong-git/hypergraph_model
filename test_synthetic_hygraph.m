tic
clear
set(groot,'defaultFigureVisible','off') % 'on' to turn back on.  

n = 50;
K = 5; % number of clusters 
m = n/K; % number of nodes per cluster


c2 = 1; % weight of simple edge
c3 = 1/3; % weight of triangles
gamma_array = 0:1:10; % gamma for likelihood plot
gamma_input = 0:1:10; % gamma for generating graph

for input_shape = "periodic"
    
for input = 2
    rand_linear = [];
    rand_periodic = [];
    triangle_density = [];
    edge_density = [];
    max_lnP_linear = [];
    max_lnP_linear_scaled = [];
    max_lnP_periodic = [];
    for gamma = gamma_input
        % add more iterations
        if input_shape == "linear" 
            a = 0.05;

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
            [W2, W3, T3] = GenerateLinearHypergraph(x, gamma, c2, c3);

        elseif input_shape == "periodic" 
            a = 0.05*pi; % noise;

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
            [W2, W3, T3] = GeneratePeriodicHypergraph(x, gamma, c2, c3);

        end
        
        
        %calculate edge density
        n_nodes = size(W2,2);
        n_edge = sum(W2,'all')/2;
        n_triangle = sum(T3, 'all')/6;
        edge_density(end+1) = 2*n_edge/(n_nodes*(n_nodes-1));
        triangle_density(end+1) = 6*n_triangle/(n_nodes*(n_nodes-1)*(n_nodes-2));
        
            
        %extract largest connected component
        [x, W2, W3, T3] = MaxConnectedSubgraph(x, c2, c3, W2, W3, T3);    
        n_nodes = size(W2,2);

        
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
        
        
        %reorder nodes according to the embedding
        [~,idx_linear] = sort(x_est_linear);
        W2_reorder_linear = W2(idx_linear,idx_linear);
        W3_reorder_linear = W3(idx_linear,idx_linear);

        %reorder nodes according to the embedding
        [~,idx_periodic] = sort(x_est_periodic);
        W2_reorder_periodic = W2(idx_periodic,idx_periodic);
        W3_reorder_periodic = W3(idx_periodic,idx_periodic);


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
        lnP_linear = [];  lnP_linear_scaled = []; 
        lnP_periodic = []; 
        
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
        
        %maximum likelihood of gamma
        [max_lnP_linear(end+1), max_linear_idx] = max(lnP_linear);
        gamma_max_linear = gamma_array(max_linear_idx); 
        [max_lnP_linear_scaled(end+1), max_linear_scaled_idx] = max(lnP_linear_scaled);
        gamma_max_linear_scaled = gamma_array(max_linear_scaled_idx); 
        [max_lnP_periodic(end+1), max_periodic_idx] = max(lnP_periodic);
        gamma_max_periodic = gamma_array(max_periodic_idx); 



        if data_type == "cluster"
            
            cluster_input = kmeans(transpose(x), K);
            cluster_est_linear = kmeans(x_est_linear, K);
            cluster_est_periodic = kmeans(x_est_periodic, K);
            rand_linear(end+1) = CalculateRandIndex(cluster_input, cluster_est_linear);
            rand_periodic(end+1) = CalculateRandIndex(cluster_input, cluster_est_periodic);

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
        
    end
    
             %perform K-mean and compare with input




cla(gca,'reset')
yyaxis left
plt = plot(gamma_input, max_lnP_linear, '-b', 'LineWidth',1.5);
hold on;
plot(gamma_input, max_lnP_periodic, '-r', 'LineWidth',1.5);
yyaxis right
plot(gamma_input, edge_density, '--k', 'LineWidth',1.5);
legend({'Linear','Periodic', 'Edge Density'},'FontSize', 20,'Location','best');
xlabel('\gamma_{in}','FontSize', 13);
yyaxis left
ylabel('Max log-likelihood','FontSize', 13);
yyaxis right
ylabel('Edge density','FontSize', 13);
set(gca,'fontsize',25);
%set(gca,'YLim',[0 1.1])
plt.LineWidth = 2;
ax = gca;
exportgraphics(ax,strcat('plots/',input_shape,'_',data_type,'_maxlnP.eps'),'Resolution',300) 
hold off;

% plot rand_index
cla(gca,'reset')
yyaxis left
plt = plot(gamma_input, max_lnP_linear-max_lnP_periodic, '-b', 'LineWidth',1.5);
yline(0);
yyaxis right
plot(gamma_input, edge_density, '--k', 'LineWidth',1.5);
legend({'Linear-Periodic', 'Edge Density'},'FontSize', 20,'Location','best');
xlabel('\gamma_{in}','FontSize', 13);
yyaxis left
ylabel('Max Log-likelihood Diff','FontSize', 13);
yyaxis right
ylabel('Edge density','FontSize', 13);
set(gca,'fontsize',25);
%set(gca,'YLim',[0 1.1])
plt.LineWidth = 2;
ax = gca;
exportgraphics(ax,strcat('plots/',input_shape,'_',data_type,'_maxlnPDiff.eps'),'Resolution',300) 
hold off;

% plot rand_index
if data_type == "cluster"
    cla(gca,'reset')

    plt = plot(gamma_input, rand_linear, 'b', 'LineWidth',1.5);
    hold on;
    plot(gamma_input, rand_periodic, '-r', 'LineWidth',1.5);
    plot(gamma_input, edge_density, '--k', 'LineWidth',1.5);
    plot(gamma_input, triangle_density, '-*k', 'LineWidth',1.5);
    legend({'Linear','Periodic', 'Edge Density', 'Triangle Density'},'FontSize', 20,'Location','east');
    xlabel('\gamma','FontSize', 13);
    ylabel('Rand Index','FontSize', 13);
    set(gca,'fontsize',30);
    set(gca,'YLim',[0 1.1])
    plt.LineWidth = 2;
    ax = gca;
    exportgraphics(ax,strcat('plots/',input_shape,'_',data_type,'_rand.eps'),'Resolution',300) 
    hold off;
end
cla(gca,'reset')


end
end
toc
load handel
sound(y,Fs)