
%n = 250;
%K = 5; % number of clusters 
%ntrial = 20;

%c2 = 1; % weight of simple edge
%gamma_array = 1; % gamma for likelihood plot
gamma_input = 1; % gamma for generating graph
c3_inp_array = 0:0.1:1; % weight of triangles for input graph
c3_inv_array = 0:0.1:1; % weight of triangles for inverse problem
input_shape = "linear";
data_type = "cluster";
load(strcat(input_shape,'_',data_type, '_gamma=', num2str(gamma_input,2),'.mat'), 'rand_linear', 'rand_periodic',  'max_lnP_linear', 'max_lnP_periodic', 'edge_density', 'triangle_density', 'correlation_lin', 'correlation_per');

%%%%%%%%%%%%%%%% PLOTS %%%%%%%%%%%%%%%%%

    for select_gamma_in = 1:length(gamma_input)
        gamma_plot = gamma_input(select_gamma_in);

        % plot heatmap of linear rand
        if data_type == "cluster"
        figure 
        h = heatmap(round(rand_linear(:,:,select_gamma_in), 2));
        %caxis([0, 1]);
        %h.Title = strcat('Linear Rand, \gamma_{in}=', num2str(gamma_plot) );
        h.XLabel = 'c3_{inv}';
        h.YLabel = 'c3_{inp}';
        h.XData = round(c3_inv_array,2);
        h.YData = round(c3_inp_array,2);
        ax = gca;% Requires R2020a or later
        exportgraphics(ax,strcat('plots/',input_shape,'_',data_type,'_gamma=', num2str(gamma_plot),'_rand_lin.eps'),'Resolution',300) 
        
        % plot heatmap of periodic rand
        h = heatmap(round(rand_periodic(:,:,select_gamma_in), 2));
        %caxis([0, 1]);
        %h.Title = strcat('Periodic Rand, \gamma_{in}=', num2str(gamma_plot) );
        h.XLabel = 'c3_{inv}';
        h.YLabel = 'c3_{inp}';
        h.XData = round(c3_inv_array,2);
        h.YData = round(c3_inp_array,2);
        ax = gca;% Requires R2020a or later
        exportgraphics(ax,strcat('plots/',input_shape,'_',data_type,'_gamma=', num2str(gamma_plot),'_rand_per.eps'),'Resolution',300) 
        
        % plot heatmap of linear - periodic rand
        h = heatmap(round(rand_linear(:,:,select_gamma_in)-rand_periodic(:,:,select_gamma_in), 2));
        %caxis([-0.5, 0.5]);
        %h.Title = strcat('Linear Rand - Periodic Rand, \gamma_{in}=', num2str(gamma_plot) );
        h.XLabel = 'c3_{inv}';
        h.YLabel = 'c3_{inp}';
        h.XData = round(c3_inv_array,2);
        h.YData = round(c3_inp_array,2);
        ax = gca;% Requires R2020a or later
        exportgraphics(ax,strcat('plots/',input_shape,'_',data_type,'_gamma=', num2str(gamma_plot),'_rand_lin_vs_per.eps'),'Resolution',300) 
        
        end

        % plot heatmap of linear max log-likelihood
        h = heatmap(round(max_lnP_linear(:,:,select_gamma_in), 2),'Colormap',parula,'CellLabelColor','none');
        h.CellLabelFormat = '%.1e';
        %h.Title = strcat('Max LnP Linear, \gamma_{in}=', num2str(gamma_plot) );
        h.XLabel = 'c3_{inv}';
        h.YLabel = 'c3_{inp}';
        h.XData = round(c3_inv_array,2);
        h.YData = round(c3_inp_array,2);
        ax = gca;% Requires R2020a or later
        exportgraphics(ax,strcat('plots/',input_shape,'_',data_type,'_gamma=', num2str(gamma_plot),'_maxlnP_lin.eps'),'Resolution',300) 
        
        % plot heatmap of periodic max log-likelihood
        h = heatmap(round(max_lnP_periodic(:,:,select_gamma_in), 2),'Colormap',parula,'CellLabelColor','none');
        %h.Title = strcat('Max LnP Periodic, \gamma_{in}=', num2str(gamma_plot) );
        h.CellLabelFormat = '%.1e';
        h.XLabel = 'c3_{inv}';
        h.YLabel = 'c3_{inp}';
        h.XData = round(c3_inv_array,2);
        h.YData = round(c3_inp_array,2);
        ax = gca;% Requires R2020a or later
        exportgraphics(ax,strcat('plots/',input_shape,'_',data_type,'_gamma=', num2str(gamma_plot),'_maxlnP_per.eps'),'Resolution',300) 
        
        % plot heatmap of linear - periodic log-likelihood
        h = heatmap(round(max_lnP_linear(:,:,select_gamma_in)-max_lnP_periodic(:,:,select_gamma_in), 2),'Colormap',parula,'CellLabelColor','none');
        %h.Title = strcat('Max LnP Linear-Periodic, \gamma_{in}=', num2str(gamma_plot) );
        h.XLabel = 'c3_{inv}';
        h.YLabel = 'c3_{inp}';
        h.XData = round(c3_inv_array,2);
        h.YData = round(c3_inp_array,2);
        ax = gca;% Requires R2020a or later
        exportgraphics(ax,strcat('plots/',input_shape,'_',data_type,'_gamma=', num2str(gamma_plot),'_maxlnP_lin_vs_per.eps'),'Resolution',300) 
      
        % plot heatmap of linear correlation
        h = heatmap(round(correlation_lin(:,:,select_gamma_in), 2),'Colormap',parula);
        %h.Title = strcat('Correlation betweeen estimation and truth (linear) \gamma_{in}=', num2str(gamma_plot) );
        h.XLabel = 'c3_{inv}';
        h.YLabel = 'c3_{inp}';
        h.XData = round(c3_inv_array,2);
        h.YData = round(c3_inp_array,2);
        ax = gca;% Requires R2020a or later
        exportgraphics(ax,strcat('plots/',input_shape,'_',data_type,'_gamma=', num2str(gamma_plot),'_corr_lin.eps'),'Resolution',300) 
        
        % plot heatmap of periodic correlation
        h = heatmap(round(correlation_per(:,:,select_gamma_in), 2),'Colormap',parula);
        %h.Title = strcat('Correlation betweeen estimation and truth (periodic) \gamma_{in}=', num2str(gamma_plot) );
        h.XLabel = 'c3_{inv}';
        h.YLabel = 'c3_{inp}';
        h.XData = round(c3_inv_array,2);
        h.YData = round(c3_inp_array,2);
        ax = gca;% Requires R2020a or later
        exportgraphics(ax,strcat('plots/',input_shape,'_',data_type,'_gamma=', num2str(gamma_plot),'_corr_per.eps'),'Resolution',300) 
     
        % plot heatmap of linear - periodic correlation
        h = heatmap(round(correlation_lin(:,:,select_gamma_in)-correlation_per(:,:,select_gamma_in), 2),'Colormap',parula);
        %h.Title = strcat('Difference between correlations \gamma_{in}=', num2str(gamma_plot) );
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