%%%% Plot Figure 1(e)-(f) and Figure 2(e)-(f) in the paper

clear all 
clc
close all
set(groot,'defaultFigureVisible','on') % 'on' to turn back on.  

gamma_input = 1; % gamma for generating graph
input_shape = "linear"; %% Plot Figure 1(e)-(f)
%input_shape = "periodic";  %% Plot Figure 2(e)-(f)
data_type = "cluster";


load(strcat(input_shape,'_',data_type, '_gamma=', num2str(gamma_input,2),'.mat'));
c3_inv_array = c3_inp_array;
%%%%%%%%%%%%%%%% PLOTS %%%%%%%%%%%%%%%%%

    for select_gamma_in = 1:length(gamma_input)
        gamma_plot = gamma_input(select_gamma_in);

        % plot heatmap of linear rand
        if data_type == "cluster"
        figure 
        h = heatmap(round(mean(rand_linear,3), 2));
        %caxis([0, 1]);
        %h.Title = strcat('Linear Rand, \gamma_{in}=', num2str(gamma_plot) );
        h.XLabel = 'c_3*';
        h.YLabel = 'c_3';
        h.XData = round(c3_inv_array,2);
        h.YData = round(c3_inp_array,2);
        ax = gca;% Requires R2020a or later
        set(gca,'fontsize', 15);
        h.YDisplayData = flipud(h.YDisplayData);  % equivalent to 'YDir', 'Reverse'
        exportgraphics(ax,strcat('plots/',input_shape,'_',data_type,'_gamma=', num2str(gamma_plot),'_rand_lin.eps'),'Resolution',300) 
        
        % plot heatmap of periodic rand
        h = heatmap(round(mean(rand_periodic,3), 2));
        %caxis([0, 1]);
        %h.Title = strcat('Periodic Rand, \gamma_{in}=', num2str(gamma_plot) );
        h.XLabel = 'c_3*';
        h.YLabel = 'c_3';
        h.XData = round(c3_inv_array,2);
        h.YData = round(c3_inp_array,2);
        set(gca,'fontsize', 15);
        h.YDisplayData = flipud(h.YDisplayData);  % equivalent to 'YDir', 'Reverse'
        ax = gca;% Requires R2020a or later
        exportgraphics(ax,strcat('plots/',input_shape,'_',data_type,'_gamma=', num2str(gamma_plot),'_rand_per.eps'),'Resolution',300) 
        
        
        end

        % plot heatmap of linear max log-likelihood
        h = heatmap(round(mean(max_lnP_linear, 3)/10^6, 3),'Colormap',parula);
        h.CellLabelFormat = '%.2f';
        %h.Title = strcat('Max LnP Linear, \gamma_{in}=', num2str(gamma_plot) );
        h.XLabel = 'c_3*';
        h.YLabel = 'c_3';
        h.XData = round(c3_inv_array,2);
        h.YData = round(c3_inp_array,2);
        ax = gca;% Requires R2020a or later
        set(gca,'fontsize', 15);
        h.YDisplayData = flipud(h.YDisplayData);  % equivalent to 'YDir', 'Reverse'
        exportgraphics(ax,strcat('plots/',input_shape,'_',data_type,'_gamma=', num2str(gamma_plot),'_maxlnP_lin.eps'),'Resolution',300) 
        
        % plot waterfall plot of linear max log-likelihood
        h = waterfall(round(mean(max_lnP_linear, 3)/10^6, 3));
        %h.Title = strcat('Max LnP Linear, \gamma_{in}=', num2str(gamma_plot) );
        %h.XLabel = 'c_3*';
        %h.YLabel = 'c_3';
        h.XData = round(c3_inv_array,2);
        h.YData = round(c3_inp_array,2);
        ax = gca;% Requires R2020a or later
        set(gca,'fontsize', 15);
        exportgraphics(ax,strcat('plots/waterfall_',input_shape,'_',data_type,'_gamma=', num2str(gamma_plot),'_maxlnP_lin.eps'),'Resolution',300) 
        
        
        % plot map of periodic max log-likelihood
        h = heatmap(round(mean(max_lnP_periodic,3)/10^6, 3),'Colormap',parula);
        h.CellLabelFormat = '%.2f';
        %h.Title = strcat('Max LnP Periodic, \gamma_{in}=', num2str(gamma_plot) );
        h.XLabel = 'c_3*';
        h.YLabel = 'c_3';
        h.XData = round(c3_inv_array,2);
        h.YData = round(c3_inp_array,2);
        ax = gca;% Requires R2020a or later
        set(gca,'fontsize', 15);
        h.YDisplayData = flipud(h.YDisplayData);  % equivalent to 'YDir', 'Reverse'
        exportgraphics(ax,strcat('plots/',input_shape,'_',data_type,'_gamma=', num2str(gamma_plot),'_maxlnP_per.eps'),'Resolution',300) 
        
        % plot waterfall plot of periodic max log-likelihood
        h = waterfall(round(mean(max_lnP_periodic,3)/10^6, 3));
        %h.Title = strcat('Max LnP Periodic, \gamma_{in}=', num2str(gamma_plot) );
        %h.CellLabelFormat = '%.2f';
        %h.XLabel = 'c_3*';
        %h.YLabel = 'c_3';
        h.XData = round(c3_inv_array,2);
        h.YData = round(c3_inp_array,2);
        ax = gca;% Requires R2020a or later
        set(gca,'fontsize', 15);
        exportgraphics(ax,strcat('plots/waterfall_',input_shape,'_',data_type,'_gamma=', num2str(gamma_plot),'_maxlnP_per.eps'),'Resolution',300) 
        
        
    end