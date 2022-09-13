tic
clear
set(groot,'defaultFigureVisible','off') % 'on' to turn back on.  

n = 6;
K = 5; % number of clusters 
m = n/K; % number of nodes per cluster
%a = 0.1*2*pi/K; % noise;


c2 = 1; % weight of simple edge
c3 = 1/3; % weight of triangles
%gamma_array = gamma_input;
gamma_array = 0:0.1:5; % gamma for likelihood plot
gamma_input = 0.5; % gamma for generating graph
no_triangle = 0; 
for input_shape = "linear"
    
for input = 1
    rand_linear = [];
    rand_periodic = [];
    triangle_density = [];
    edge_density = [];
    max_lnP = [];
    for gamma = gamma_input
        
        if input_shape == "linear" 
            a = 0.2;

            switch input

            case 1 
                x = linspace(1,n,n); % uniform from 0 to 2pi
                data_type = "uniform";
            case 2 
                x = sort(repmat(linspace(0,2, K),1,m)+(2*a*rand(1,n)-a));
                %x = sort(repmat(linspace(-pi,pi,K),1,m)+(2*a*rand(1,n)-a)); % angles from -pi to pi
                data_type = "cluster";
            case 3
                x = linspace(0,0,n); %overlapping x
                data_type = "overlap";
            end
            [W2, W3, T3] = GenerateLinearHypergraph(x, gamma, c2, c3, data_type, no_triangle);

        elseif input_shape == "periodic" 
            a = 0.1*pi; % noise;

            switch input
    
            case 1 
                x = linspace(2*pi/n,2*pi,n);%uniform from 0 to 2pi
                data_type = "uniform";
            case 2 
                x = sort(repmat(linspace(-pi,pi,K),1,m)+(2*a*rand(1,n)-a)); % angles from -pi to pi
                data_type = "cluster";
            case 3
                x = linspace(0,0,n); %overlapping x
                data_type = "overlap";
            end
            [W2, W3, T3] = GeneratePeriodicHypergraph(x, gamma, c2, c3, data_type);

        end
        
        
        %calculate likelihood
        lnP = [];
        
        %calculate eta
        [~, eta] = CalculateModelLikelihood(x, W2, T3, c2, c3, 2, input_shape, no_triangle);
    
        
        for test_gamma = gamma_array
            [lnP(end+1), test_eta] = CalculateModelLikelihood(x, W2, T3, c2, c3, test_gamma, input_shape, no_triangle);
        end
        
        %maximum likelihood of gamma
        [max_lnP(end+1), max_idx] = max(lnP);
        gamma_max = gamma_array(max_idx);


       
        
        % plot likelihood
        plt = plot(gamma_array, lnP, 'b', 'LineWidth',1.5);
        hold on;
        plot(gamma_array, -1*eta*gamma_array, ':b', 'LineWidth',1);
        plot(gamma_max, lnP(max_idx), 'ob', 'MarkerSize',10, 'LineWidth',2);
        xline(gamma,'-',{'True \gamma'},'fontsize',20)
        legend({'Log-likelihood','Slope = -Incoherence','MLE'},'FontSize', 20,'Location','best');
        xlabel('\gamma','FontSize', 13);
        ylabel('Log-likelihood','FontSize', 13);
        set(gca,'fontsize',30);
        set(gca,'XLim',[0 max(gamma_array)])
        plt.LineWidth = 2;
        ax = gca;
        exportgraphics(ax,strcat('plots/',input_shape,'_',data_type, '_model_comparison_gamma_', num2str(round(gamma,2)),'.eps'),'Resolution',300) 
        hold off;
    
        
    end
    
end
cla(gca,'reset')


end

toc
load handel
sound(y,Fs)