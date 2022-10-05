%%%%%% Plot Figure 1(c)-(d) and 2(c)-(d) in the paper

clear

c3 = 1/3; % weight of triangles
input_shape = "linear"; % Plot Figure 1(c)-(d)
%input_shape = "periodic"; % Plot Figure 2(c)-(d)
data_type = "cluster";

%load data
load(strcat(input_shape,'_',data_type, '_c3=', num2str(c3,2) ,'.mat'), ...
    'c2','c3','rand_linear', 'rand_periodic',  'max_lnP_linear', 'max_linear_idx', ...
'max_lnP_linear_scaled', 'max_linear_scaled_idx', 'max_lnP_periodic', 'max_periodic_idx', ...
'n_edge', 'n_nodes', 'n_triangle', 'gamma_array', 'gamma_input');


disp('************************')
disp('max_lnP_linear')
disp(['mean = ', num2str( mean(max_lnP_linear))])
disp(['std = ', num2str( std(max_lnP_linear))])
disp(['coefficient of variation = ', num2str(100*std(max_lnP_linear)./mean(max_lnP_linear),2), '%'])

disp('************************')
disp('max_lnP_periodic')
disp(['mean = ', num2str( mean(max_lnP_periodic))])
disp(['std = ', num2str( std(max_lnP_periodic))])
disp(['coefficient of variation = ', num2str(100*std(max_lnP_periodic)./mean(max_lnP_periodic),2), '%'])

disp('************************')
disp('rand_linear')
disp(['mean = ', num2str( mean(rand_linear))])
disp(['std = ', num2str( std(rand_linear))])
disp(['coefficient of variation = ', num2str(100*std(rand_linear)./mean(rand_linear),2), '%'])

disp('************************')
disp('rand_periodic')
disp(['mean = ', num2str( mean(rand_periodic))])
disp(['std = ', num2str( std(rand_periodic))])
disp(['coefficient of variation = ', num2str(100*std(rand_periodic)./mean(rand_periodic),2), '%'])



%%%%%%%%%%% plots %%%%%%%
%%%%%% likelihood plot %%%%%%

gamma_max_linear = gamma_array(max_linear_idx); 
gamma_max_linear_scaled = gamma_array(max_linear_scaled_idx); 
gamma_max_periodic = gamma_array(max_periodic_idx); 

%%%%%  RAND plot %%%%%%
% edge and triangle density
edge_density =  mean(2*n_edge./(n_nodes.*(n_nodes-1)),1);
triangle_density = mean(6*n_triangle./(n_nodes.*(n_nodes-1).*(n_nodes-2)),1);


if data_type == "cluster"
    % Rand Index
    cla(gca,'reset')
    yyaxis left
    % Linear ARI
    plt = plot(gamma_input, mean(rand_linear,1), '-*k', 'LineWidth',1.5);
    hold on
    rand_linear_s = sort(rand_linear);
    %CI
    xconf = [gamma_input gamma_input(end:-1:1)] ;         
    yconf = [rand_linear_s(5,:) rand_linear_s(size(rand_linear_s,1)-4, end:-1:1)];
    fill(xconf,yconf,'black','FaceAlpha', 0.3, 'EdgeColor', 'none');
    
    % Periodic ARI
    left_color = [0 0 0];
    right_color = [0 0 0];
    set(plt,'defaultAxesColorOrder',[left_color; right_color]);
    hold on;
    plot(gamma_input, mean(rand_periodic,1), '-or', 'LineWidth',1.5);
    rand_periodic_s = sort(rand_periodic);
    % CI
    xconf = [gamma_input gamma_input(end:-1:1)] ;         
    yconf = [rand_periodic_s(5,:) rand_periodic_s(size(rand_periodic_s,1)-4, end:-1:1)];
    fill(xconf,yconf,'red','FaceAlpha', 0.3, 'EdgeColor', 'none');
    
    % Edge density
    plot(gamma_input, edge_density, '--*k', 'LineWidth',1.5, 'Color', [0 0.4470 0.7410]);
    yyaxis right
    % Triangle density
    plot(gamma_input, triangle_density, '--ok', 'LineWidth',1.5, 'Color', [0.8500 0.3250 0.0980]);
    legend({'Linear','80% CI', 'Periodic','80% CI', 'Edge Density', 'Triangle Density'},'FontSize', 20,'Location','east');
    xlabel('\gamma_0','FontSize', 13);
    yyaxis left
    ylabel('ARI','FontSize', 13);
    yyaxis right
    ylabel('Triangle Density')
    set(gca,'fontsize',30);
    %set(gca,'YLim',[0 1.1])
    set(gca,'XLim',[0 max(gamma_input)])
    plt.LineWidth = 2;
    pos = get(gca, 'OuterPosition');
    set(gca,'OuterPosition',[pos(1) pos(2) pos(3) pos(4)]);
    ax = gca;
    exportgraphics(ax,strcat('plots/',input_shape,'_',data_type,'_rand.eps'),'Resolution',300) 
    hold off;
end
cla(gca,'reset')

%%% maximum likelihood
%%% plot confidence bounds
x = gamma_input;
y_array = max_lnP_linear;
ys = sort(y_array);
y = mean(y_array,1);


xconf = [x x(end:-1:1)] ;         
yconf = [ys(5,:) ys(size(y_array,1)-4, end:-1:1)];

figure
plot(x,y,'-*k','LineWidth',0.8, 'Color', [0 0 0])
hold on

p = fill(xconf,yconf,'black');
p.FaceColor = 	[0 0 0];      
p.EdgeColor = 'none';           
p.FaceAlpha = 0.3;



y_array = max_lnP_periodic;
ys = sort(y_array);
y = mean(y_array,1);

xconf = [x x(end:-1:1)] ;         
yconf = [ys(5,:) ys(size(rand_linear,1)-4, end:-1:1)];

plot(x,y,'--or','LineWidth',0.8, 'Color', [1 0 0])

p = fill(xconf,yconf,'red');
p.FaceColor = [1 0.8 0.8];      
p.EdgeColor = 'none';
p.FaceAlpha = 0.8;

xlabel('\gamma_0','FontSize', 13);
ylabel('Maximum LnP','FontSize', 13);
legend({'Linear','80% CI','Periodic','80% CI'},'FontSize', 20,'Location','best');
set(gca,'fontsize',30);
pos = get(gca, 'OuterPosition');
set(gca,'OuterPosition',[pos(1) pos(2) pos(3) pos(4)]);
%set(gca,'OuterPosition',[pos(1) pos(2)+0.05 pos(3) pos(4)-0.05]);
ax = gca;
exportgraphics(ax,strcat('plots/',input_shape,'_',data_type,'_lnP.eps'),'Resolution',300) 


hold off


