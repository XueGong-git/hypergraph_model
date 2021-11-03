tic
clear


c2 = 1/2; % weight of simple edge
n = 327;
E2 = []; % edge list
W2 = zeros(n); % 2nd order adjacency matrix  

E3 = []; % edge list
W3 = zeros(n); % 2nd order adjacency matrix  


%read hyper edges in a cell array

fid = fopen('highschool.txt');
line1 = fgetl(fid);
res=line1;
while ischar(line1)
line1 = fgetl(fid);
res = char(res,line1);
end
fclose(fid);


for k=1:size(res,1)
  E{k} = str2num(res(k,:));
  %update edge list
  if length(E{k}) == 2
      E2(end+1,:) = [E{k}(1),E{k}(2)];
      E2(end+1,:) = [E{k}(2),E{k}(1)];
      W2(E{k}(1),E{k}(2))=1; % only fill uppder triangle
  elseif length(E{k}) == 3
      E3(end+1,:) = E{k};
      E3(end+1,:) = [E{k}(1),E{k}(2),E{k}(3)];
      E3(end+1,:) = [E{k}(1),E{k}(3),E{k}(2)];
      E3(end+1,:) = [E{k}(2),E{k}(1),E{k}(3)];
      E3(end+1,:) = [E{k}(2),E{k}(3),E{k}(1)];
      E3(end+1,:) = [E{k}(3),E{k}(1),E{k}(2)];
      E3(end+1,:) = [E{k}(3),E{k}(2),E{k}(1)];
        
      W2(E{k}(1),E{k}(2))=1; %only fill uppder triangle
      W3(E{k}(1),E{k}(2))= W3(E{k}(2),E{k}(2))+1;
      W3(E{k}(1),E{k}(3))= W3(E{k}(1),E{k}(3))+1;
      W3(E{k}(2),E{k}(3))= W3(E{k}(2),E{k}(3))+1;
  end
end

%symmetryze matrix
W2 = W2 + W2'; W3 = W3 + W3'; 

data_type = "highschool";

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


% plot L2 eigenvalues
% eigenvalues of L2 and L3
d2 = sum(W2, 2); D2 = diag(d2);
d3 = sum(W3, 2); D3 = diag(d3);
d2_inv = 1./d2; d2_inv(isinf(d2_inv)|isnan(d2_inv)) = 0;
d3_inv = 1./d3; d3_inv(isinf(d3_inv)|isnan(d3_inv)) = 0;
L2 = D2 - W2;
L3 = D3 - W3;
%normalize Laplacian
L2 = diag(d2_inv)*L2;
L3 = diag(d3_inv)*L3;
[V2, lambda2] = eigs(L2,n,'smallestabs'); % all eigenvectors ranging from smallest eigenvalue to largest eigenvalue
[V3, lambda3] = eigs(L3,n,'smallestabs'); % all eigenvectors ranging from smallest eigenvalue to largest eigenvalue


eigenvalues2 = diag(lambda2);
scatter(1:40, eigenvalues2(1:40), 100, 'MarkerFaceColor','black');
xlabel('i','FontSize', 13);
ylabel('Eigenvalue \lambda_i','FontSize', 13);
set(gca,'fontsize',30);
ax = gca;
exportgraphics(ax,strcat('plots/periodic_',data_type,'_L2_eigenvalues.eps'),'Resolution',300) 

% plot L3 eigenvalues
eigenvalues3 = diag(lambda3);
scatter(1:40, eigenvalues3(1:40), 100, 'MarkerFaceColor','black');
xlabel('i','FontSize', 13);
ylabel('Eigenvalue \lambda_i','FontSize', 13);
set(gca,'fontsize',30);
ax = gca;
exportgraphics(ax,strcat('plots/periodic_',data_type,'_L3_eigenvalues.eps'),'Resolution',300) 

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
exportgraphics(t,strcat('plots/periodic_',data_type,'_L2_eigenvectors.eps'),'Resolution',300) 


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
exportgraphics(t,strcat('plots/periodic_',data_type,'_L3_eigenvectors.eps'),'Resolution',300) 



%shuffle input adjacency matrix
idx_rand = randperm(n);% shuffle the nodes
W2 = W2(idx_rand,idx_rand);
W3 = W3(idx_rand,idx_rand);

for c3 = [0, 1/3, 2/3, 100/3] % weight of triangles

[x_est, V, lambda] = PeriodicHypergraphEmbedding(W2, W3, c2, c3, "norm");

%reorder nodes according to the embedding
[~,idx] = sort(x_est);
W2_reorder = W2(idx,idx);
W3_reorder = W3(idx,idx);


% Generate plots
% 2-d scatter plot for edges
figure
s = scatter(x_est(E2(:,1)),x_est(E2(:,2)),'MarkerFaceColor','black','MarkerEdgeColor','none');
xlim([0 n])
ylim([0 n])
alpha(s,0.3) % transparent color
set(gca,'fontsize',30);
ax = gca;
exportgraphics(ax,strcat('plots/periodic_',data_type,'_edges_c3=',num2str(round(c3,1)),'.eps'),'Resolution',300) 

% 3-d scatter plot for triangle
s = scatter3(E3(:,1),E3(:,2),E3(:,3),10,'MarkerFaceColor','black','MarkerEdgeColor','none');
xlim([0 n])
xlim([0 n])
zlim([0 n])
alpha(s,0.2) % transparent color
set(gca,'fontsize',30);
ax = gca;
exportgraphics(ax,strcat('plots/periodic_',data_type,'_triangles_c3=',num2str(round(c3,1)),'.eps'),'Resolution',300) 


% plot reordered W2
imagesc(W2_reorder,[0,1]); %plot color map of original matrix
colormap(flipud(gray(2)));
set(gca,'FontSize',30) ;
ax = gca;% Requires R2020a or later
exportgraphics(ax,strcat('plots/periodic_',data_type,'_W2_c3=',num2str(round(c3,1)),'.eps'),'Resolution',300) 


% plot reordered W3
imagesc(W3_reorder); % plot color map of original matrix
colormap(flipud(gray(256)));colorbar
set(gca,'FontSize',30) ;
set(gca,'ColorScale','log')
ax = gca;% Requires R2020a or later
exportgraphics(ax,strcat('plots/periodic_',data_type,'_W3_c3=',num2str(round(c3,1)),'.eps'),'Resolution',300) 


% plot eigenvalues of L
eigenvalues = diag(lambda);
scatter(1:40, eigenvalues(1:40), 100, 'MarkerFaceColor','black');
xlabel('i','FontSize', 13);
ylabel('Eigenvalue \lambda_i','FontSize', 13);
set(gca,'fontsize',30);
ax = gca;
exportgraphics(ax,strcat('plots/periodic_',data_type,'_eigenvalues_c3=',num2str(round(c3,1)),'.eps'),'Resolution',300) 

% plot top eigenvectors
t = tiledlayout(4,1);
ax1 = nexttile;
plot(V(:,1), 'Color', 'black');
ylabel('v_0')
ax2 = nexttile;
plot(V(:,2), 'Color', 'black');
ylabel('v_1')
ax3 = nexttile;
plot(V(:,3), 'Color', 'black');
ylabel('v_2')
ax4 = nexttile;
plot(V(:,4), 'Color', 'black');
ylabel('v_3')
% Link the axes
linkaxes([ax1,ax2,ax3,ax4],'x');
linkaxes([ax1,ax2,ax3,ax4],'y');
% Add shared title and axis labels
title(t,'Top Eigenvectors')
% Move plots closer together
xticklabels(ax1,{})
t.TileSpacing = 'compact';
axis([0 n -1 1])
exportgraphics(t,strcat('plots/periodic_',data_type,'_eigenvectors_c3=',num2str(round(c3,1)),'.eps'),'Resolution',300) 



end

toc
