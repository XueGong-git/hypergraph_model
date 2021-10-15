% generate_linear_hypergraph  Generate unweighted hypergraph that shows a linear
% structure
%
% INPUTS
% 
% - n   Number of nodes
% - x   Node embedding
% - gamma   Decay parameter for the linear hypergraph model
%
% OUTPUTS
% - (A2, A3)       Hypergraph -- with a matrix + a tensor
tic
clear


c2 = 1/2; % weight of simple edge
c3 = 1/3; % weight of simple edge
n = 327
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
res =char(res,line1);
end
fclose(fid);

n = 327
E2 = []; % edge list
W2 = zeros(n); % 2nd order adjacency matrix  

E3 = []; % edge list
W3 = zeros(n); % 2nd order adjacency matrix  


for k=1:size(res,1)
  E{k}=str2num(res(k,:));
  %update edge list
  if length(E{k}) == 2
      E2(end+1,:) = [E{k}(1),E{k}(2)];
      E2(end+1,:) = [E{k}(2),E{k}(1)];
      W2(E{k}(1),E{k}(2))=1; %only fill uppder triangle
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

c2 = 0;

%shuffle input adjacency matrix
idx_rand = randperm(size(W2,1));% shuffle the nodes
[x_est] = LinearHypergraphEmbedding(W2, W3, c2, c3);


%reorder nodes
[~,idx] = sort(x_est);
W2_reorder = W2(idx,idx);
W3_reorder = W3(idx,idx);


% 2-d scatter plot for edges
s = scatter(E2(:,1),E2(:,2),'MarkerFaceColor','black','MarkerEdgeColor','none')
alpha(s,0.3) % transparent color
set(gca,'fontsize',30);
ax = gca;
exportgraphics(ax,strcat('plots/',data_type,'_edges_c3=',num2str(round(c3,1)),'.eps'),'Resolution',300) 

% 3-d scatter plot for triangle
ax1 = nexttile;
s=scatter3(ax1,E3(:,1),E3(:,2),E3(:,3),10,'MarkerFaceColor','black','MarkerEdgeColor','none')
alpha(s,0.2) % transparent color
set(gca,'fontsize',30);
ax = gca;
exportgraphics(ax,strcat('plots/',data_type,'_triangles_c3=',num2str(round(c3,1)),'.eps'),'Resolution',300) 


%plot reordered W2
imagesc(W2_reorder,[0,1]); %plot color map of original matrix
colormap(flipud(gray(2)));
set(gca,'FontSize',30) ;
ax = gca;% Requires R2020a or later
exportgraphics(ax,strcat('plots/',data_type,'_W2_c3=',num2str(round(c3,1)),'.eps'),'Resolution',300) 


% plot 3rd order adjacency matirx
imagesc(W3_reorder); % plot color map of original matrix
colormap(flipud(gray(256)));colorbar
set(gca,'FontSize',30) ;
ax = gca;% Requires R2020a or later
exportgraphics(ax,strcat('plots/',data_type,'_W3_c3=',num2str(round(c3,1)),'.eps'),'Resolution',300) 


toc
