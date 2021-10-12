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

clear
%input parameters
n = 25;
x = linspace(0,0,n); %overlapping x
%x = linspace(1,n,n);
K = 5; a= 1; m = n/K;
%x = sort(repmat(linspace(n/K,n,K),1,m)+(2*a*rand(1,n)-a)); % trophic levels

gamma = 0.5;
I = zeros(n);
f = zeros(n);

%calculate pairwise incoherence
for i = 1:n-1 % smallest node index
    for j = i+1:n % largest node index
        %calculate incoherene of nodes
        I(i,j) = (x(i)-x(j))^2;        
    end
end

edgeList = [];
%generate simple edges
c2 = 1/2; % weight of simple edge
for i = 1:n-1 % smallest node index
    for j = i+1:n % largest node index
        %calculate likelihood
        f(i,j) = 1/(1+exp(gamma*c2*I(i,j)));
          if rand() < f(i,j)
                %edgeList(end+1,:) = [x(i),x(j)];
                %edgeList(end+1,:) = [x(j),x(i)];
                edgeList(end+1,:) = [i,j];
                edgeList(end+1,:) = [j,i];
            end
    end
end

%generate hyperedges
randM = rand(n,n);
A2 = randM < f; A2 = A2 + A2';

triangleList = [];
%generate  hyperedges that connect 3 nodes
c3 = 1/3; % weight of simple edge
for i = 1:n-2 % smallest node index
    for j = i+1:n-1 % second smallest index
        for k = j+1:n % largest node index
            %calculate incoherene of nodes
            I_R = I(i,j)+ I(i,k) + I(j,k);
            f(i,j,k) = 1/(1+exp(gamma*c3*I_R));
            if rand() < f(i,j,k)
                %triangleList(end+1,:) = [x(i),x(j),x(k)];
                %triangleList(end+1,:) = [x(j),x(k),x(i)];
                %triangleList(end+1,:) = [x(k),x(i),x(j)];
                %triangleList(end+1,:) = [x(k),x(j),x(i)];
                %triangleList(end+1,:) = [x(j),x(i),x(k)];
                %triangleList(end+1,:) = [x(i),x(k),x(j)];
                
                triangleList(end+1,:) = [i,j,k];
                triangleList(end+1,:) = [j,k,i];
                triangleList(end+1,:) = [k,i,j];
                triangleList(end+1,:) = [k,j,i];
                triangleList(end+1,:) = [j,i,k];
                triangleList(end+1,:) = [i,k,j];
            end
        end
    end
end

%generate hyperedges
randM = rand(n,n,n);
A3 = randM < f; %only upper triangles are filled


%2-d scatter plot for edges
scatter(edgeList(:,1),edgeList(:,2),'black','filled')
ax = gca;
exportgraphics(ax,strcat('plots/linear_hypergraph_edges_overlap_gamma=', num2str(round(gamma,2)),'.eps'),'Resolution',300) 

%3-d scatter plot for triangle
scatter3(triangleList(:,1),triangleList(:,2),triangleList(:,3),'black','filled')
ax = gca;
exportgraphics(ax,strcat('plots/linear_hypergraph_triangle_overlap_gamma=', num2str(round(gamma,2)),'.eps'),'Resolution',300) 

