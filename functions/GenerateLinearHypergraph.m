% GenerateLinearHypergraph  Generate unweighted hypergraph that shows a linear
% structure
%
% INPUTS
% 
% - x   Node embedding
% - c2  Weight of edge
% - c3  Weight of triangles 
% - gamma   Decay parameter for the linear hypergraph model
%
% OUTPUTS
% - (A2, A3)       Hypergraph -- with a matrix + a tensor

function [E2, W2, W3, T3] = GenerateLinearHypergraph(x, gamma, c2, c3, no_triangle, data_name)

if ~exist('no_triangle','var')
 % third parameter does not exist, so default it to 0
  no_triangle = 0;
end
n = length(x);

size_lcc = 0;

while size_lcc < n
    

    I = zeros(n);
    f = zeros(n);

    E2 = []; % edge list
    W2 = zeros(n); % 2nd order adjacency matrix  

    E3 = []; % edge list
    W3 = zeros(n); % 2nd order adjacency matrix  
    T3 = zeros(n,n,n);

    times = [];
    times_idx = 0;
    nverts = [];
    simplices = [];



    %calculate pairwise incoherence
    for i = 1:n-1 % smallest node index
        for j = i+1:n % largest node index
            %calculate incoherene of nodes
            I(i,j) = (x(i)-x(j))^2;        
        end
    end


    %generate simple edges
    for i = 1:n-1 % smallest node index
        for j = i+1:n % largest node index
            %calculate likelihood
            f(i,j) = 1/(1+exp(gamma*c2*I(i,j)));
              if rand() < f(i,j)
                    %E2(end+1,:) = [x(i),x(j)];
                    %E2(end+1,:) = [x(j),x(i)];
                    W2(i,j)=1; %only fill uppder triangle

                    E2(end+1,:) = [i,j];
                    E2(end+1,:) = [j,i];

                    % update the simplex list
                    times_idx = times_idx + 1;
                    nverts(end + 1) = 2;
                    simplices(end + 1) = i;
                    simplices(end + 1) = j;

             end
        end
    end
    W2 = W2 + W2'; % get a symmetric matrix

    if no_triangle ~= 0
    %generate  hyperedges that connect 3 nodes
    for i = 1:n-2 % smallest node index
        for j = i+1:n-1 % second smallest index
            for k = j+1:n % largest node index
                %calculate incoherene of nodes
                I_R = I(i,j)+ I(i,k) + I(j,k);
                f(i,j,k) = 1/(1+exp(gamma*c3*I_R));
                if rand() < f(i,j,k)
                    T3(i,j,k) = 1; %tensor
                    T3(i,k,j) = 1; %tensor
                    T3(j,i,k) = 1; %tensor
                    T3(j,k,i) = 1; %tensor
                    T3(k,i,j) = 1; %tensor
                    T3(k,j,i) = 1; %tensor

                    % 3rd order adjacency matrix, only fill the upper triangle
                    W3(i,j)= W3(i,j)+1;
                    W3(i,k)= W3(i,k)+1;
                    W3(j,k)= W3(j,k)+1;


                    %E3(end+1,:) = [i,j,k];
                    %E3(end+1,:) = [j,k,i];
                    %E3(end+1,:) = [k,i,j];
                    %E3(end+1,:) = [k,j,i];
                    %E3(end+1,:) = [j,i,k];
                    %E3(end+1,:) = [i,k,j];

                    % update the simplex list
                    times_idx = times_idx + 1;
                    nverts(end + 1) = 3;
                    simplices(end + 1) = i;
                    simplices(end + 1) = j;
                    simplices(end + 1) = k;


                end
            end
        end
    end
    end
    %get a symmetric adjacency matirx
    W3 = W3 + W3';
    
    % check if the hypergraph is connected
    [x, ~] = MaxConnectedSubgraph(x, c2, c3, W2, W3, T3);   
    size_lcc = size(x, 2);
end

% output lists of hyperedges

%{
%%%%%%%%%%%%%%%%%% PLOTS %%%%%%%%%%%%%%%%%%%%%
% 2-d scatter plot for edges
if ~isempty(E2)
    s=scatter(E2(:,1),E2(:,2),'MarkerFaceColor','black','MarkerEdgeColor','none');
    alpha(s,0.3) % transparent color
    set(gca,'fontsize',30);
    ax = gca;
    exportgraphics(ax,strcat('plots/linear_hygraph_edges_', data_type,'_gamma=', num2str(round(gamma,2)),'.eps'),'Resolution',300) 
end

% 3-d scatter plot for triangle
if ~isempty(E3)
    s=scatter3(E3(:,1),E3(:,2),E3(:,3),'MarkerFaceColor','black','MarkerEdgeColor','none');
    alpha(s,0.3) % transparent color
    set(gca,'fontsize',30);
    ax = gca;
    exportgraphics(ax,strcat('plots/linear_hygraph_triangles_', data_type,'_gamma=', num2str(round(gamma,2)),'.eps'),'Resolution',300) 
end
% plot 2nd order adjacency matrix
imagesc(W2,[0,1]); %plot color map of original matrix
colormap(flipud(gray(2)));
set(gca,'FontSize',30) ;
ax = gca;% Requires R2020a or later
exportgraphics(ax,strcat('plots/linear_hygraph_W2_', data_type,'_gamma=', num2str(round(gamma,2)),'.eps'),'Resolution',300) 


% plot 3rd order adjacency matirx
imagesc(W3); % plot color map of original matrix
colormap(flipud(gray(256)));colorbar
set(gca,'FontSize',30) ;
ax = gca;% Requires R2020a or later
exportgraphics(ax,strcat('plots/linear_hygraph_W3_', data_type,'_gamma=', num2str(round(gamma,2)),'.eps'),'Resolution',300) 

%%%%%%%%%%%%%%%%%%%%%%%%% END OF PLOTS %%%%%%%%%%%%%%%%%%%%%%%%
%}
if exist('data_name','var')
    %time
    times = randperm(times_idx); % randomly assign a timestamp
    fid = fopen(strcat('raw_data/', data_name,'/', data_name, '-times.txt'),'w+');
    fprintf(fid,'%d\n',times');
    fclose(fid);
    
    %nverts
    fid = fopen(strcat('raw_data/', data_name,'/', data_name, '-nverts.txt'),'w+');
    fprintf(fid,'%d\n',nverts');
    fclose(fid);
    
    %nverts
    fid = fopen(strcat('raw_data/', data_name,'/', data_name, '-simplices.txt'),'w+');
    fprintf(fid,'%d\n',simplices');
    fclose(fid);
end

end