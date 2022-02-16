% GenerateModelLikelihood  Using the given embedding, calculate the
% likelihood of unweighted hyper graph
%
% INPUTS
% 
% - x   Node embedding
% - (W2, T3)       Hypergraph -- with a matrix and a tensor. !!Make sure x,
% W2, T3 are in the same order.
% - c3  Weight of triangles 
% - c2  Weight of edge
% - gamma   Decay parameter for the linear hypergraph model
% - structure Wether it is linear or periodic
%
% OUTPUTS
% - lnP  log-Likelihood of the grpah
% - lnP0 log-likelihood of null graph

function [lnP, eta, P_edge, P_triangles] = CalculateModelLikelihood(x, W2, T3, c2, c3, gamma, structure, no_triangle)

if ~exist('no_triangle','var')
 % third parameter does not exist, so default it to something
  no_triangle = 0;
end
n = length(x);
I = zeros(n);
lnP = 0;
eta = 0;
lnP0 = 0;
P_edge = zeros(n, n);
P_triangles = zeros(n,n,n);
I3 = zeros(n,n,n); %incoherence between tripples of nodes

%sort x, W2 and T3 for visualization
%[x, idx] = sort(x);
%W2 = W2(idx, idx);
%T3 = T3(idx, idx, idx);

if structure == "linear"
%calculate pairwise incoherence
    for i = 1:n-1 % smallest node index
        for j = i+1:n % largest node index
             %calculate incoherene of nodes
            I(i,j) = (x(i)-x(j))^2;        
            I(j,i) = (x(i)-x(j))^2;        

        end
    end
elseif structure == "periodic"
    for i = 1:n-1 % smallest node index
        for j = i+1:n % largest node index
            %calculate incoherene of nodes
            I(i,j) = abs(exp(1i*x(i))-exp(1i*x(j)))^2;        
            I(j,i) = abs(exp(1i*x(i))-exp(1i*x(j)))^2;        

        end
    end
end

% calculate probability of edges
for i = 1:n-1
    for j = i+1:n
        lnP0 = lnP0 + c2*gamma*I(i,j) - log(1+exp(c2*gamma*I(i,j)));
        P_edge(i,j) = 1/(1+exp(c2*gamma*I(i,j)));
        if W2(i,j) == 1
            lnP = lnP - log(1+exp(c2*gamma*I(i,j))); 
            eta = eta + c2*I(i,j);
        else
            lnP = lnP + c2*gamma*I(i,j) - log(1+exp(c2*gamma*I(i,j)));
        end
    end
end

P_edge = P_edge + P_edge';

% calculate probability of triangles
if ~no_triangle
for i = 1:n-2 % smallest node index
    for j = i+1:n-1 % second smallest index
        for k = j+1:n % largest node index
            %fprintf('%d\n',round(i))
            %fprintf('%d\n',round(j))
            %fprintf('%d\n',round(k))

            %calculate incoherene of nodes
            I_R = I(i,j)+ I(i,k) + I(j,k);
            I3(i,j,k) = I_R;
            I3(i,k,j) = I_R;
            I3(j,i,k) = I_R;
            I3(j,k,i) = I_R;
            I3(k,i,j) = I_R;
            I3(k,j,i) = I_R;
            
            P_triangle = 1/(1+exp(c3*gamma*I_R));
            P_triangles(i,j,k) = P_triangle; %tensor
            P_triangles(i,k,j) = P_triangle; %tensor
            P_triangles(j,i,k) = P_triangle; %tensor
            P_triangles(j,k,i) = P_triangle; %tensor
            P_triangles(k,i,j) = P_triangle; %tensor
            P_triangles(k,j,i) = P_triangle; %tensor
            lnP0 = lnP0 + c3*gamma*I_R - log(1+exp(c3*gamma*I_R)); % log-likelihood of null graph
            if T3(i,j,k) == 1
                lnP = lnP - log(1+exp(c3*gamma*I_R));
                eta = eta + c3*I_R;
            else
                lnP = lnP + c3*gamma*I_R - log(1+exp(c3*gamma*I_R));

            end
        end
    end
end
end

%{
% plot I(i,j)
imagesc(I); %plot color map of original matrix
%colormap(flipud(gray()));
set(gca,'FontSize',30) ;
ax = gca;% Requires R2020a or later
colormap parula
colorbar
exportgraphics(ax,strcat('plots/',structure,'_I2_gamma=', num2str(round(gamma,2)),'.eps'),'Resolution',300) 


% plot I3(i,j)
nonzero_idx = find(I3);
[px,py,pz] = ind2sub(size(I3),nonzero_idx);
scatter3(px,py,pz,30, nonzeros(I3),'filled');
set(gca,'XLim',[1 n],'YLim',[1 n],'ZLim',[1 n])
colormap parula
colorbar
ax = gca;% Requires R2020a or later
exportgraphics(ax,strcat('plots/',structure,'_I3_gamma=', num2str(round(gamma,2)),'.eps'),'Resolution',300) 

% plot probability of edges
imagesc(P_edge); %plot color map of original matrix
%colormap(flipud(gray()));
set(gca,'FontSize',30) ;
ax = gca;% Requires R2020a or later
colormap parula
colorbar
%set(gca,'ColorScale','log')
exportgraphics(ax,strcat('plots/',structure,'_P2_gamma=', num2str(round(gamma,2)),'.eps'),'Resolution',300) 


% plot probability of triangles
nonzero_idx = find(P_triangles);
[px,py,pz] = ind2sub(size(P_triangles),nonzero_idx);
scatter3(px,py,pz,30, nonzeros(P_triangles),'filled');
set(gca,'XLim',[1 n],'YLim',[1 n],'ZLim',[1 n])
set(gca,'FontSize',30) ;
colormap parula
colorbar
%set(gca,'ColorScale','log')
ax = gca;% Requires R2020a or later
exportgraphics(ax,strcat('plots/',structure,'_P3_gamma=', num2str(round(gamma,2)),'.eps'),'Resolution',300) 


%plot probability W2ijlnPij + (1-W2ij)ln(1-Pij)
PG2 = W2.*log(P_edge) + (1-W2).*log(1-P_edge);
PG2(isnan(PG2)) = 0; PG2(isinf(PG2)) = 0;
imagesc(PG2); %plot color map of original matrix
set(gca,'FontSize',30) ;
ax = gca;% Requires R2020a or later
colormap parula
colorbar
caxis(ax,[-3 0]);
exportgraphics(ax,strcat('plots/',structure,'_PG2_gamma=', num2str(round(gamma,2)),'.eps'),'Resolution',300) 


%plot triangles
T3(isnan(T3)) = 0; T3(isinf(T3)) = 0;
nonzero_idx = find(T3);
[px,py,pz] = ind2sub(size(T3),nonzero_idx);
scatter3(px,py,pz,30, nonzeros(T3),'filled');
set(gca,'XLim',[1 n],'YLim',[1 n],'ZLim',[1 n])
set(gca,'FontSize',30) ;
colormap parula
colorbar
ax = gca;% Requires R2020a or later
exportgraphics(ax,strcat('plots/',structure,'_T3_gamma=', num2str(round(gamma,2)),'.eps'),'Resolution',300) 


%plot probability T3ijklnPijk + (1-T3ijk)ln(1-Pij)
PG3 = T3.*log(P_triangles) + (1-T3).*log(1-P_triangles);
PG3(isnan(PG3)) = 0; PG3(isinf(PG3)) = 0;
nonzero_idx = find(PG3);
[px,py,pz] = ind2sub(size(PG3),nonzero_idx);
scatter3(px,py,pz,30, nonzeros(PG3),'filled');
set(gca,'XLim',[1 n],'YLim',[1 n],'ZLim',[1 n])
colormap parula
colorbar
ax = gca;% Requires R2020a or later
exportgraphics(ax,strcat('plots/',structure,'_PG3_gamma=', num2str(round(gamma,2)),'.eps'),'Resolution',300) 

checklnG = sum(PG2, 'all') + sum(PG3,'all');
%}

end