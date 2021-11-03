% GenerateModelLikelihood  Using the given embedding, calculate the
% likelihood of unweighted hyper graph
%
% INPUTS
% 
% - x   Node embedding
% - (W2, T3)       Hypergraph -- with a matrix and a tensor
% - c3  Weight of triangles 
% - c2  Weight of edge
% - gamma   Decay parameter for the linear hypergraph model
% - structure Wether it is linear or periodic
%
% OUTPUTS
% - lnP  log-Likelihood of the grpah
function [lnP] = CalculateModelLikelihood(x, W2, T3, c3, gamma, structure)
n = length(x);
I = zeros(n);
lnP = 0;

if structure == "linear"
%calculate pairwise incoherence
    for i = 1:n-1 % smallest node index
        for j = i+1:n % largest node index
             %calculate incoherene of nodes
            I(i,j) = (x(i)-x(j))^2;        
        end
    end
elseif structure == "periodic"
    for i = 1:n-1 % smallest node index
        for j = i+1:n % largest node index
            %calculate incoherene of nodes
            I(i,j) = abs(exp(1i*x(i))-exp(1i*x(j)))^2;        
        end
    end
end

% calculate probability of edges
for i = 1:n-1
    for j = i+1:n
        if W2(i,j) == 1
            lnP = lnP - log(1+exp(gamma*I(i,j)));
        end
    end
end



% calculate probability of triangles
for i = 1:n-2 % smallest node index
    for j = i+1:n-1 % second smallest index
        for k = j+1:n % largest node index
            %calculate incoherene of nodes
            I_R = I(i,j)+ I(i,k) + I(j,k);
            if T3(i,j,k) == 1
            lnP = lnP - log(1+exp(gamma*c3*I_R));
            end
        end
    end
end


end