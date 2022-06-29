% Calculate generalized means
%
% INPUTS
% 
% 
%
% OUTPUTS
% - score for edges and triangles
%


function S_triangle = CalculateHarmMean(W3)

n = size(W3, 1);

S_triangle = zeros(n,n,n);


% calculate score of triangles
for i = 1:n-2 % smallest node index
    for j = i+1:n-1 % second smallest index
        for k = j+1:n % largest node index            
            mean_triangle = harmmean([ W3(i,j), W3(i,k), W3(j,k)]);
            S_triangle(i,j,k) = mean_triangle; %tensor
            S_triangle(i,k,j) = mean_triangle; %tensor
            S_triangle(j,i,k) = mean_triangle; %tensor
            S_triangle(j,k,i) = mean_triangle; %tensor
            S_triangle(k,i,j) = mean_triangle; %tensor
            S_triangle(k,j,i) = mean_triangle; %tensor
        end
    end
end
S_triangle = S_triangle(:);

end
