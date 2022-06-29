function [W2, W3, T3, E2, E3] = LoadSimplex(simplex_list, n)


E2 = []; % edge list
W2 = zeros(n); % 2nd order adjacency matrix  

E3 = []; % edge list
W3 = zeros(n); % 2nd order adjacency matrix  
T3 = zeros(n,n,n);


for k=1:size(simplex_list,2)
  E{k} = simplex_list{k};
  %update edge list
  if length(E{k}) == 2
      E2(end+1,:) = [E{k}(1),E{k}(2)];
      E2(end+1,:) = [E{k}(2),E{k}(1)];
      W2(E{k}(1),E{k}(2))=1; % only fill uppder triangle
  elseif length(E{k}) == 3
      E3(end+1,:) = [E{k}(1),E{k}(2),E{k}(3)];
      E3(end+1,:) = [E{k}(1),E{k}(3),E{k}(2)];
      E3(end+1,:) = [E{k}(2),E{k}(1),E{k}(3)];
      E3(end+1,:) = [E{k}(2),E{k}(3),E{k}(1)];
      E3(end+1,:) = [E{k}(3),E{k}(1),E{k}(2)];
      E3(end+1,:) = [E{k}(3),E{k}(2),E{k}(1)];
      
      T3(E{k}(1),E{k}(2),E{k}(3)) = 1; %tensor
      T3(E{k}(1),E{k}(3),E{k}(2)) = 1; %tensor
      T3(E{k}(2),E{k}(1),E{k}(3)) = 1; %tensor
      T3(E{k}(2),E{k}(3),E{k}(1)) = 1; %tensor
      T3(E{k}(3),E{k}(1),E{k}(2)) = 1; %tensor
      T3(E{k}(3),E{k}(2),E{k}(1)) = 1; %tensor
      
      
      W3(E{k}(1),E{k}(2))= W3(E{k}(2),E{k}(2))+1;
      W3(E{k}(1),E{k}(3))= W3(E{k}(1),E{k}(3))+1;
      W3(E{k}(2),E{k}(3))= W3(E{k}(2),E{k}(3))+1;
  end
end

%symmetryze matrix
W2 = W2 + W2'; W3 = W3 + W3'; 


end