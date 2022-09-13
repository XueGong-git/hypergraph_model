function [W2, W3, T3] = Edge2Matrix(n, E2, E3)


W2 = zeros(n); % 2nd order adjacency matrix  
W3 = zeros(n); % 2nd order adjacency matrix  
T3 = zeros(n,n,n);


%read edges 


for k=1:length(E2)
  E = E2(k,:);
  %update edge list
  W2(E(1),E(2))=1; % only fill uppder triangle
end

for j = 1:length(E3)
  
  T = E3(j,:);

  T3(T(1),T(2),T(3)) = 1; %tensor
  T3(T(1),T(3),T(2)) = 1; %tensor
  T3(T(2),T(1),T(3)) = 1; %tensor
  T3(T(2),T(3),T(1)) = 1; %tensor
  T3(T(3),T(1),T(2)) = 1; %tensor
  T3(T(3),T(2),T(1)) = 1; %tensor
      
      
  W3(T(1),T(2))= W3(T(2),T(2))+1;
  W3(T(1),T(3))= W3(T(1),T(3))+1;
  W3(T(2),T(3))= W3(T(2),T(3))+1;
end
%symmetryze matrix
W2 = W2 + W2'; W3 = W3 + W3'; 


end


