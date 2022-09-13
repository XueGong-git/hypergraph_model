function  [x, W2, W3, T3, data_type, label] = GenerateHygraph(n, K, gamma, c2, c3, input_shape)
        m = n/K; % number of nodes per cluster
        label = sort(repmat(1:1:K,1, m)); % angles from -pi to pi
        idx_rand = randperm(n);% shuffle the nodes
        label = label(idx_rand);
        if input_shape == "linear" 
            a = 0.05;

            x = sort(repmat(linspace(0,2, K),1,m)+(2*a*rand(1,n)-a));
            x = x(idx_rand);
            %x = sort(repmat(linspace(-pi,pi,K),1,m)+(2*a*rand(1,n)-a)); % angles from -pi to pi
            data_type = "cluster";

            [~, W2, W3, T3] = GenerateLinearHypergraph(x, gamma, c2, c3, 'linear');

        elseif input_shape == "periodic" 
            a = 0.05*pi; % noise;

    
            x = sort(repmat(linspace(-pi,pi,K),1,m)+(2*a*rand(1,n)-a)); % angles from -pi to pi
            x = x(idx_rand);
            data_type = "cluster";
  
            [W2, W3, T3] = GeneratePeriodicHypergraph(x, gamma, c2, c3);

        end
        
        
end