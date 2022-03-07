function  [x, W2, W3, T3, data_type] = GenerateHygraph(n, K, gamma, c2, c3, input_shape, input)
        m = n/K; % number of nodes per cluster

        if input_shape == "linear" 
            a = 0.05;

            switch input

            case 1 
                x = linspace(0,2,n);%uniform from 0 to 2pi
                data_type = "uniform";
            case 2 
                x = sort(repmat(linspace(0,2, K),1,m)+(2*a*rand(1,n)-a));
                %x = sort(repmat(linspace(-pi,pi,K),1,m)+(2*a*rand(1,n)-a)); % angles from -pi to pi
                data_type = "cluster";
            case 3
                x = linspace(0,0,n); %overlapping x
                data_type = "overlap";
            end
            [W2, W3, T3] = GenerateLinearHypergraph(x, gamma, c2, c3);

        elseif input_shape == "periodic" 
            a = 0.05*pi; % noise;

            switch input
    
            case 1 
                x = linspace(0,2*pi,n);%uniform from 0 to 2pi
                data_type = "uniform";
            case 2 
                x = sort(repmat(linspace(-pi,pi,K),1,m)+(2*a*rand(1,n)-a)); % angles from -pi to pi
                data_type = "cluster";
            case 3
                x = linspace(0,0,n); %overlapping x
                data_type = "overlap";
            end
            [W2, W3, T3] = GeneratePeriodicHypergraph(x, gamma, c2, c3);

        end
end