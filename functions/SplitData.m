% GenerateModelLikelihood  Using the given embedding, calculate the
% likelihood of unweighted hyper graph
%
% INPUTS
% 
%
%
% OUTPUTS
% 

function [train_list, test_list, n] = SplitData(data_name,train_ratio)

times = int64(importdata(strcat('raw_data/', data_name,'/', data_name, '-times.txt'))); 
nverts = importdata(strcat('raw_data/', data_name,'/', data_name, '-nverts.txt'));
simplices = importdata(strcat('raw_data/', data_name,'/', data_name, '-simplices.txt'));
n = max(unique(simplices));
train_list = {};
test_list = {};

curr_ind = 1;
train_end = round(quantile(times,train_ratio));
for i = 1:length(times)
    nvert = nverts(i); time = times(i);
    if time <= train_end
        train_list{end+1} = simplices(curr_ind:(curr_ind + nvert-1));
    else
        test_list{end+1} = simplices(curr_ind:(curr_ind + nvert-1));
    end
    curr_ind  = curr_ind+nvert;
end


end

