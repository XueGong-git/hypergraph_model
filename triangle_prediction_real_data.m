clear all 
clc
close all

% folder that contains algorithm scripts
addpath 'functions';
set(groot,'defaultFigureVisible','off') % 'on' to turn back on.  
tic
gamma_array = 1:50:500;
c3_array = 1:1:10;
n_eig = 3;
for train_ratio = [0.8]
    for data_name = "contact-high-school"
        PredictTriangles(data_name, train_ratio, c3_array, gamma_array, "true", n_eig)      
    end
end
toc
load train.mat
sound(y)