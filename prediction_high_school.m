clear all 
clc
close all

% folder that contains algorithm scripts
addpath 'functions';
set(groot,'defaultFigureVisible','off') % 'on' to turn back on.  
tic
gamma_array = 1:50:500;
c3_array = 0:1:10;
n_eig = 3;
for train_ratio = [0.2, 0.6, 0.8]
    for data_name = "contact-high-school"
        PredictTriangles(data_name, train_ratio, c3_array, gamma_array, "false", n_eig)      
    end
end
toc
load chirp.mat
sound(y)