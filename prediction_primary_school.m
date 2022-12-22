clear all 
clc
close all

% folder that contains algorithm scripts
addpath 'functions';
set(groot,'defaultFigureVisible','off') % 'on' to turn back on.  
tic
c3_array = 0:0.1:2;
gamma_array = logspace(0,3.3,10);
n_eig = 4;
for train_ratio = [0.4 0.8]
    for data_name = "contact-primary-school"
        PredictTriangles(data_name, train_ratio, c3_array, gamma_array, "false", n_eig)      
    end
end
toc
load chirp.mat
sound(y)