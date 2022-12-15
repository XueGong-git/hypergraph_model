clear all 
clc
close all

% folder that contains algorithm scripts
addpath 'functions';
set(groot,'defaultFigureVisible','off') % 'on' to turn back on.  
tic
gamma_array = 1:10:100;
c3_array = 1:1:10;
for train_ratio = [0.1, 0.2, 0.6, 0.8]
    for data_name = "email-Eu"
        PredictTriangles(data_name, train_ratio, c3_array, gamma_array, "true")      
    end
end