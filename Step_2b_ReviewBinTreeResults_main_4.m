% Step 1: View HDF5 data and choose time slice:
% This script focuses on the results produced by main_2.cpp and main_3.cpp
% To get the resutls from main_4.cpp you need to run the
% Step_2b_ReviewBinTreeResults_main_4.m

clear all
close all
clc

saveFig  = 1;
saveData = 1;

% Select time slices:
tt = 50;

% Select species:
ss = 2;

% Open figure with raw data:
% =========================================================================
folderName = "figures/";
caseName = "t_" + num2str(tt);
baseName = "Step_1_ChoosingData_";
figureName = [baseName + caseName + "_species_" + num2str(ss)];

% .fig file
open([folderName + figureName + ".fig"])

% Import dataset produced by PICOS++:
% =========================================================================
% Open and load data from Step_1:
target_file = "./input_files/";
x_p = readmatrix(target_file + "Step_1_x_p.csv");
v_p = readmatrix(target_file + "Step_1_v_p.csv");
a_p = readmatrix(target_file + "Step_1_a_p.csv");

% Run binary tree c++ code to produce index vector:
% =========================================================================
% ./run.sh 2
% This produces a series of outputs labelled ip_main_xx

% Load data produced by binary tree c++ executable:
% =========================================================================
folderName = "input_files/";
ip_a  = csvread(folderName + "ip_main_4a.csv") + 1;
ip_b  = csvread(folderName + "ip_main_4b.csv") + 1;
ip_c  = csvread(folderName + "ip_main_4c.csv") + 1;
ip_d  = csvread(folderName + "ip_main_4d.csv") + 1;
ip_e  = csvread(folderName + "ip_main_4e.csv") + 1;
ip_f  = csvread(folderName + "ip_main_4f.csv") + 1;

hold on
hip(1) = plot3(x_p(ip_a) ,v_p(ip_a,1) ,v_p(ip_a,2) ,'r.','MarkerSize',5);
hip(2) = plot3(x_p(ip_b) ,v_p(ip_b,1) ,v_p(ip_b,2) ,'g.','MarkerSize',5);
hip(3) = plot3(x_p(ip_c) ,v_p(ip_c,1) ,v_p(ip_c,2) ,'bl.','MarkerSize',5);
hip(4) = plot3(x_p(ip_d) ,v_p(ip_d,1) ,v_p(ip_d,2) ,'m.','MarkerSize',5);
hip(5) = plot3(x_p(ip_e) ,v_p(ip_e,1) ,v_p(ip_e,2) ,'c.','MarkerSize',5);
hip(6) = plot3(x_p(ip_f) ,v_p(ip_f,1) ,v_p(ip_f,2) ,'y.','MarkerSize',5);

hL = legend(hip,'a','b','c','d','e','f');

% Save figure:
if saveFig
    folderName = "figures/";
    caseName = "t_" + num2str(tt);
    baseName = "Step_2_BinTreeData_";
    figureName = [baseName + caseName + "_species_" + num2str(ss)];

    % TIFF figure
    % exportgraphics(gcf,[folderName,figureName,'.tiff'],'Resolution',600);

    % PDF figure
    exportgraphics(gcf,[folderName + figureName + ".pdf"],'Resolution',600);

    % .fig file
    savefig([folderName + figureName + ".fig"])
end