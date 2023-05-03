% Step 1: View HDF5 data and choose time slice:
% This script focuses on the results produced by main_2.cpp and main_3.cpp
% To get the resutls from main_4.cpp you need to run the
% Step_2b_ReviewBinTreeResults_main_4.m

clear all
close all
clc

saveFig  = 0;
saveData = 0;

% Select time slices:
tt = 50;

% Select species:
ss = 2;

% Open figure with raw data:
% =========================================================================
folderName = "";
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
ip_2D  = csvread('ip_main_2.csv') + 1;
ip_3D  = csvread('ip_main_3.csv') + 1;
ip_3Db = csvread('ip_main_3b.csv') + 1;
ip_3Dc = csvread('ip_main_3c.csv') + 1;

hold on
hip(1) = plot3(x_p(ip_2D) ,v_p(ip_2D,1) ,v_p(ip_2D,2) ,'r.','MarkerSize',5);
hip(2) = plot3(x_p(ip_3D) ,v_p(ip_3D,1) ,v_p(ip_3D,2) ,'g.','MarkerSize',5);
hip(3) = plot3(x_p(ip_3Db),v_p(ip_3Db,1),v_p(ip_3Db,2),'bl.','MarkerSize',5);
hip(4) = plot3(x_p(ip_3Dc),v_p(ip_3Dc,1),v_p(ip_3Dc,2),'m.','MarkerSize',5);

hL = legend(hip,'2','3','3b','3c');

% Save figure:
if saveFig
    folderName = "";
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