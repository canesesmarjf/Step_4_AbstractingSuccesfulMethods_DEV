% Step 2c: Review the results produced by the binary tree.
% Here we focus on the results produed by main_2.cpp

clear all
close all
clc

saveFig  = 1;
saveData = 0;

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
ip_a  = csvread(folderName + "ip_main_2a.csv") + 1;
ip_b  = csvread(folderName + "ip_main_2b.csv") + 1;
ip_c  = csvread(folderName + "ip_main_2c.csv") + 1;

hold on
ip = 10318;
hip(1) = plot3(+1   ,v_p(ip,1)   ,v_p(ip,2)   ,'mo','MarkerSize',5);
plot3(+1   ,v_p(ip,1)   ,v_p(ip,2)   ,'m.','MarkerSize',10);
hip(2) = plot3(x_p(ip_a) ,v_p(ip_a,1) ,v_p(ip_a,2) ,'r.','MarkerSize',5);
hip(3) = plot3(x_p(ip_b) ,v_p(ip_b,1) ,v_p(ip_b,2) ,'g.','MarkerSize',5);
hip(4) = plot3(x_p(ip_c) ,v_p(ip_c,1) ,v_p(ip_c,2) ,'bl.','MarkerSize',5);

hL = legend(hip,'index = 10318','index, dimensionality = 1',...
    'index, dimensionality = 2','xq = 0.25e6 (dimensionality = 1)');
hL.Location = 'southoutside';
view([90,0])

if saveFig
    folderName = "figures/";
    caseName = "t_" + num2str(tt);
    baseName = "Step_2c_main_2_";
    figureName = [baseName + caseName + "_species_" + num2str(ss)];

    % TIFF figure
    % exportgraphics(gcf,[folderName,figureName,'.tiff'],'Resolution',600);

    % PDF figure
    exportgraphics(gcf,[folderName + figureName + ".pdf"],'Resolution',600);

    % .fig file
    savefig([folderName + figureName + ".fig"])
end