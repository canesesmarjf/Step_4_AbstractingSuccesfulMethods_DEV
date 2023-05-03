% Step 1: View HDF5 data and choose time slice: In this script, we extract
% the data from the HDF5 files produced by PICOS++ and plot it so that we
% may select the time slice and species that We need. Once we are happy
% with a particular data set, it can be saved to .csv data which becomes
% the input for the C++ binary tree code. In the binary tree code we can
% select the location of phase space to obtain the data which is then
% reviewed in Step_2.

clear all
close all
clc

saveFig  = 0;
saveData = 0;

% Select time slices:
tt = 50;

% Select species:
ss = 2;

% Import data from HDF5 dataset produced by PICOS++:
% =========================================================================
% x_p, v_p and a_p:
% It important to note that he we are only extracting data from a SINGLE
% MPI process.
% This is relevant since the binary tree algorithm will be used in
% precisely this manner: for each MPI process individually

fileName = './case_0/HDF5/PARTICLES_FILE_0.h5';
info = hdf5info(fileName);

dataset = ["/" + num2str(tt) + "/ions/species_" + num2str(ss) + "/x_p"];
x_p = hdf5read(fileName,dataset);

dataset = ["/" + num2str(tt) + "/ions/species_" + num2str(ss) + "/v_p"];
v_p = hdf5read(fileName,dataset);

dataset = ["/" + num2str(tt) + "/ions/species_" + num2str(ss) + "/a_p"];
a_p = hdf5read(fileName,dataset);

% Create 3D grid:
% =========================================================================
% Calculate maximum speed:
% This is needed to define grid size
v_max = max(max(sqrt( v_p(:,1).^2 + v_p(:,2).^2 )));
v_max = ceil(v_max/1e5)*1e5

% Select grid parameters:
kx = 4;
grid.x = linspace(-1,+1,2^kx + 1);
ky = 5;
grid.y = linspace(-v_max,v_max,2^ky +1);
kz = 4;
grid.z = linspace(0,v_max,2^kz + 1);

% Total number of nodes:
disp("Total nodes: " + num2str(2^(kx+ky+kz)))

% Plot data to choose a point from:
if 0
    figure; 
    plot3(v_p(:,1),v_p(:,2),1:numel(x_p),'k.','MarkerSize',2)
end

% Plot raw data:
% =========================================================================
figure('color','w')
plot3(x_p,v_p(:,1),v_p(:,2),'k.','MarkerSize',2)
ax = gca;
ylim([-1,+1]*v_max)
zlim([-0,+1]*v_max)
ax.XTick = grid.x;
ax.YTick = grid.y;
ax.ZTick = grid.z;
ax.PlotBoxAspectRatio = [1 2 1];
xlabel(ax,'$x$','Interpreter','latex','FontSize',15);
ylabel(ax,'$v_{\parallel}$','Interpreter','latex','FontSize',15);
zlabel(ax,'$v_{\perp}$','Interpreter','latex','FontSize',15);
grid on

% Save figure:
if saveFig
    folderName = "";
    caseName = "t_" + num2str(tt);
    baseName = "Step_1_ChoosingData_";
    figureName = [baseName + caseName + "_species_" + num2str(ss)];

    % TIFF figure
    % exportgraphics(gcf,[folderName,figureName,'.tiff'],'Resolution',600);

    % PDF figure
    exportgraphics(gcf,[folderName + figureName + ".pdf"],'Resolution',600);

    % .fig file
    savefig([folderName + figureName + ".fig"])
end

% Save data so it can be used in the binary tree c++ code:
% =========================================================================
if saveData
    % Save data to tables:
    target_file = "./input_files/"
    writematrix(x_p,target_file + "Step_1_x_p.csv")
    writematrix(v_p,target_file + "Step_1_v_p.csv")
    writematrix(a_p,target_file + "Step_1_a_p.csv")
end