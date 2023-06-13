% Step_1: Preview data produced by PICOS++
% We extract all the data. This includes from all MPI processes

clear all
close all
clc

save_fig  = 1;
save_data = 1;

% Choose PICOS++ case:
% =========================================================================
% picos_case = 'PICOS_case_0/'; 
% picos_case = 'PICOS_case_1/'; 
picos_case = 'PICOS_case_2/'; 

% Programatically define output data destination:
% =========================================================================
script_path = which(mfilename);
[~,script_name,~] = fileparts(script_path);
end_token = strfind(script_name,'_');
root_output = [script_name(1:end_token(2)),'output'];

% Create directory to hold script output data if it doesnt exist:
output_dir = [root_output,'/',picos_case];
if (isempty(dir(output_dir)) == 1)
    mkdir(output_dir);
end

% Print information to CLI:
disp("Script running: " + newline + string(script_name));
disp(newline)
disp("Input data selected for this script: " + newline + string(picos_case))
disp(newline);
disp("Output data destination:" + newline + string(output_dir))
disp(newline)

% Physical constants:
% =========================================================================
e_c = 1.6020e-19;
k_B = 1.3806e-23;
m_p = 1.6726e-27;
m_e = 9.1094e-31;
mu0 = 4*pi*1e-7;
c   = 299792458;
E_0 = m_p*c^2;

% Import data from HDF5 dataset produced by PICOS++:
% =========================================================================
% Includes:
addpath('/home/jfcm/Repos/PICOS/picosFILES/matlabFiles');
       
% Extract data:
root = picos_case;  
fileName = [root,'HDF5/main.h5'];
extractDataFromH5;
clear root

%% Preview the data:

% Select species and time slice::
ss = 2;
tt = 12;
tt = 25;
tt = 37;
tt = 50;

% Plasma density:
figure('color','w')
plot_increase_size(2,1)

subplot(1,2,1)
box on
mesh(x_m,t_p*1e3,n_m{ss}')
view([150,45])
xlim([-2,2])
xlabel('x [m]')
ylabel('t [ms]')
title("$n_m$ [m$^{-3}$]",'Interpreter','latex','FontSize',14)

% Computational particle density:
subplot(1,2,2)
box on
mesh(x_m,t_p*1e3,ncp_m{ss}')
view([150,45])
xlim([-2,2])
xlabel('x [m]')
ylabel('t [ms]')
title("$n_m^{cp}$ [m$^{-1}$]",'Interpreter','latex','FontSize',14)

if save_fig
    path = "./" + string(output_dir);
    figure_name = "mesh" + "_ss_" + num2str(ss)+ "_tt_" + num2str(tt);
    save_figure(path,figure_name,"pdf",250)
end

mean_ncp = mean(ncp_m{ss}(:,tt));

% Plot main results
figure('color','w')
plot_increase_size(2,1)

subplot(1,2,1)
box on
hold on
plot(x_m,n_m{ss}(:,tt),'k','LineWidth',2)
plot(x_m,Bx_m(:,tt)*4e17,'r','LineWidth',2)
xlabel("$x_m$ [m]",'Interpreter','latex','FontSize',14)
title("$n_m$ [m$^{-3}$]",'Interpreter','latex','FontSize',14)
xlim([-2,2])

subplot(1,2,2)
box on
hold on
plot(x_m,ncp_m{ss}(:,tt),'k','LineWidth',2)
plot(x_m,Bx_m(:,tt)*4e3,'r','LineWidth',2)
line([-2,2],[1,1]*mean_ncp,'LineWidth',2)
xlabel("$x_m$ [m]",'Interpreter','latex','FontSize',14)
title("$n_m^{cp}$ [m$^{-1}$]",'Interpreter','latex','FontSize',14)
xlim([-2,2])

if save_fig
    path = "./" + string(output_dir);
    figure_name = "plot" + "_ss_" + num2str(ss)+ "_tt_" + num2str(tt);
    save_figure(path,figure_name,"pdf",250)
end

%% Save data for C++ binary-quad tree code:
% x_p, v_p and a_p:
% It important to note that he we are only extracting data from a SINGLE
% MPI process.
% This is relevant since the binary tree algorithm will be used in
% precisely this manner: for each MPI process individually

% Get particle data:
fileName = "./" + picos_case + "HDF5/PARTICLES_FILE_0.h5";
info = hdf5info(fileName);

dataset = ["/" + num2str(tt) + "/ions/species_" + num2str(ss) + "/x_p"];
x_p = hdf5read(fileName,dataset);

dataset = ["/" + num2str(tt) + "/ions/species_" + num2str(ss) + "/v_p"];
v_p = hdf5read(fileName,dataset);

dataset = ["/" + num2str(tt) + "/ions/species_" + num2str(ss) + "/a_p"];
a_p = hdf5read(fileName,dataset);

% save data to tables:
file_name = "_ss_" + num2str(ss)+ "_tt_" + num2str(tt) + ".csv";

if save_data
    % Save data to tables:
    writematrix(x_p,output_dir + "x_p" + file_name)
    writematrix(v_p,output_dir + "v_p" + file_name)
    writematrix(a_p,output_dir + "a_p" + file_name)
end

%% Plot phase space of computational particles:
% =========================================================================
figure('color','w');
plot_increase_size(1.5,1.5);
kx = 5;
ky = 6;
kz = 5;

x_norm = double(main.mesh.Lx_max);
v_norm = max(max(v_p(:,2)));

% Set up the first subplot (3D view)
subplot(2,2,1);
box on;
plot3(x_p/x_norm,v_p(:,1)/v_norm,v_p(:,2)/v_norm,'k.','MarkerSize',1);
xlabel('$x/x_{norm}$','Interpreter','latex','FontSize',18);
ylabel('$v_{\parallel}$','Interpreter','latex','FontSize',18);
zlabel('$v_{\perp}$','Interpreter','latex','FontSize',18);
view([-30,30]);
title("3D view, $x_{norm}$ = " + num2str(x_norm),'Interpreter','latex','FontSize',14);
set(gca,'FontSize',9)

% Set up the second subplot (xy view)
subplot(2,2,2);
box on;
plot(x_p/x_norm,v_p(:,1)/v_norm,'k.','MarkerSize',1);
xlabel('$x/x_{norm}$','Interpreter','latex','FontSize',18);
ylabel('$v_{\parallel}$','Interpreter','latex','FontSize',18);
title('xy view','Interpreter','latex','FontSize',14);
plot_binary_tree_grid(gca,kx,ky,kz);
set(gca,'FontSize',7)
xlim([-1,1])

% Set up the third subplot (xz view)
subplot(2,2,3);
box on;
plot(x_p/x_norm,v_p(:,2)/v_norm,'k.','MarkerSize',1);
xlabel('$x/x_{norm}$','Interpreter','latex','FontSize',18);
ylabel('$v_{\perp}$','Interpreter','latex','FontSize',18);
title('xz view','Interpreter','latex','FontSize',14);
plot_binary_tree_grid(gca,kx,ky,kz);
set(gca,'FontSize',7)
xlim([-1,1])

% Set up the fourth subplot (yz view)
subplot(2,2,4);
box on;
plot(v_p(:,1)/v_norm,v_p(:,2)/v_norm,'k.','MarkerSize',1);
xlabel('$v_{\parallel}$','Interpreter','latex','FontSize',18);
ylabel('$v_{\perp}$','Interpreter','latex','FontSize',18);
title('yz view','Interpreter','latex','FontSize',14);
plot_binary_tree_grid(gca,kx,ky,kz);
set(gca,'FontSize',7)
xlim([-1,1])

if save_fig
    path = "./" + string(output_dir);
    figure_name = "phase_space" + "_ss_" + num2str(ss)+ "_tt_" + num2str(tt);
    save_figure(path,figure_name,"pdf",250)
end

%% Functions:
function [] = plot_increase_size(sfx,sfy)
    set(gcf,'Position',get(gcf,'Position').*[1 1 sfx sfy]);
end

function [] = save_figure(path,figure_name,format,resolution)
    switch lower(format)
        case 'tiff'
            exportgraphics(gcf,[path + figure_name + ".tiff"],'Resolution',resolution);
        case 'pdf'
            exportgraphics(gcf,[path + figure_name + ".pdf"],'Resolution',resolution);
        case 'fig'
            savefig([path + figure_name + ".fig"])
        otherwise
            error('Invalid format. Choose "pdf", "tiff" or "fig');
    end
end

function [] = plot_binary_tree_grid(ax,kx,ky,kz)
    % Produce grid:
    grids.x = linspace(-1,+1,2^kx + 1);
    grids.y = linspace(-1,+1,2^ky +1);
    grids.z = linspace(0,+1,2^kz + 1);
    
    disp("Total nodes: " + num2str(2^(kx+ky+kz)))
    
    % Draw grid:
    ax.XTick = grids.x(1:1:end); ax.XTickLabel = grids.x(1:1:end);
    ax.YTick = grids.y(1:1:end); ax.YTickLabel = grids.y(1:1:end);
    ax.ZTick = grids.z(1:1:end); ax.ZTickLabel = grids.z(1:1:end);
    grid on;
end
