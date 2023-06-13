% Main_1 post process:
clear all
close all
clc

% NEED TO HAVE FIGURES SAVED AND ALSO CLEAN UP THE SCRIPT SO THAT IS FLOWS
% MORE LOGICALLY

save_fig  = 0;
save_data = 0;

% Choose time slice and species:
% =========================================================================
ss = 2;
tt = 12;
tt = 25;
tt = 37;
tt = 50;
scenario = "ss_" + string(ss) + "_tt_" + string(tt);

% Choose PICOS++ case:
% =========================================================================
picos_case = "PICOS_case_2/"; 

% Programatically define output data destination:
% =========================================================================
script_path = which(mfilename);
[~,script_name,~] = fileparts(script_path);
end_token = strfind(script_name,'_');
root_output = script_name(1:end_token(2)) + "output";

% Create directory to hold script output data if it doesnt exist:
output_dir = root_output + "/" + picos_case + scenario;
if (isempty(dir(output_dir)) == 1)
    mkdir(output_dir);
end

% Load particle data PRIOR to resampling (original):
% =========================================================================
data_folder = "./Step_1_output/";
x_p = readmatrix(data_folder + picos_case + "x_p_" + scenario + ".csv");
v_p = readmatrix(data_folder + picos_case + "v_p_" + scenario +".csv");
a_p = readmatrix(data_folder + picos_case + "a_p_" + scenario + ".csv");
clear data_folder;

% Normalize data:
x_norm = ceil(max(x_p));
v_norm = max(max(v_p));

if (ss == 1)
    v_norm = 2e6;
end

x_p = x_p/x_norm;
v_p = v_p/v_norm;

% Load particle data resampled by C++ code:
% =========================================================================
data_folder = "./Step_2_output/" + picos_case + scenario + "/";
x_pn = readmatrix(data_folder + "x_p_new" + ".csv");
v_pn = readmatrix(data_folder + "v_p_new" + ".csv");
a_pn = readmatrix(data_folder + "a_p_new" + ".csv");
clear data_folder;

% Normalize data:
x_pn = x_pn/x_norm;
v_pn = v_pn/v_norm;

% Load quadtree data from RESAMPLED data produced by C++ code:
% =========================================================================
particle_count = cell(1,64);
node_center    = cell(1,64);
node_dim       = cell(1,64);

data_folder = "./Step_2_output/" + picos_case + scenario + "/";
for xx = 1:64
    try
        particle_count{xx} = readmatrix(data_folder + "leaf_v_" + "p_count" ...
            + "_xx_" + string(xx) + ".csv");
        node_center{xx} = readmatrix(data_folder + "leaf_v_" + "node_center" ...
            + "_xx_" + string(xx) + ".csv");
        node_dim{xx} = readmatrix(data_folder + "leaf_v_" + "node_dim" ...
            + "_xx_" + string(xx) + ".csv");
    catch
        continue;
    end
end
x_q = readmatrix(data_folder + "x_q" + ".csv");  
clear data_folder;

% Computational particle density profile:
% -------------------------------------------------------------------------
data_folder = "./Step_2_output/" + picos_case + scenario + "/";
p_count = readmatrix(data_folder + "leaf_x_p_count" + ".csv");
% p_count_new = readmatrix(data_folder + "leaf_x_p_count_new" + ".csv");
clear data_folder;

mean_p_count = ceil(mean(p_count));

% Checking conservation of mass, momentum and energy:
% =========================================================================
% Here we perform a check conservation laws over the entire distribution
% before and after the resampling:

% Get data from C++ code:
data_folder = "./Step_2_output/" + picos_case + scenario + "/";

m_t = readmatrix(data_folder + "m_profile" + ".csv");
p_x = readmatrix(data_folder + "p_x_profile" + ".csv");
p_r = readmatrix(data_folder + "p_r_profile" + ".csv");
KE  = readmatrix(data_folder + "KE_profile" + ".csv");

m_tn = readmatrix(data_folder + "m_new_profile" + ".csv");
p_xn = readmatrix(data_folder + "p_x_new_profile" + ".csv");
p_rn = readmatrix(data_folder + "p_r_new_profile" + ".csv");
KEn  = readmatrix(data_folder + "KE_new_profile" + ".csv");

path = "./" + output_dir + "/";
file_name = "Conservation.txt";
diary off;
! rm Conservation.txt
diary(path + file_name);

% Global calculation:
% Conservation of mass:
mass_0 = sum(a_p);
mass_1 = sum(a_pn);
disp(" mass of original distribution is  " + string(mass_0));
disp(" mass of resampled distribution is " + string(mass_1));
disp("Error: " + string(100*abs(mass_1 - mass_0)/mass_0) + " %")
disp(newline)

% Conservation of momentum:
px_0 = sum(a_p.*v_p(:,1));
pr_0 = sum(a_p.*v_p(:,2));
px_1 = sum(a_pn.*v_pn(:,1));
pr_1 = sum(a_pn.*v_pn(:,2));

disp(" p_x of original distribution is  " + string(px_0));
disp(" p_x of resampled distribution is " + string(px_1));
disp("Error: " + string(100*abs(px_1 - px_0)/px_0) + " %")
disp(newline)

disp(" p_r of original distribution is  " + string(pr_0));
disp(" p_r of resampled distribution is " + string(pr_1));
disp("Error: " + string(100*abs(pr_1 - pr_0)/pr_0) + " %")
disp(newline)

% Conservation of energy:
K_0 = sum( a_p.*v_p(:,1).^2  +  a_p.*v_p(:,2).^2);
K_1 = sum(a_pn.*v_pn(:,1).^2 + a_pn.*v_pn(:,2).^2);
disp(" KE of original distribution is  " + string(K_0));
disp(" KE of resampled distribution is " + string(K_1));
disp("Error: " + string(100*abs(K_1 - K_0)/K_0) + " %")
disp(newline)

diary off;

% Plot data:
% =========================================================================

% Derived quantities:
dx = mean(diff(x_q));

% Set the minimum particle count to use for reampling a node:
min_count = 144;

% Grid depths:
kx = 6;
ky = 6;
kz = 6;
n_skip = 3;

% Compare Quadtree data with particle data:
% -------------------------------------------------------------------------
for xx = 1:numel(particle_count)

    if (isempty(particle_count{xx}))
        continue
    end
    
    % Calculate particle rng belonging to this node:
    nn = xx + 1;
    z1 = x_q(nn) - dx/2;
    z2 = x_q(nn) + dx/2;
    rng = find(x_p > z1 & x_p < z2);
    rng_new = find(x_pn > z1 & x_pn < z2);
    
    figure('color','w')
    plot_increase_size(2,2)
    hold on
    box on
    for vv = 1:numel(particle_count{xx})
        hsq(vv) = plotSquare(node_center{xx}(vv,:), node_dim{xx}(vv,:), particle_count{xx}(vv),min_count);
        set(hsq(vv),'lineWidth',1);
    end
    hold on
    plot(v_p(rng,1),v_p(rng,2),'k.','MarkerSize',15);
    plot(v_pn(rng_new,1),v_pn(rng_new,2),'g.','MarkerSize',8);
    title("xx = " + string(xx) + ", x = " + string(x_q(nn)) + ...
        ", Counts = " + string(numel(rng_new)) + ...
        ", surplus = " + string(numel(rng_new) - mean_p_count));
    plot_binary_tree_grid(gca,kx,ky,kz,n_skip);
    axis image
    xlim([-1,1])
    ylim([0,1])

    % Save figure:
    if save_fig
        path = "./" + output_dir + "/";
        figure_name = "QuadTree_xx_" + string(xx);
        save_figure(path,figure_name,"pdf",600)
    end
end

figure('color','w');
plot_increase_size(2,2)
box on
hold on
hb(1) = bar(x_q,p_count);
hb(2) = line([min(x_q),max(x_q)],[1,1]*mean_p_count,'LineWidth',3,'color','k');
hb(3) = plot(x_q,p_count_new,'r.-','LineWidth',3);
plot(x_q(nn),p_count(nn),'ro');
xlabel('$x/x_{norm}$ [m]','Interpreter','latex','FontSize',14)
ylabel('Counts','Interpreter','latex','FontSize',14)
legendText{1} = "Original";
legendText{2} = "Mean";
legendText{3} = "Resampled";
hL = legend(hb,legendText);
set(hL,'Interpreter','Latex')

% Save figure:
if save_fig
    path = "./" + output_dir + "/";
    figure_name = "ncp_m";
    save_figure(path,figure_name,"pdf",300)
end

% %% Checking conservation of mass, momentum and energy:
% % Here we perform a check conservation laws over the entire distribution
% % before and after the resampling:
% 
% % Get data from C++ code:
% data_folder = "./Step_2_output/" + picos_case + scenario + "/";
% 
% m_t = readmatrix(data_folder + "m_profile" + ".csv");
% p_x = readmatrix(data_folder + "p_x_profile" + ".csv");
% p_r = readmatrix(data_folder + "p_r_profile" + ".csv");
% KE  = readmatrix(data_folder + "KE_profile" + ".csv");
% 
% m_tn = readmatrix(data_folder + "m_new_profile" + ".csv");
% p_xn = readmatrix(data_folder + "p_x_new_profile" + ".csv");
% p_rn = readmatrix(data_folder + "p_r_new_profile" + ".csv");
% KEn  = readmatrix(data_folder + "KE_new_profile" + ".csv");
% 
% path = "./" + output_dir + "/";
% file_name = "Conservation.txt";
% diary off;
% ! rm Conservation.txt
% diary(path + file_name);
% 
% % Global calculation:
% % Conservation of mass:
% mass_0 = sum(a_p);
% mass_1 = sum(a_pn);
% disp(" mass of original distribution is  " + string(mass_0));
% disp(" mass of resampled distribution is " + string(mass_1));
% disp("Error: " + string(100*abs(mass_1 - mass_0)/mass_0) + " %")
% disp(newline)
% 
% % Conservation of momentum:
% px_0 = sum(a_p.*v_p(:,1));
% pr_0 = sum(a_p.*v_p(:,2));
% px_1 = sum(a_pn.*v_pn(:,1));
% pr_1 = sum(a_pn.*v_pn(:,2));
% 
% disp(" p_x of original distribution is  " + string(px_0));
% disp(" p_x of resampled distribution is " + string(px_1));
% disp("Error: " + string(100*abs(px_1 - px_0)/px_0) + " %")
% disp(newline)
% 
% disp(" p_r of original distribution is  " + string(pr_0));
% disp(" p_r of resampled distribution is " + string(pr_1));
% disp("Error: " + string(100*abs(pr_1 - pr_0)/pr_0) + " %")
% disp(newline)
% 
% % Conservation of energy:
% K_0 = sum( a_p.*v_p(:,1).^2  +  a_p.*v_p(:,2).^2);
% K_1 = sum(a_pn.*v_pn(:,1).^2 + a_pn.*v_pn(:,2).^2);
% disp(" KE of original distribution is  " + string(K_0));
% disp(" KE of resampled distribution is " + string(K_1));
% disp("Error: " + string(100*abs(K_1 - K_0)/K_0) + " %")
% disp(newline)
% 
% diary off;

return;

%% Comparing phase space data:

% Plot results:
% =========================================================================
% rng = find(x_pn ~= -1);

IONS.x_p = x_p;
IONS.v_p = v_p;
IONS_new.x_p = x_pn;
IONS_new.v_p = v_pn;

var{1} = IONS;
var{2} = IONS_new;

kx = 5;
ky = 6;
kz = 6;
n_skip = 3;

for rr = 1:2
    figure('color','w');
    plot_increase_size(1.5,1.5);

    % Set up the first subplot (3D view)
    subplot(2,2,1);
    box on;
    plot3(var{rr}.x_p,var{rr}.v_p(:,1),var{rr}.v_p(:,2),'k.','MarkerSize',1);
    xlabel('$x$','Interpreter','latex','FontSize',18);
    ylabel('$v_{\parallel}$','Interpreter','latex','FontSize',18);
    zlabel('$v_{\perp}$','Interpreter','latex','FontSize',18);
    view([-30,30]);
    xlim([-1,1])
    ylim([-1,1])
    zlim([0,1])
    title('3D view','Interpreter','latex','FontSize',14);
    plot_binary_tree_grid(gca,kx,ky,kz,n_skip,false);

    % Set up the second subplot (xy view)
    subplot(2,2,2);
    box on;
    plot(var{rr}.x_p,var{rr}.v_p(:,1),'k.','MarkerSize',1);
    xlabel('$x$','Interpreter','latex','FontSize',18);
    ylabel('$v_{\parallel}$','Interpreter','latex','FontSize',18);
    title('xy view','Interpreter','latex','FontSize',14);
    plot_binary_tree_grid(gca,kx,ky,kz,n_skip,false);
    
    % Set up the third subplot (xz view)
    subplot(2,2,3);
    box on;
    plot(var{rr}.x_p,var{rr}.v_p(:,2),'k.','MarkerSize',1);
    xlabel('$x$','Interpreter','latex','FontSize',18);
    ylabel('$v_{\perp}$','Interpreter','latex','FontSize',18);
    title('xz view','Interpreter','latex','FontSize',14);
    plot_binary_tree_grid(gca,kx,ky,kz,n_skip,false);
    
    % Set up the fourth subplot (yz view)
    subplot(2,2,4);
    box on;
    plot(var{rr}.v_p(:,1),var{rr}.v_p(:,2),'k.','MarkerSize',1);
    xlabel('$v_{\parallel}$','Interpreter','latex','FontSize',18);
    ylabel('$v_{\perp}$','Interpreter','latex','FontSize',18);
    title('yz view','Interpreter','latex','FontSize',14);
    plot_binary_tree_grid(gca,kx,ky,kz,n_skip,true);
end

%% Functions:
% =========================================================================
function [] = plot_binary_tree_grid(ax,kx,ky,kz,n_skip,verbose)
    % Default value:
    if nargin < 6
        verbose = false;
    end

    % Produce grid:
    grids.x = linspace(-1,+1,2^kx + 1);
    grids.y = linspace(-1,+1,2^ky +1);
    grids.z = linspace(0,+1,2^kz + 1);
    
    % Execute if verbose is true
    if verbose
        disp("Total nodes: " + num2str(2^(kx+ky+kz)))
    end

    % Draw grid:
    ax.XTick = grids.x(1:1:end);
    ax.YTick = grids.y(1:1:end); ax.YTickLabel = grids.y(1:1:end);
    ax.ZTick = grids.z(1:1:end); ax.ZTickLabel = grids.z(1:1:end);
    grid on;

    decimalPlaces = 3;

    ticks = get(gca,'Xtick');
    str = num2str(ticks',['%.' num2str(decimalPlaces) 'f']);
    ticklabels = cellstr(str);
    rng = mod(1:numel(ticklabels),n_skip) ~= 0;
    ticklabels(rng) = {''};
    set(gca,'XTickLabel',ticklabels);

    ticks = get(gca,'Ytick');
    str = num2str(ticks',['%.' num2str(decimalPlaces) 'f']);
    ticklabels = cellstr(str);    rng = mod(1:numel(ticklabels),n_skip) ~= 0;
    ticklabels(rng) = {''};
    set(gca,'YTickLabel',ticklabels);
   
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

function [] = plot_increase_size(sfx,sfy)
    set(gcf,'Position',get(gcf,'Position').*[1 1 sfx sfy]);
end

function h = plotSquare(center, dimensions, particle_count, min_count)
    % center: a vector [x, y] representing the coordinates of the center
    % dimensions: a vector [dx, dy] representing the width and height of the square

    % Extract center coordinates
    x = center(1);
    y = center(2);

    % Extract dimensions
    dx = dimensions(1);
    dy = dimensions(2);

    % Min dim:
    dx_min = min(dx);
    dy_min = min(dy);

    % Calculate the coordinates of the square's vertices
    verticesX = [x-dx/2, x+dx/2, x+dx/2, x-dx/2, x-dx/2];
    verticesY = [y-dy/2, y-dy/2, y+dy/2, y+dy/2, y-dy/2];

    % Plot the square
    hold on;  % Optional: maintain existing plot
    if (particle_count > min_count)
       h = plot(verticesX, verticesY, 'r-', 'LineWidth', 4);
       fill(verticesX, verticesY, 'r', 'LineWidth',1);
    else
       h = plot(verticesX, verticesY, 'b-', 'LineWidth', 1);
    end

    % Annotate particle_count{xx}:
    text(x,y,num2str(particle_count),'HorizontalAlignment',...
        'center','VerticalAlignment','middle','FontSize',6)

    hold off; % Optional: release hold on existing plot

    % Optional: Adjust figure axes if necessary
    axis equal; % Set equal scaling for x and y axes
    % Adjust xlim and ylim as needed to fit the square within the figure
    % xlim([xmin, xmax]);
    % ylim([ymin, ymax]);

    % Optional: Add labels or title to the plot
    % xlabel('X');
    % ylabel('Y');
    % title('Square Plot');
end