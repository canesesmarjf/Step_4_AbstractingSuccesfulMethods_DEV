% Step 3: Use data produced by main_5
clear all
close all
clc

saveFig  = 0;
saveData = 0;

% Import dataset produced by PICOS++:
% =========================================================================
% Save data to tables:
target_file = "./input_files/";
x_p = readmatrix(target_file + "Step_1_x_p.csv");
v_p = readmatrix(target_file + "Step_1_v_p.csv");
a_p = readmatrix(target_file + "Step_1_a_p.csv");

% Get M set calculated by c++ code:
% =========================================================================
target_file = "./input_files/";
cwi_M = readmatrix(target_file + "wi_M.csv");
cxi_M = readmatrix(target_file + "xi_M.csv");
cyi_M = readmatrix(target_file + "yi_M.csv");
czi_M = readmatrix(target_file + "zi_M.csv");
ip_free = 1 + readmatrix(target_file + "ip__free_main_5a.csv");

% Select figure:
% =========================================================================
% Time slices:
tt = 50;

% Species:
ss = 2;

% Open figure with raw data:
% =========================================================================
folderName = "figures/";
caseName = "t_" + num2str(tt);
baseName = "Step_1_ChoosingData_";
figureName = [baseName + caseName + "_species_" + num2str(ss)];

% .fig file
open([folderName + figureName + ".fig"])

% Update grid:
v_max = max(max(sqrt( v_p(:,1).^2 + v_p(:,2).^2 )));
v_max = ceil(v_max/1e5)*1e5;
kx = 4+1-1;
grid.x = linspace(-1,+1,2^kx + 1);
ky = 5+1-1;
grid.y = linspace(-v_max,v_max,2^ky +1);
kz = 4+1-1;
grid.z = linspace(0,v_max,2^kz + 1);
ax = gca;
ax.XTick = grid.x;
ax.YTick = grid.y;
ax.ZTick = grid.z;

% Total number of nodes:
disp("Total nodes: " + num2str(2^(kx+ky+kz)))

% Test of 3D vranic method:
% =========================================================================
% Get index data:
ip = csvread(target_file + "ip_main_5" + "a" + ".csv") + 1;

% Select subset N:
N = round(numel(ip)/1);

% Plot merge cell data:
hold on
plot3(x_p(ip(1:N)) ,v_p(ip(1:N),1) ,v_p(ip(1:N),2) ,'r.','MarkerSize',5)

% Plot points to be repurposed:
plot3(x_p(ip_free) ,v_p(ip_free,1) ,v_p(ip_free,2) ,'ro','MarkerSize',5)

% Get data:
wi = a_p(ip);
xi = x_p(ip);
yi = v_p(ip,1);
zi = v_p(ip,2);

% Calculate set M:
[wi_M, xi_M,yi_M,zi_M] = down_sample(wi,xi,yi,zi);

% Calculate merge cell statistics:
cs = cell_stats(wi,xi,yi,zi,1);
cs_M = cell_stats(wi_M,xi_M,yi_M,zi_M,0);

% Calculate number of vranic cycles:
% N0 = 12;
% v_cycles = mod(numel(ip),N0);
% 
% k = 1;
% l = 1;
% for vv = 1:v_cycles
%     % Assemble N set data:
%     rng = k:(k+N0-1);
%     wi = a_p(ip(rng));
%     xi = x_p(ip(rng));
%     yi = v_p(ip(rng),1);
%     zi = v_p(ip(rng),2);
%     
%     % Calculate M set data:
%     rng_M = l:(l+6-1);
%     [wi_M(rng_M),xi_M(rng_M),yi_M(rng_M),zi_M(rng_M)] = down_sample(wi,xi,yi,zi);
%     
%     % Calculate merge cell statistics:
%     cs = cell_stats(wi,xi,yi,zi,1);
%     cs_M = cell_stats(wi_M(rng_M),xi_M(rng_M),yi_M(rng_M),zi_M(rng_M),0);
% 
%     % Increment:
%     k = rng(end) + 1;
%     l = rng_M(end) + 1;
% end

% Plot data:
hold on
plot3(xi_M,yi_M,zi_M,'g.','MarkerSize',5)
plot3(xi_M,yi_M,zi_M,'go','MarkerSize',5)

% Check conservation of mass, momentum and energy:
disp('Number of particles:');
disp("Set N = " + num2str(N))
disp("Set N = " + num2str(6))

disp('Conservation of mass:');
disp("Set N, w_t = " + num2str(cs.w_t))
disp("Set M, w_t = " + num2str(cs_M.w_t))

disp('Conservation of momentum:');
disp("Set N, E_x = " + num2str(cs.E_x))
disp("Set M, E_x = " + num2str(cs_M.E_x))
disp("Set N, E_y = " + num2str(cs.E_y))
disp("Set M, E_y = " + num2str(cs_M.E_y))
disp("Set N, E_z = " + num2str(cs.E_z))
disp("Set M, E_z = " + num2str(cs_M.E_z))
disp("Set N, E_r = " + num2str(cs.E_r))
disp("Set M, E_r = " + num2str(cs_M.E_r))

disp('Conservation of energy:');
disp("Set N, sigma_r = " + num2str(cs.sigma_r))
disp("Set M, sigma_r = " + num2str(cs_M.sigma_r))

plot3(cxi_M,cyi_M,czi_M,'bl.','MarkerSize',5)
%% Plot merge cell with coordinate system:
ii = 1;

figure('color','w');
hold on

% Plot merge-cell data:
dx = xi - cs.E_x;
dy = yi - cs.E_y;
dz = zi - cs.E_z;
plot3(dx,dy,dz,'k.')

% Plot particles to be repurposed:
dx = x_p(ip_free)   - cs.E_x;
dy = v_p(ip_free,1) - cs.E_y;
dz = v_p(ip_free,2) - cs.E_z;
plot3(dx,dy,dz,'ro')

% Plot coordinate system:
e_prime = basis_prime(cs);
basis_plot(e_prime(:,1),e_prime(:,2),e_prime(:,3),[1,1,1]*4e4)

% Plot M set particles:
plot3(xi_M - cs_M.E_x ,yi_M - cs_M.E_y, zi_M - cs_M.E_z,...
    'g.','MarkerSize',20)
axis image;
box on

%% Functions:
function [] = basis_plot(xhat,yhat,zhat,scale)
    hold on
    line([0,xhat(1)]*scale(1),[0,xhat(2)]*scale(2),[0,xhat(3)]*scale(3),'color','k')
    line([0,yhat(1)]*scale(1),[0,yhat(2)]*scale(2),[0,yhat(3)]*scale(3),'color','bl')
    line([0,zhat(1)]*scale(1),[0,zhat(2)]*scale(2),[0,zhat(3)]*scale(3),'color','r')
end

function [cs] = cell_stats(wi,xi,yi,zi,skew_flag)
    % Calculate probablity:
    cs.w_t = sum(wi);
    cs.p_i = wi/cs.w_t;
    
    % Calculate expectation values:
    cs.E_x  = dot(cs.p_i,xi);
    cs.E_y  = dot(cs.p_i,yi);
    cs.E_z  = dot(cs.p_i,zi);
    cs.E_r   = sqrt(cs.E_x^2 + cs.E_y^2 + cs.E_z^2);
    
    % Calculate deltas:
    % xi = E_x + dx; thus dx = xi - E_x 
    cs.dx = xi - cs.E_x;
    cs.dy = yi - cs.E_y;
    cs.dz = zi - cs.E_z;
    cs.dr = sqrt(cs.dx.^2 + cs.dy.^2 + cs.dz.^2);

    % Standard deviation: (Only used for energy conservation check)
    cs.sigma_x = sqrt(dot(cs.p_i,cs.dx.^2));
    cs.sigma_y = sqrt(dot(cs.p_i,cs.dy.^2));
    cs.sigma_z = sqrt(dot(cs.p_i,cs.dz.^2));
    cs.sigma_r = sqrt(cs.sigma_x^2 + cs.sigma_y^2 + cs.sigma_z^2);

    if skew_flag
        % Skewness of delta vector: (Used to setup merge cell coordinate system)
        n = 2;
        cs.mu3_dx = dot(cs.p_i,cs.dx.*cs.dr.^n); 
        cs.mu3_dy = dot(cs.p_i,cs.dy.*cs.dr.^n); 
        cs.mu3_dz = dot(cs.p_i,cs.dz.*cs.dr.^n); 
        cs.mu3_dr = sqrt(cs.mu3_dx^2 + cs.mu3_dy^2 + cs.mu3_dz^2);
    end
end

function [e_prime] = basis_prime(cs) 
    % Unit vectors in standard coordinate:
    x_hat = [1,0,0]';
    y_hat = [0,1,0]';
    z_hat = [0,0,1]';
    
    % Unit vectors in merge-cell coordinate:
    uvec_type = 2;
    switch uvec_type
        case 1
            x_hat_prime = [+cs.E_x,+cs.E_y,+cs.E_z]'/cs.E_r;
        case 2
            x_hat_prime = [+cs.mu3_dx,+cs.mu3_dy,+cs.mu3_dz]'/cs.mu3_dr;
    end
    z_hat_prime = cross(x_hat,x_hat_prime);
    y_hat_prime = cross(z_hat_prime,x_hat_prime);
    
    % Coordinate basis matrix:
    e_prime = [x_hat_prime,y_hat_prime,z_hat_prime];
end

function [xi_prime,yi_prime,zi_prime] = data_prime(e_prime,xi,yi,zi)
    % Rotation vector:
    R = transpose(e_prime);
    
    % Values in merge cell coord system:
    xi_prime = R(1,1)*xi + R(1,2)*yi + R(1,3)*zi;
    yi_prime = R(2,1)*xi + R(2,2)*yi + R(2,3)*zi;
    zi_prime = R(3,1)*xi + R(3,2)*yi + R(3,3)*zi;
end

function [wi_M,xi_M,yi_M,zi_M] = down_sample_vranic(e_prime,cs_prime,cs)

    % Calculate the standard deviations:
    sigma_x_prime = cs_prime.sigma_x;
    sigma_y_prime = cs_prime.sigma_y;
    sigma_z_prime = cs_prime.sigma_z;
    
    % Calculate deltas of new set M:
    M = 6;
    dx_prime_M(1,1) = + sqrt(M/2)*sigma_x_prime;
    dx_prime_M(2,1) = - sqrt(M/2)*sigma_x_prime;
    dx_prime_M(3,1) = 0;
    dx_prime_M(4,1) = 0;
    dx_prime_M(5,1) = 0;
    dx_prime_M(6,1) = 0;
    
    dy_prime_M(1,1) = 0;
    dy_prime_M(2,1) = 0;
    dy_prime_M(3,1) = + sqrt(M/2)*sigma_y_prime;
    dy_prime_M(4,1) = - sqrt(M/2)*sigma_y_prime;
    dy_prime_M(5,1) = 0;
    dy_prime_M(6,1) = 0;
    
    dz_prime_M(1,1) = 0;
    dz_prime_M(2,1) = 0;
    dz_prime_M(3,1) = 0;
    dz_prime_M(4,1) = 0;
    dz_prime_M(5,1) = + sqrt(M/2)*sigma_z_prime;
    dz_prime_M(6,1) = - sqrt(M/2)*sigma_z_prime;
    
    % Rotation vector:
    R = transpose(e_prime);

    % Convert deltas of M set to standard frame:
    dx_M = R(1,1)*dx_prime_M + R(2,1)*dy_prime_M + R(3,1)*dz_prime_M;
    dy_M = R(1,2)*dx_prime_M + R(2,2)*dy_prime_M + R(3,2)*dz_prime_M;
    dz_M = R(1,3)*dx_prime_M + R(2,3)*dy_prime_M + R(3,3)*dz_prime_M;
    
    % Get vectors for M set:
    xi_M = cs.E_x + dx_M;
    yi_M = cs.E_y + dy_M;
    zi_M = cs.E_z + dz_M;
    wi_M = ones(size(xi_M))*cs.w_t/M;
end

function [wi_M,xi_M,yi_M,zi_M] = down_sample(wi,xi,yi,zi)
    % Calculate merge cell statistics:
    cs = cell_stats(wi,xi,yi,zi,1);

    % Calculate coordinate system for merge-cell:
    e_prime = basis_prime(cs);

    % Express data in merge-cell basis:
    [xi_prime,yi_prime,zi_prime] = data_prime(e_prime,xi,yi,zi);

    % Calculate merge cell statistics:
    cs_prime = cell_stats(wi,xi_prime,yi_prime,zi_prime,0);

    % Calculate set M:
    [wi_M,xi_M,yi_M,zi_M] = down_sample_vranic(e_prime,cs_prime,cs);
end