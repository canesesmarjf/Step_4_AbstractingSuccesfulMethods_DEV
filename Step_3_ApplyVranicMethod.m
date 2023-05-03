% Step 3: Use data produced by main_4
clear all
close all
clc

saveFig  = 0;
saveData = 0;

% Select data:
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

% Import dataset produced by PICOS++:
% =========================================================================
% Save data to tables:
target_file = "./input_files/";
x_p = readmatrix(target_file + "Step_1_x_p.csv");
v_p = readmatrix(target_file + "Step_1_v_p.csv");
a_p = readmatrix(target_file + "Step_1_a_p.csv");

% Test of 3D vranic method:
% =========================================================================
% Choose dataset to apply method:
file_id = ["a", "b" , "c", "d" , "e", "f"];

for ii = 1:numel(file_id)
% Get index data:
ip = csvread(target_file + "ip_main_4" + file_id{ii} + ".csv") + 1;

% Select subset N:
N(ii) = round(numel(ip)/1);

% Plot merge cell data:
hold on
plot3(x_p(ip(1:N(ii))) ,v_p(ip(1:N(ii)),1) ,v_p(ip(1:N(ii)),2) ,'r.','MarkerSize',5)

% Assemble N set data:
wi{ii} = a_p(ip(1:N(ii)));
xi{ii} = x_p(ip(1:N(ii)));
yi{ii} = v_p(ip(1:N(ii)),1);
zi{ii} = v_p(ip(1:N(ii)),2);

% Calculate M set data:
[wi_M{ii},xi_M{ii},yi_M{ii},zi_M{ii}] = down_sample(wi{ii},xi{ii},yi{ii},zi{ii});

% Calculate merge cell statistics:
cs{ii} = cell_stats(wi{ii},xi{ii},yi{ii},zi{ii},1);
cs_M{ii} = cell_stats(wi_M{ii},xi_M{ii},yi_M{ii},zi_M{ii},0);

% Plot data:
hold on
plot3(xi_M{ii},yi_M{ii},zi_M{ii},'g.','MarkerSize',5)
plot3(xi_M{ii},yi_M{ii},zi_M{ii},'go','MarkerSize',5)

% Check conservation of mass, momentum and energy:
disp('Number of particles:');
disp("Set N = " + num2str(N(ii)))
disp("Set N = " + num2str(6))

disp('Conservation of mass:');
disp("Set N, w_t = " + num2str(cs{ii}.w_t))
disp("Set M, w_t = " + num2str(cs_M{ii}.w_t))

disp('Conservation of momentum:');
disp("Set N, E_r = " + num2str(cs{ii}.E_r))
disp("Set M, E_r = " + num2str(cs_M{ii}.E_r))

disp('Conservation of energy:');
disp("Set N, sigma_r = " + num2str(cs{ii}.sigma_r))
disp("Set M, sigma_r = " + num2str(cs_M{ii}.sigma_r))
end

%% Plot merge cell with coordinate system:
ii = 3;

figure('color','w');

% Plot merge-cell data:
dx = xi{ii} - cs{ii}.E_x;
dy = yi{ii} - cs{ii}.E_y;
dz = zi{ii} - cs{ii}.E_z;
plot3(dx,dy,dz,'k.')

% Plot coordinate system:
e_prime = basis_prime(cs{ii});
basis_plot(e_prime(:,1),e_prime(:,2),e_prime(:,3),[1,1,1]*4e4)

% Plot M set particles:
plot3(xi_M{ii} - cs_M{ii}.E_x ,yi_M{ii} - cs_M{ii}.E_y, zi_M{ii} - cs_M{ii}.E_z,...
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
    cs.dx = cs.E_x - xi;
    cs.dy = cs.E_y - yi;
    cs.dz = cs.E_z - zi;
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