%--------------------------------------------------------------------------
% Beam gain over distance.
%--------------------------------------------------------------------------
clear;clc;close all;
addpath('UCAfunction\');
%% Cylindrical array settings
Nt     = 512;
fc     = 30e9;
c      = 3e8;
lambda = c/fc;
d      = lambda/2;

r_radius  = Nt*d/2/pi; % half-wavelength spacing
D_fres    = (2*r_radius)^(4/3)/2/(lambda)^(1/3);
plot_line = 0:0.01:1;
%% plot the beam gain over distance
theta_0 = pi/2;
r_0     = 20;
phi_0   = 0;

num_point = 100000;
num_r     = 5000;
r_list  = linspace(D_fres, 100, num_r);
M_list  = [150,300,600];
num_m   = length(M_list);
g_r     = zeros(num_m, num_r);
g_r_est = zeros(num_m, num_r);

for i_m = 1:num_m
    M = M_list(i_m);
    at = generate_cylindrical(fc, Nt, M, r_radius, d, r_0, theta_0, phi_0);
    for idx_r = 1:num_r
        r_iter = r_list(idx_r);
        h0 = generate_cylindrical(fc, Nt, M, r_radius, d, r_iter, theta_0, phi_0);
        g_r(i_m, idx_r)=  abs(at*h0');
        
        beta = sqrt(M^2*d^2/2/lambda*abs(1/r_0-1/r_iter));
        [f_sin, f_cos] = fresnel(beta, num_point);
        G = abs(f_cos+1i*f_sin)/beta;
        xi = pi*r_radius^2/2/lambda*(1/r_0-1/r_iter);
        g_r_est(i_m, idx_r) = abs(besselj(0, xi)*G);
    end
end
figure;
hold on;
plot(r_list, g_r(1,:),'-','Color',1.00*[0.80 0 1],'LineWidth',1.5);
plot(r_list, g_r(2,:),'-','Color',0.85*[1 0 0],'LineWidth',1.5);
plot(r_list, g_r(3,:),'-','Color',0.80*[0 1 0],'LineWidth',1.5);
plot(r_list, g_r_est,'k--','LineWidth',1.5);
plot(D_fres*ones(1, length(plot_line)), plot_line,'-.','Color',0.75*[1 1 1],'LineWidth',1.2);
grid on;
box on;
xlabel('Distance (m)');
ylabel('Beamforming Gain');
legend('M=150','M=300','M=600','Bessel fitting','FontSize',16,'Location','NorthEast');
set(gca,'YTick',[0 0.2 0.4 0.6 0.8 1],'FontName','Times New Roman','FontSize',16,'GridLineStyle','-.');
