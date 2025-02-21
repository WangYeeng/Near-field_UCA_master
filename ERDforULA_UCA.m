%--------------------------------------------------------------------------
% Effective rayleigh distance (ERD) for different angles
%--------------------------------------------------------------------------
clear;clc;close all;
addpath('UCAfunction\');
%% UCA settings
Nt        = 512;                     % element for UCA     
fc        = 30e9;                    % operating frequency
c         = 3e8;
lambda    = c/fc;                    % wavelength
d         = lambda/2;
r_radius  = Nt*d/2/pi;               % half-wavelength spacing
%% effective rayleigh distance for different angles
threshold = 0.95;
num_bessel = 3000;

num_phi = 10000;
phi_list = linspace(-pi/2, pi/2, num_phi);
g = zeros(1, num_phi);
r_uca = zeros(1, num_phi);
r_ula = zeros(1, num_phi);

% decide the value for UCA
r_list = linspace(0, 1, num_bessel);
g_r = abs(besselj(0, r_list));
[~,idx] = min(abs(g_r-threshold));
epsilon_c = r_list(idx);

%decide the epsilon value for ULA
num_point = 100000;
beta = linspace(0.1, 1, num_bessel);
[f1_sin, f1_cos] = fresnel(beta, num_point);
g_r_ula = abs(f1_cos+1i*f1_sin)./beta;

[~,idx] = min(abs(g_r_ula-threshold));
epsilon_l = beta(idx);
for i_phi = 1:num_phi
    phi = phi_list(i_phi);
    r_uca(i_phi) = pi*r_radius^2/2/lambda/epsilon_c;
    r_ula(i_phi) = (2*r_radius)^2*cos(phi)^2/2/epsilon_l^2/lambda;
end

%% plot in 2D space
[x_ula,y_ula] = pol2cart(phi_list, r_ula);
[x_uca,y_uca] = pol2cart(phi_list, r_uca);

figure;
hold on;grid on;box on;
plot(x_uca,y_uca,'-.','Color',0.80*[1 0 0],'Linewidth',1.5)
plot(x_ula,y_ula,'- ','Color',0.80*[0 0 1],'Linewidth',1.5)
xlabel('x-axis');
ylabel('y-axis');
legend('ERD UCA','ERD ULA','FontSize',16);
set(gca,'YTick',[-60 -30 0 30 60],'XTick',[0 20 40 60],'FontName','Times New Roman','FontSize',16,'GridLineStyle','-.');
