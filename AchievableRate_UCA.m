%--------------------------------------------------------------------------
% Achievable Rate.
%--------------------------------------------------------------------------
clear;clc; close all;
addpath('UCAfunction\');
%% parameters
Nt_ULA = 512; % element for UCA
fc = 30e9;    % 
c = 3e8;
lambda = c/fc;
d = lambda/2;

% find corresponding point for threshold
threshold = 0.5;
num_bessel = 1000;
r_list = linspace(0, 4, num_bessel);
g_r = abs(besselj(0, r_list));
[~,idx] = min(abs(g_r-threshold));
xi_delta = r_list(idx);

xi_delta_zero = 2.4048;

r_radius = Nt_ULA*d/2/pi; % half-wavelength spacing
D_frau = 2*(2*r_radius)^2/lambda;
D_fres = (2*r_radius)^(4/3)/2/(lambda)^(1/3);
Nt = 256;

%% codebook in the azimuth domain
delta_phi = asin(xi_delta*lambda/4/pi/r_radius)*2;
num_phi_beam = 2*pi/delta_phi;

delta_phi_zero = asin(xi_delta_zero*lambda/4/pi/r_radius)*2;
num_phi_zero = 2*pi/delta_phi_zero;

%% verify the codebook
num_beam_phi = floor(num_phi_beam);
phi_list_tmp = (0:num_beam_phi-1)*delta_phi;
beam_phi_list = UCA_far_beam(fc, Nt, r_radius, phi_list_tmp);
correlation_phi = beam_phi_list*beam_phi_list';
correlation_phi = abs(correlation_phi-diag(diag(correlation_phi)));
max_corr_phi = max(correlation_phi);

%% verify the non-sparse property in the phi domain
phi_idx = 1;
beam_phi_fixed = phi_list_tmp(phi_idx);
r_near = 30;
r_far  = 10000;
h_near = UCA_generate(fc, Nt, r_radius, r_near, beam_phi_fixed);
num_phi = 1e5;
g_near = zeros(1, num_phi);
g_far = zeros(1, num_phi);
phi_list = linspace(-pi/20, pi/20, num_phi);
for idx_phi = 1:num_phi
    phi_iter = phi_list(idx_phi);
    h0_near = UCA_generate(fc, Nt, r_radius, r_near, phi_iter);
    h0_far  = UCA_generate(fc, Nt, r_radius, r_far, phi_iter);
    g_near(idx_phi)=  abs(h_near*h0_near');
    g_far(idx_phi) =  abs(h_near*h0_far');
end
figure;
hold on;box on;
plot(phi_list, g_near,'-','Color',0.80*[1 0 0],'LineWidth',1.2);
plot(phi_list, g_far, '-','Color',0.80*[0 0 0],'LineWidth',1.2);
xlabel('Angle \phi');
ylabel('Beamforming Gain');
legend('r = 30 m','r = 10000 m','FontSize',16);
set(gca,'YTick',[0 0.2 0.4 0.6 0.8 1],'FontName','Times New Roman','FontSize',16);

% plot the multi-beam in azimuth domain
num_phi_plot = 1e4;
plot_idx = 10;
phi_list_plot = linspace(0, pi/80, num_phi_plot);
gain_phi = zeros(plot_idx, num_phi_plot);
for idx_phi_plot = 1:num_phi_plot
    phi_iter = phi_list_plot(idx_phi_plot);
    h0 = UCA_generate(fc, Nt, r_radius, r_far, phi_iter);
    gain_phi(:, idx_phi_plot)=  abs(beam_phi_list(1:plot_idx,:)*h0');
end

figure;
plot(phi_list_plot, gain_phi,'LineWidth',1.5);
xlabel('\phi', 'FontSize', 13);
ylabel('Beamforming Gain');


%% codebook in the distance domain
r_min = 4;
delta_distance = xi_delta*2*lambda/r_radius^2/pi;
num_distance_beam = 1/r_min/delta_distance;

delta_distance_zero = xi_delta_zero*2*lambda/r_radius^2/pi;
num_distance_zero = 1/r_min/delta_distance_zero;

%% verify the codebook
num_beam_distance = floor(num_distance_beam);
dis_list_tmp = zeros(1, num_beam_distance);
for i_dis = 1:num_beam_distance
    if i_dis == 1
        dis_list_tmp(i_dis) = 10000;
    else
        dis_list_tmp(i_dis) = 1/delta_distance/(i_dis-1);
    end
end
phi_0 = 0;
phi_list_for_dis = phi_0*ones(1, num_beam_distance);
beam_dis_list = UCA_generate(fc, Nt, r_radius, dis_list_tmp, phi_list_for_dis);
correlation_dis = beam_dis_list*beam_dis_list';
correlation_dis = abs(correlation_dis-diag(diag(correlation_dis)));
max_corr_dis = max(correlation_dis);

% plot the multi-beam in the distance domain
num_r = 10000;
r_list = linspace(r_min, 150, num_r);
gain_beams = zeros(num_beam_distance, num_r);
for idx_r = 1:num_r
    r_iter = r_list(idx_r);
    h0 = UCA_generate(fc, Nt, r_radius, r_iter, phi_0);
    gain_beams(:, idx_r) = abs(beam_dis_list*h0');
end

figure;
plot(r_list, gain_beams,'LineWidth',1.5);
xlabel('Distance (m)');
ylabel('Beamforming Gain');

%% polar domain codebook
num_code = num_beam_distance*num_beam_phi;
polar_codebook = zeros(num_code, Nt);
cnt = 1;
for i_dis = 1:num_beam_distance
    r_value = dis_list_tmp(i_dis);
    for j_phi = 1:num_beam_phi
        phi_value = phi_list_tmp(j_phi);
        polar_codebook(cnt, :) = UCA_generate(fc, Nt, r_radius, r_value, phi_value);
        cnt = cnt + 1;
    end
end

% far codebook
far_codebook = zeros(num_beam_phi, Nt);
for j_phi = 1:num_beam_phi
    phi_value = phi_list_tmp(j_phi);
    far_codebook(j_phi, :) = UCA_generate(fc, Nt, r_radius, r_far, phi_value);
end

%% channel with random locations
N_user = 100;
r_min = 5;
r_max = 20;

cnt = 1;
x_list = zeros(1,N_user);
y_list = zeros(1,N_user);
while(cnt <= N_user)
    x_1 = abs(rand()*r_max);
    y_1 = (rand()-0.5)*2*r_max;
    [theta_1, r_1]  = cart2pol(x_1, y_1);
    if (x_1^2+y_1^2<r_max^2) && (x_1^2+y_1^2>r_min^2)
        x_list(cnt) = x_1;
        y_list(cnt) = y_1;
        cnt = cnt+1;
    end
end
[theta_user, r_user]  = cart2pol(x_list, y_list);
H_user = UCA_generate(fc, Nt, r_radius, r_user, theta_user);

SNR_dB = -10:5:30;
SNR_linear = 10.^(SNR_dB/10.);
sum_rate_near = zeros(1, length(SNR_linear));
sum_rate_ideal = zeros(1, length(SNR_linear));
sum_rate_far = zeros(1, length(SNR_linear));
for i_snr=1:length(SNR_linear)
    SNR = SNR_linear(i_snr);
    gain = max(abs(H_user*polar_codebook'),[],2);
    gain_far = max(abs(H_user*far_codebook'),[],2);
    sum_rate_near(i_snr) = mean(log2(1+SNR*gain));
    sum_rate_far(i_snr) = mean(log2(1+SNR*gain_far));
    sum_rate_ideal(i_snr) = mean(log2(1+SNR));
end
figure;
hold on;grid on;box on;
plot(SNR_dB, sum_rate_ideal,'k--*','LineWidth',1.2,'MarkerSize',9);
plot(SNR_dB, sum_rate_near, 'r-o','LineWidth',1.2,'MarkerSize',9);
plot(SNR_dB, sum_rate_far,  'b-.s','LineWidth',1.2,'MarkerSize',10);
xlabel('SNR (dB)');
ylabel('Achievable Rate (bps/Hz)');
legend('Optimal','Near-field','Far-field','Location','NorthWest','FontSize',16);
set(gca,'FontName','Times New Roman','FontSize',16);
