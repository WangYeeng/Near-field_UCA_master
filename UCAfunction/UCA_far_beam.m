function [H] = UCA_far_beam(fc, Nt, r_radius, theta_list)

c = 3e8;
lambda = c/fc;
[~, d_len] = size(theta_list);

%% UCA settings
pp = 1:1:Nt;
phi_p = (pp-1)*2*pi/Nt;

%% channel generation
H = zeros(d_len, Nt);
for i_d = 1:d_len
    theta0 = theta_list(i_d);
    dis_matrix = -r_radius*cos(theta0-phi_p);
    at = exp(-1j*2*pi/lambda*dis_matrix)/sqrt(Nt);
    w = at;
    H(i_d, :) = at;
end
end

