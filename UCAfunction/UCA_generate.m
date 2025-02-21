function [H, H_approx, dis_matrix, dis_approx] = UCA_generate(fc, Nt, r_radius, r_list, theta_list)

c = 3e8;
lambda = c/fc;
[~, d_len] = size(r_list);

%% UCA settings
pp = 1:1:Nt;
phi_p = (pp-1)*2*pi/Nt;
[x_p, y_p] = pol2cart(phi_p, r_radius);

%% channel generation
H = zeros(d_len, Nt);
H_approx = zeros(d_len, Nt);
for i_d = 1:d_len
    r0 = r_list(i_d);
    theta0 = theta_list(i_d);
    [x0, y0] = pol2cart(theta0, r0);
    
    dis_matrix = sqrt((x_p-x0).^2+(y_p-y0).^2);
    dis_approx = r0-r_radius*cos(theta0-phi_p)+r_radius^2/2/r0*sin(theta0-phi_p).^2;
    
    at = exp(-1i*2*pi/lambda*dis_matrix)/sqrt(Nt);
    H(i_d, :) = at;
    H_approx(i_d, :) = exp(-1i*2*pi/lambda*dis_approx)/sqrt(Nt);
end

