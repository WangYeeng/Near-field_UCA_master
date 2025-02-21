function [H] = generate_cylindrical(fc, Nt, M, r_radius, d, r_list, theta_list, phi_list)

c = 3e8;
lambda = c/fc;
[~, d_len] = size(r_list);

%% UCA settings
pp = 1:1:Nt;
phi_p = (pp-1)*2*pi/Nt;
[x_p, y_p] = pol2cart(phi_p, r_radius);

mm = (-(M-1)/2:1:(M-1)/2).';


%% channel generation
H = zeros(d_len, M*Nt);
for i_d = 1:d_len
    theta0 = theta_list(i_d);
    r0 = r_list(i_d);
    r0_xy = r0*sin(theta0);
    r0_z = r0*cos(theta0);
    phi0 = phi_list(i_d);
    [x0, y0] = pol2cart(phi0, r0_xy);
        
    dis_matrix = sqrt((x_p-x0).^2+(y_p-y0).^2+(r0_z-mm*d).^2);
    at = exp(-1i*2*pi/lambda*dis_matrix)/sqrt(Nt*M);
    H(i_d, :) = reshape(at,[],1);
end

end

