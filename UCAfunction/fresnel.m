function [f_sin, f_cos] = fresnel(x, num_point)

[~,num_list] = size(x);
f_sin = zeros(1, num_list);
f_cos = zeros(1, num_list);

for i_iter = 1:num_list
    x_iter = x(i_iter);
    x_list = linspace(0,x_iter,num_point);

    y_sin = sin(pi/2*x_list.^2);
    y_cos = cos(pi/2*x_list.^2);
    f_sin(i_iter) = trapz(x_list, y_sin);
    f_cos(i_iter) = trapz(x_list, y_cos);
end

end

