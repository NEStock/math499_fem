function [u_vec_r,u_vec_th,u_vec_z,s_vec_r,s_vec_th,s_vec_z,f_vec_r,f_vec_th,f_vec_z] = get_data_1(n)
%GET_DATA1 Get u vector and s equation exact solution
%   data1
%   u = [ 3r(r-1)
%         -3r^2 + 2r
%         0          ]
%   s = [ 0
%         0
%         r^2(r-1) ]
%   f = [ 3n^2 - (3n^2)/r + 6n - (4n)/r - 9
%         -3n^2 + (2n^2)/r - 6n + (6n)/r + 9
%         0      ]
% Author: Nicole Stock
% Date: Fall 2020

s_vec_r = @(r,z) 3.*r.*(r-1); % 3.*r.^2 - 3.*r
s_vec_th = @(r,z) -3.*r.^2 + 2.*r;
s_vec_z = @(r,z) 0;
u_vec_r = @(r,z) 0;
u_vec_th = @(r,z) 0;
u_vec_z = @(r,z) (r.^2).*(r-1); % r.^3 - r.^2
f_vec_r = @(r,z) 0;
f_vec_th = @(r,z) 0;
f_vec_z = @(r,z) (n.^2).*r - (n.^2) - 9.*r + 4;
end