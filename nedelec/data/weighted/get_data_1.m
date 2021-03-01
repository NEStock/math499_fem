function [u_vec_r,u_vec_z,f_vec_r,f_vec_z] = get_data_1()
%GET_DATA_1 Get u vector and f vector
%   u = [ sin(pi*r) ; cos(pi*r) ]
%   F = u
% Author: Nicole Stock
% Date: Fall 2020

u_vec_r = @(r,z) sin(pi.*r);
u_vec_z = @(r,z) cos(pi.*z);

f_vec_r = @(r,z) sin(pi.*r);
f_vec_z = @(r,z) cos(pi.*z);
end