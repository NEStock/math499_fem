function [u_vec_r,u_vec_z,f_vec_r,f_vec_z] = get_data_1()
%GET_DATA_1 Get u vector and f vector
%   u = [ (1/3)r^3z^2 - (1/2)r^2z^2 
%         -(1/2)r^2z^2 + (1/2)rz^2  ]
%   F = [ (1/3)r^3z^2 - (1/2)r^2z^2 - 2rz^2 + z^2 + 2rz - z
%         -(1/2)r^2z^2 + (1/2)rz^2 - 2r^2z + 2rz + r^2 - r   ]
% Author: Nicole Stock
% Date: Fall 2020

u_vec_r = @(r,z) (1./3).*(r.^3).*(z.^2) - (1./2).*(r.^2).*(z.^2);
u_vec_z = @(r,z) -(1./2).*(r.^2).*(z.^2) + (1./2).*r.*z.^2;

f_vec_r = @(r,z) (1./3).*(r.^3).*(z.^2) - (1./2).*(r.^2).*(z.^2) ...
    - 2.*r.*z.^2 + z.^2 + 2.*r.*z - z;
f_vec_z = @(r,z) -(1./2).*(r.^2).*(z.^2) + (1./2).*r.*z.^2 ...
    - 2.*z.*r.^2 + r.^2 + 2.*r.*z - r;
end