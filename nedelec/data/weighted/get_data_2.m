function [u_vec_r,u_vec_z,f_vec_r,f_vec_z] = get_data_2()
%GET_DATA_2 Get u vector and f vector
%   u = [ (1/3)(r-1)z^3
%         (1/2)r^2z - rz ]
%   F = [ (1/3)rz^3 - (1/3)z^3 - 2rz + 2z + r - 1
%         2z^2 - (z^2)/r + (1/2)r^2z - rz - 2z + z/r ]
% Author: Nicole Stock
% Date: Fall 2020

u_vec_r =@(r,z) (1/3).*(r-1).*z.^3;
u_vec_z =@(r,z) (1/2).*r.^2.*z - r.*z;
f_vec_r =@(r,z) (1/3).*r.*z.^3 - (1/3).*z.^3 - 2.*r.*z + 2.*z + r - 1;
f_vec_z =@(r,z) 2.*z.^2 - (z.^2)./r + (1/2).*r.^2.*z - r.*z - 2.*z + z./r;
end