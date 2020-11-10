function [u_vec_r,u_vec_z,f_vec_r,f_vec_z] = get_data_2()
%GET_DATA_1 Get u vector and f vector
%   u = [ (1/3)r^3z^2 - (1/3)r^3z 
%         -(1/3)rz^3 + (1/2)rz^2  ]
%   F = [ (1/3)r^3z^2 - (1/3)r^3z - 2rz^2 + z^2 + 2rz - z
%         -(1/3)rz^3 + (1/2)rz^2 - 2r^2z + 2rz + r^2 - r  ]

u_vec_r = @(r,z) (1./3).*(r.^3).*(z.^2) - (1./3).*z.*r.^3;
u_vec_z = @(r,z) -(1./3).*r.*z.^3 + (1./2).*r.*z.^2;

f_vec_r = @(r,z) (1./3).*(r.^3).*(z.^2) - (1./3).*z.*r.^3 ...
    - 2.*r.*z.^2 + z.^2 + 2.*r.*z - z;
f_vec_z = @(r,z) -(1./3).*r.*z.^3 + (1./2).*r.*z.^2 ...
    - 2.*z.*r.^2 + r.^2 + 2.*r.*z - r;
end