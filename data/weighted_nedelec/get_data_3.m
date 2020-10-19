function [u_vec_r,u_vec_z,f_vec_r,f_vec_z] = get_data_3()
%GET_DATA_3 Get u vector and f vector
%   u = [ (-1/2)rz^2 + (1/2)z^2
%         (-1/2)r^2z^2 + rz^2 ]
%   F = [ (-1/2)rz^2 + (1/2)z^2 - 2rz + 2z + r - 1 
%         (-1/2)r^2z^2 + rz^2 + 2z^2 - z^2/r - 2z + z/r ]

u_vec_r =@(r,z) (-1/2).*r.*z.^2 + (1./2).*z.^2;
u_vec_z =@(r,z) (-1/2).*r.^2.*z.^2 + r.*z.^2;
f_vec_r =@(r,z) (-1/2).*r.*z.^2 + (1/2).*z.^2 - 2.*r.*z + 2.*z + r - 1;
f_vec_z =@(r,z) (-1/2).*r.^2.*z.^2 + r.*z.^2 + 2.*z.^2 - z.^2./r - 2.*z + z./r;
end