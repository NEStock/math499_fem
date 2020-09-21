function [u,grad_u_r,grad_u_z] = get_data6_u()
%GET_U_DATA1 Get u, grad_rz^1(u) wrt r and grad_rz^1(u) wrt z for data1
%   data3
%   u = r*sin(z)

u = @(r,z) r.*sin(z);
grad_u_r =@(r,z) sin(z);
grad_u_z =@(r,z) -r.*cos(z);
end

