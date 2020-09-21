function [u,grad_u_r,grad_u_z] = get_data1_u()
%GET_U_DATA1 Get u, grad_rz^1(u) wrt r and grad_rz^1(u) wrt z for data1
%   data1
%   u = (1/3)r^3 - (1/2)r^2
%   f = -(8/3)r + (5/2)

u = @(r,z) (1./3).*r.^3 - (1./2).*r.^2;
grad_u_r =@(r,z) r.^2 - r;
grad_u_z =@(r,z) 0;
end

