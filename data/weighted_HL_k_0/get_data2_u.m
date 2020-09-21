function [u,grad_u_r,grad_u_z] = get_u_data2()
%GET_U_DATA1 Get u, grad_rz^1(u) wrt r and grad_rz^1(u) wrt z for data1
%   data1: f =

u = @(r,z) r.^3 - (3./2).*r.^2;
grad_u_r =@(r,z) 3.*r.^2 - 3.*r;
grad_u_z =@(r,z) 0;
end

