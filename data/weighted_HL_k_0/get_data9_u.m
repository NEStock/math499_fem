function [u,grad_u_r,grad_u_z] = get_data9_u()
%GET_U_DATA1 Get u, grad_rz^1(u) wrt r and grad_rz^1(u) wrt z for data1
%   data3
%   u = r^(5/3)

u = @(r,z) r.^(5./3);
grad_u_r =@(r,z) (5./3).*r.^(2./3);
grad_u_z =@(r,z) 0;
end
