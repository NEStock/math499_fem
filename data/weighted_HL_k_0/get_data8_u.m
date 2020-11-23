function [u,grad_u_r,grad_u_z] = get_data8_u()
%GET_U_DATA1 Get u, grad_rz^1(u) wrt r and grad_rz^1(u) wrt z for data1
%   data3
%   u = r^(4/3)

u = @(r,z) r.^(4./3);
grad_u_r =@(r,z) (4./3).*r.^(1./3);
grad_u_z =@(r,z) 0;
end
