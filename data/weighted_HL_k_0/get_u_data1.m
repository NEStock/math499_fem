function [u,grad_u_r,grad_u_z] = get_u_data1()
%GET_U_DATA1 Get u, grad(u) wrt r and grad(u) wrt z for data1
%   data1: u = r^(1/2)

u = @(r,z) r.^(1/2);
grad_u_r =@(r,z) (1/2).*r^(-1/2);
grad_u_z =@(r,z) 0;
end

