function [f,grad_f_r,grad_f_z] = get_f_data1()
%GET_U_DATA1 Get u, grad(u) wrt r and grad(u) wrt z for data1
%   data1: f = -(8/3)r + (5/2)

f = @(r,z) -(8./3).*r + (3./2);
grad_f_r =@(r,z) -8./3;
grad_f_z =@(r,z) 0;
end

