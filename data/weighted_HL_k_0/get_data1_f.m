function [f,grad_f_r,grad_f_z] = get_data1_f()
%GET_U_DATA1 Get u, grad_rz^1(f) wrt r and grad_rz^1(f) wrt z for data1
%   data1
%   u = (1/3)r^3 - (1/2)r^2
%   f = -(8/3)r + (5/2)

f = @(r,z) -(8./3).*r + (3./2);
grad_f_r =@(r,z) -8./3;
grad_f_z =@(r,z) 0;
end

