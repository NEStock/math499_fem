function [f,grad_f_r,grad_f_z] = get_data2_f()
%GET_U_DATA1 Get f, grad_rz^1(f) wrt r and grad_rz^1(f) wrt z for data1
%   data2
%   u = r^3 - (3/2)r^2
%   f = -8r + 9/2

f = @(r,z) -8.*r + (9./2);
grad_f_r =@(r,z) -8;
grad_f_z =@(r,z) 0;
end

