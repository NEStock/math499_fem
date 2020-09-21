function [f,grad_f_r,grad_f_z] = get_f_data2()
%GET_U_DATA1 Get f, grad_rz^1(f) wrt r and grad_rz^1(f) wrt z for data1
%   data2

f = @(r,z) -8.*r + (9./2);
grad_f_r =@(r,z) -8;
grad_f_z =@(r,z) 0;
end

