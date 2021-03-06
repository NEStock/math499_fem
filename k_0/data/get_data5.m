function [f,grad_f_r,grad_f_z,u,grad_u_r,grad_u_z] = get_data5()
%GET_U_DATA5 Get u and f function handles
%
%   u = r^(5/6)
% Author: Nicole Stock
% Date: Fall 2020

u = @(r,z) r.^(5./6);
grad_u_r =@(r,z) (5./6).*r.^(-1./6);
grad_u_z =@(r,z) 0;

f = u;
grad_f_r = grad_u_r;
grad_f_z = grad_u_z;
end

