function [u_vec_r,u_vec_th,u_vec_z,s,f_vec_r,f_vec_th,f_vec_z] = get_data_3(n)
%GET_DATA_3 Get u vector and s equation exact solution
%   data3
%   u = [ r^4 - r^3
%         0
%         (z^2 - z)((1/3)r^3 - (1/2)r^2) ]
%   s = -(5r^3 - 4r^2 + (2z - 1)((1/3)r^3 - (1/2)r^2))
%   f = [ n^2(r^2 - r) + (2z - 1)(r^2 - r) - 15r^2 + 8r - (2z - 1)(r^2 - r)
%         -n(2z - 1)((1/3)r^2 - (1/2)r) - n(3r^2 - 2r) + n(5r^2 - 4r + (2z - 1)((1/3)r^2 - (1/2)r))
%         n^2(z^2 - z)((1/3)r - (1/2)) - (z^2 - z)(r - 1) - (z^2 - z)(2r - 1) - (2/3)r^3 + r^2; ]
% Author: Nicole Stock
% Date: Spring 2021

u_vec_r = @(r,z) r.^4 - r.^3;
u_vec_th = @(r,z) 0;
u_vec_z = @(r,z) (z.^2 - z).*((1./3).*r.^3 - (1./2).*r.^2);
s = @(r,z) -(5.*r.^3 - 4.*r.^2 + (2.*z - 1).*((1./3).*r.^3 - (1./2).*r.^2));
f_vec_r = @(r,z) n.^2.*(r.^2 - r) + (2.*z - 1).*(r.^2 - r) ...
    - 15.*r.^2 + 8.*r - (2.*z - 1).*(r.^2 - r);
f_vec_th = @(r,z) -n.*(2.*z - 1).*((1./3).*r.^2 - (1./2).*r) ...
    - n.*(3.*r.^2 - 2.*r) ...
    + n.*(5.*r.^2 - 4.*r + (2.*z - 1).*((1./3).*r.^2 - (1./2).*r));
f_vec_z = @(r,z) n.^2.*(z.^2 - z).*((1./3).*r - (1./2)) ...
    - (z.^2 - z).*(r - 1) - (z.^2 - z).*(2.*r - 1) ...
    - (2./3).*r.^3 + r.^2;
end