function [err] = errors_exact_weighted_nedelec(p,t,t_ed,basis,x,u_vec_r,u_vec_z)
% ERRORS_EXACT_WEIGHTED_NEDELEC - Calculate the errors of our solution u_h
% compared to the exact solution u.
%
% Syntax:
%     [err,grad_err,max_err] = errors_exact_weighted_nedelec(p,t,ed,t_ed,basis,x,u_vec_r,u_vec_z)
% Inputs:
%     p - a 2xNumNodes matrix representing nodal coordinates.
%     t - a 4xNumTriangles matrix representing the element connectivity in 
%         terms of node IDs. The end row of T represents the geometry face ID 
%         to which the element belongs
%     basis - a 3x3xNumTriangles matrix representing piece-wise basis 
%         functions for each node in each triangle. basis(i,:,k) represents 
%         the pieceiwise basis function for the ith node in triangle k. 
%
% Outputs:
%    err - L2 error
%    grad_err - L2 gradient error
%    max_err - max error
%
% Author: Nicole Stock
% Date: Spring 2020

[~,triangles] = size(t);

integral = 0;

for T = 1:triangles
    
    % get coordinates of triangle T
    coordinates = zeros(3,2);
    for n = 1:3
        node = t(n,T);
        % get x,y coordinates
        coordinates(n,:) = p(:,node);
    end
        
    [X,Y,Wx,Wy] = triquad(7, coordinates);

    I = basis(:,1,T);
    J = basis(:,2,T);
    K = basis(:,3,T);
    
    ai = I(1); bi = I(2); ci = I(3);
    aj = J(1); bj = J(2); cj = J(3);
    ak = K(1); bk = K(2); ck = K(3);
    
    ei = t_ed(1,T);
    ej = t_ed(2,T);
    ek = t_ed(3,T);
    
    approx_r =@(r,z) x(ei).*(bi - ai.*z) + x(ej).*(bj - aj.*z) ...
        + x(ek).*(bk - ak.*z);
    approx_z =@(r,z) x(ei).*(ci + ai.*r) + x(ej).*(cj + aj.*r) ...
        + x(ek).*(ck + ak.*r);
    
    % find L2 Error
    integrand =@(r,z) ((u_vec_r(r,z) - approx_r(r,z)).^2 ...
        + (u_vec_z(r,z) - approx_z(r,z)).^2).*r;
    
    integral = integral + Wx'*feval(integrand,X,Y)*Wy;
end

err = sqrt(integral);
% end