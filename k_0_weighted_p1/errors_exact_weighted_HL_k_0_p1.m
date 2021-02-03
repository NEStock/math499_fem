function [err] = errors_exact_weighted_HL_k_0_p1(p,t,basis,u_h,n,u,grad_u_r,grad_u_z)                                                             
% ERRORS_EXACT_WEIGHTED_HL_L_0_p1 - Calculate the errors of our solution u_h
% compared to the exact solution u.
%   Hodge Laplacian k = 0 case, P1
%
% Syntax:
%     [err,grad_err,max_err] = errors_exact_weighted_HL_k_0_p1(p,e,t,u_h_km1,u_h_k)
%
% Inputs:
%     p - a 2xNumNodes matrix representing nodal coordinates.
%     t - a 4xNumTriangles matrix representing the element connectivity in 
%         terms of node IDs. The end row of T represents the geometry face ID 
%         to which the element belongs.
%     basis - a 6x6xNumTriangles matrix representing piece-wise basis 
%         functions for each node in each triangle. basis(i,:,k) represents 
%         the pieceiwise basis function for the ith node in triangle T.
%     u_h - approximtate solution for u
%     n - Hodge Laplacian on Axisymmetrix Domain and its discretization
%     weight
%     u - exact solution
%     grad_u_r - gradient(u) w.r.t r
%     grad_u_z - gradient(u) w.r.t z
%
% Outputs:
%    err - L2 error
%    grad_err - L2 gradient error
%    max_err - max error
%
% Author: Nicole Stock
% Date: Fall 2020

[~,triangles] = size(t);

integral = 0;

for T = 1:triangles
    
    % get coordinates of triangle T
    coordinates = zeros(3,2);
    for i = 1:3
        node = t(i,T);
        % get x,y coordinates
        coordinates(i,:) = p(:,node);
    end
    
    [R,Z,Wr,Wz] = triquad(7, coordinates);

    I = basis(:,1,T);
    J = basis(:,2,T);
    K = basis(:,3,T);
    
    ai = I(1); aj = J(1); ak = K(1);
    bi = I(2); bj = J(2); bk = K(2); 
    ci = I(3); cj = J(3); ck = K(3);
        
    n1 = t(1,T);
    n2 = t(2,T);
    n3 = t(3,T);
    
    % find L2 Error for u_h_km1 & u_h_k
    approx =@(r,z) u_h(n1).*(1./n).*(ai.*r.^2 + bi.*r.*z + ci.*r) ...
        + u_h(n2).*(1./n).*(aj.*r.^2 + bj.*r.*z + cj.*r) ...
        + u_h(n3).*(1./n).*(ak.*r.^2 + bk.*r.*z + ck.*r);
    
    integrand =@(r,z) ((u(r,z) - approx(r,z)).^2).*r;
    
    integral = integral + Wr'*feval(integrand,R,Z)*Wz;
    
end

err = sqrt(integral);

% end