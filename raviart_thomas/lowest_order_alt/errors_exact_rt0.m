function [err] = errors_exact_rt0(p,t,ed,t_ed,basis,x,u_vec_r,u_vec_z)
% ERRORS_EXACT_RT0 - Calculate the errors of our solution x
% compared to the exact solution u.
%
% Syntax:
%     [err,grad_err,max_err] = 
%         errors_exact_weighted_rt0(p,t,t_ed,basis,x,u_vec_r,u_vec_z)
% Inputs:
%     p - a 2xNumNodes matrix representing nodal coordinates.
%     t - a 4xNumTriangles matrix representing the element connectivity in 
%         terms of node IDs. The end row of T represents the geometry face ID 
%         to which the element belongs
%     t_ed - a 3xNumTriangles matrix representing the which edges
%         correspond to which triangles. t_ed(i,T) represents the ith edge
%         in triangle T.
%     basis - a 3x2xNumTriangles matrix representing basis functions for
%         each node in each triangle. basis(i,:,T) represents the basis 
%         function for the ith node in triangle T. 
%     x - approximated solution
%     u_vec_r - exact solution vector r component
%     u_vec_z - exact solution vector z component
%
% Outputs:
%    err - L2 error
%
% Author: Nicole Stock
% Date: Fall 2020

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

    e1 = t_ed(1,T);
    e2 = t_ed(2,T);
    e3 = t_ed(3,T);
    
    approx_r =@(r,z) x(e1).*basis{1,1,T}(r,z) + x(e2).*basis{2,1,T}(r,z) ...
        + x(e3).*basis{3,1,T}(r,z);
    
    approx_z =@(r,z) x(e1).*basis{1,2,T}(r,z) + x(e2).*basis{2,2,T}(r,z) ...
        + x(e3).*basis{3,2,T}(r,z);
    
    % find L2 Error
    integrand =@(r,z) ((u_vec_r(r,z) - approx_r(r,z)).^2 ...
        + (u_vec_z(r,z) - approx_z(r,z)).^2);
    
    integral = integral + Wx'*feval(integrand,X,Y)*Wy;
end

err = sqrt(integral);
% end