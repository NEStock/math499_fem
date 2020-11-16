function [err] = errors_exact_weighted_fourier_rt(p,t,ed,t_ed,basis_edges,basis_triangles,x,u_vec_r,u_vec_th,u_vec_z,n)
% ERRORS_EXACT_WEIGHTED_FOURIER_RT - Calculate the errors of our solution x
% compared to the exact solution u.
%
% Syntax:
%     [err,grad_err,max_err] = 
%         errors_exact_weighted_fourier_rt(p,t,ed,t_ed,basis_edges,basis_triangles,x,u_vec_r,u_vec_th,u_vec_z,n)
% Inputs:
%     p - a 2xNumNodes matrix representing nodal coordinates.
%     t - a 4xNumTriangles matrix representing the element connectivity in 
%         terms of node IDs. The end row of T represents the geometry face 
%         ID to which the element belongs.
%     ed - a 2xNumEdges matrix representing each edge as a row with
%         starting node in column 1 and the ending node in column 2.
%     t_ed - a 3xNumTriangles matrix representing the which edges
%         correspond to which triangles. t_ed(i,T) represents the ith edge
%         in triangle T.
%     basis_edges - a matrix representing piece-wise basis functions for 
%         each edge in each triangle. basis(i,:,T) represents the 
%         pieceiwise basis function for the ith edge in triangle T.
%     basis_triangles - a vector representing piece-wise basis functions
%         for edge triangle. basis(1,T) represents the piecewise basis
%         function for the Tth triangle.
%     x - approximated solution
%     u_vec_r - exact solution vector r component
%     u_vec_th - exact solution vector theta component
%     u_vec_z - exact solution vector z component
%     n - Hodge Laplacian on Axisymmetrix Domain and its discretization
%     weight
%
% Outputs:
%    err - L2 error
%
% Author: Nicole Stock
% Date: Fall 2020

[~,triangles] = size(t_ed);
[edges,~] = size(ed);
integral = 0;

for T = 1:triangles
    
    % get coordinates of triangle T
    coordinates = zeros(3,2);
    for i = 1:3
        node = t(i,T);
        % get x,y coordinates
        coordinates(i,:) = p(:,node);
    end
        
    [X,Y,Wx,Wy] = triquad(7, coordinates);
 
    Ie = basis_edges(:,1,T);
    Je = basis_edges(:,2,T);
    Ke = basis_edges(:,3,T);

    Ai = Ie(1); Bi = Ie(2); Ci = Ie(3);
    Aj = Je(1); Bj = Je(2); Cj = Je(3);
    Ak = Ke(1); Bk = Ke(2); Ck = Ke(3);

    ei = t_ed(1,T);
    ej = t_ed(2,T);
    ek = t_ed(3,T);

    Di = basis_triangles(1,T);
     
    ti = T + edges;
    
    approx_r =@(r,z) x(ei).*(Bi + Ai.*r) + x(ej).*(Bj + Aj.*r) ...
        + x(ek).*(Bk + Ak.*r);
    % r component of triangle basis functions is 0
    
    approx_th =@(r,z) x(ei).*(1./n).*(Bi + Ai.*r) ...
        + x(ej).*(1./n).*(Bj + Aj.*r) + x(ek).*(1./n).*(Bk + Ak.*r) ...
        + x(ti).*(1./n).*Di.*r;
    
    approx_z =@(r,z) x(ei).*(Ci + Ai.*z) + x(ej).*(Cj + Aj.*z) ...
        + x(ek).*(Ck + Ak.*z);
    % z component of triangle basis functions is 0
    
    % find L2 Error
    integrand =@(r,z) ((u_vec_r(r,z) - approx_r(r,z)).^2 ...
        + (u_vec_th(r,z) - approx_th(r,z)).^2 ...
        + (u_vec_z(r,z) - approx_z(r,z)).^2).*r;
    
    integral = integral + Wx'*feval(integrand,X,Y)*Wy;
end

err = sqrt(integral);
% end