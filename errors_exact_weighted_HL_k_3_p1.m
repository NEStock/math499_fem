function [err_z,err_p] = errors_exact_weighted_HL_k_3_p1(p,t,ed,t_ed,basis_edges,basis_triangles,z_h,z_vec_r,z_vec_th,z_vec_z,p_h,p_exact,n)
% ERRORS_EXACT_WEIGHTED_HL_K_3_P1 - Calculate the errors of our solution x
% compared to the exact solution u.
%
% Syntax:
%     [err,grad_err,max_err] = 
%         errors_exact_weighted_HL_k_3_p1(p,t,t_ed,basis,x,u_vec_r,u_vec_z)
% Inputs:
%     p - a 2xNumNodes matrix representing nodal coordinates.
%     t - a 4xNumTriangles matrix representing the element connectivity in 
%         terms of node IDs. The end row of T represents the geometry face ID 
%         to which the element belongs.
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
%     z_h - approximated solution of z
%     z_vec_r - exact solution p vector r component
%     z_vec_r - exact solution p vector theta component
%     z_vec_z - exact solution p vector z component
%     p_h - approximated solution of p
%     p_exact - exact solution p
%     n - Hodge Laplacian on Axisymmetrix Domain and its discretization
%     weight
%
% Outputs:
%    err_z - L2 error for z approximation
%    err_p - L2 error for p approximation
%
% Author: Nicole Stock
% Date: Fall 2020

%err_z = errors_exact_weighted_fourier_rt(p,t,ed,t_ed,basis_edges,basis_triangles,z_h,z_vec_r,z_vec_th,z_vec_z,n);

[~,triangles] = size(t_ed);
integral_z = 0;
integral_p = 0;

for T = 1:triangles
    
    % get coordinates of triangle T
    coordinates = zeros(3,2);
    for i = 1:3
        node = t(i,T);
        % get x,y coordinates
        coordinates(i,:) = p(:,node);
    end
        
    [X,Y,Wx,Wy] = triquad(7, coordinates);

    % find L2 Error for z
    Ie = basis_edges(:,1,T);
    Je = basis_edges(:,2,T);
    Ke = basis_edges(:,3,T);

    Ai = Ie(1); Bi = Ie(2); Ci = Ie(3);
    Aj = Je(1); Bj = Je(2); Cj = Je(3);
    Ak = Ke(1); Bk = Ke(2); Ck = Ke(3);

    ei = t_ed(1,T) + triangles;
    ej = t_ed(2,T) + triangles;
    ek = t_ed(3,T) + triangles;

    Di = basis_triangles(1,T);
     
    ti = T; % + edges;
    
    approx_r =@(r,z) z_h(ei).*(Bi + Ai.*r) + z_h(ej).*(Bj + Aj.*r) ...
        + z_h(ek).*(Bk + Ak.*r);
    % r component of triangle basis functions is 0
    
    approx_th =@(r,z) z_h(ei).*(1./n).*(Bi + Ai.*r) ...
        + z_h(ej).*(1./n).*(Bj + Aj.*r) + z_h(ek).*(1./n).*(Bk + Ak.*r) ...
        + z_h(ti).*(1./n).*Di.*r;
    
    approx_z =@(r,z) z_h(ei).*(Ci + Ai.*z) + z_h(ej).*(Cj + Aj.*z) ...
        + z_h(ek).*(Ck + Ak.*z);
    % z component of triangle basis functions is 0
    
    % find L2 Error
    integrand =@(r,z) ((z_vec_r(r,z) - approx_r(r,z)).^2 ...
        + (z_vec_th(r,z) - approx_th(r,z)).^2 ...
        + (z_vec_z(r,z) - approx_z(r,z)).^2).*r;
    
    integral_z = integral_z + Wx'*feval(integrand,X,Y)*Wy;
    
    % find L2 Error for p
    integrand =@(r,z) ((p_exact(r,z) - p_h(T)).^2).*r;
    
    integral_p = integral_p + Wx'*feval(integrand,X,Y)*Wy;    

    
end

err_z = sqrt(integral_z);
err_p = sqrt(integral_p);
% end