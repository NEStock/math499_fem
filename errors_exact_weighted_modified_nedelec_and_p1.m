function [err] = errors_exact_weighted_modified_nedelec_and_p1(p,t,t_ed,basis_nodes,basis_edges,x,u_vec_r,u_vec_th,u_vec_z,n)
% ERRORS_EXACT_WEIGHTED_MODIFIED_NEDELEC_AND_P1 - Calculate the errors of our solution x
% compared to the exact solution u.
%
% Syntax:
%     [err,grad_err,max_err] = 
%         errors_exact_weighted_modified_nedelec_and_p1(p,t,t_ed,basis,x,u_vec_r,u_vec_z)
% Inputs:
%     p - a 2xNumNodes matrix representing nodal coordinates.
%     t - a 4xNumTriangles matrix representing the element connectivity in 
%         terms of node IDs. The end row of T represents the geometry face ID 
%         to which the element belongs
%     t_ed - a 3xNumTriangles matrix representing the which edges
%         correspond to which triangles. t_ed(i,T) represents the ith edge
%         in triangle T.
%     basis - a 3x3xNumTriangles matrix representing piece-wise basis 
%         functions for each node in each triangle. basis(i,:,T) represents 
%         the pieceiwise basis function for the ith node in triangle T. 
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
[~,nodes] = size(p);
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

    In = basis_nodes(:,1,T);
    Jn = basis_nodes(:,2,T);
    Kn = basis_nodes(:,3,T);

    ai = In(1); bi = In(2); ci = In(3);
    aj = Jn(1); bj = Jn(2); cj = Jn(3);
    ak = Kn(1); bk = Kn(2); ck = Kn(3);
     
    ni = t(1,T);
    nj = t(2,T);
    nk = t(3,T);
 
    Ie = basis_edges(:,1,T);
    Je = basis_edges(:,2,T);
    Ke = basis_edges(:,3,T);

    Ai = Ie(1); Bi = Ie(2); Ci = Ie(3);
    Aj = Je(1); Bj = Je(2); Cj = Je(3);
    Ak = Ke(1); Bk = Ke(2); Ck = Ke(3);

    ei = t_ed(1,T) + nodes;
    ej = t_ed(2,T) + nodes;
    ek = t_ed(3,T) + nodes;

    approx_r =@(r,z) x(ni).*((-1./n).*(ci + ai.*r + bi.*z)) ...
        + x(nj).*((-1./n).*(cj + aj.*r + bj.*z)) ...
        + x(nk).*((-1./n).*(ck + ak.*r + bk.*z)) ...
        + x(ei).*((1./n).*(Bi.*r - Ai.*r.*z)) ...
        + x(ej).*((1./n).*(Bj.*r - Aj.*r.*z)) ...
        + x(ek).*((1./n).*(Bk.*r - Ak.*r.*z));
    
    approx_th =@(r,z) x(ni).*(ci + ai.*r + bi.*z) ...
        + x(nj).*(cj + aj.*r + bj.*z) ...
        + x(nk).*(ck + ak.*r + bk.*z);
    % theta component of edge basis functions is 0
    
    approx_z =@(r,z) x(ei).*((1./n).*(Ci.*r + Ai.*r.^2)) ...
        + x(ej).*((1./n).*(Cj.*r + Aj.*r.^2)) ...
        + x(ek).*((1./n).*(Ck.*r + Ak.*r.^2));
    % z component of nodal basis functions is 0
    
    % find L2 Error
    integrand =@(r,z) ((u_vec_r(r,z) - approx_r(r,z)).^2 ...
        + (u_vec_th(r,z) - approx_th(r,z)).^2 ...
        + (u_vec_z(r,z) - approx_z(r,z)).^2).*r;
    
    integral = integral + Wx'*feval(integrand,X,Y)*Wy;
end

err = sqrt(integral);
% end