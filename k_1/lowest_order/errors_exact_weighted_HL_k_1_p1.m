function [err_u,err_s] = errors_exact_weighted_HL_k_1_p1(p,t,t_ed,basis_nodes,basis_edges,u_h,u_vec_r,u_vec_th,u_vec_z,s_h,s,n)
% ERRORS_EXACT_WEIGHTED_HL_K_1_P1 - Calculate the errors of our solution x
% compared to the exact solution u.
%   Hodge Laplacian k = 1 case, P1
%
% Syntax:
%     [err_u,err_s] = 
%         errors_exact_weighted_HL_k_1_p1
% Inputs:
%     p - a 2xNumNodes matrix representing nodal coordinates.
%     t - a 4xNumTriangles matrix representing the element connectivity in 
%         terms of node IDs. The end row of T represents the geometry face ID 
%         to which the element belongs.
%     t_ed - a 3xNumTriangles matrix representing the which edges
%         correspond to which triangles. t_ed(i,T) represents the ith edge
%         in triangle T.
%     basis_nodes - a matrix representing piece-wise basis functions for
%         each node in each triangle. basis(i,:,T) represents the
%         pieceiwise basis function for the ith node in triangle T.
%     basis_edges - a matrix representing piece-wise basis functions for 
%         each edge in each triangle. basis(i,:,T) represents the 
%         pieceiwise basis function for the ith edge in triangle T.
%     u_h - approximated solution of u
%     u_vec_r - exact solution u vector r component
%     u_vec_r - exact solution u vector theta component
%     u_vec_z - exact solution u vector z component
%     s_h - approximated solution of s
%     s - exact solution s
%     n - Fourier mode
%
% Outputs:
%    err_u - L2 error for u approximation
%    err_v - L2 error for v approximation
%
% Author: Nicole Stock
% Date: Fall 2020

[~,triangles] = size(t);
[~,nodes] = size(p);
integral_u = 0;
integral_s = 0;

for T = 1:triangles
    
    % get coordinates of triangle T
    coordinates = zeros(3,2);
    for i = 1:3
        node = t(i,T);
        % get x,y coordinates
        coordinates(i,:) = p(:,node);
    end
        
    [X,Y,Wx,Wy] = triquad(7, coordinates);

    % find L2 Error for s
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

    % compute approximation for s_h
    approx_s =@(r,z) s_h(ni).*(1./n).*(ci.*r + ai.*r.^2 + bi.*r.*z) ...
        + s_h(nj).*(1./n).*(cj.*r + aj.*r.^2 + bj.*r.*z) ...
        + s_h(nk).*(1./n).*(ck.*r + ak.*r.^2 + bk.*r.*z);
    
    % compute approximations for u_h
    approx_r_u =@(r,z) u_h(ni).*((-1./n).*(ci + ai.*r + bi.*z)) ...
        + u_h(nj).*((-1./n).*(cj + aj.*r + bj.*z)) ...
        + u_h(nk).*((-1./n).*(ck + ak.*r + bk.*z)) ...
        + u_h(ei).*((1./n).*(Bi.*r - Ai.*r.*z)) ...
        + u_h(ej).*((1./n).*(Bj.*r - Aj.*r.*z)) ...
        + u_h(ek).*((1./n).*(Bk.*r - Ak.*r.*z));
    
    approx_th_u =@(r,z) u_h(ni).*(ci + ai.*r + bi.*z) ...
        + u_h(nj).*(cj + aj.*r + bj.*z) ...
        + u_h(nk).*(ck + ak.*r + bk.*z);
    % theta component of edge basis functions is 0
    
    approx_z_u =@(r,z) u_h(ei).*((1./n).*(Ci.*r + Ai.*r.^2)) ...
        + u_h(ej).*((1./n).*(Cj.*r + Aj.*r.^2)) ...
        + u_h(ek).*((1./n).*(Ck.*r + Ak.*r.^2));
    % z component of nodal basis functions is 0

    
    % find L2 Error for u
    integrand =@(r,z) ((u_vec_r(r,z) - approx_r_u(r,z)).^2 ...
        + (u_vec_th(r,z) - approx_th_u(r,z)).^2 ...
        + (u_vec_z(r,z) - approx_z_u(r,z)).^2).*r;
    
    integral_u = integral_u + Wx'*feval(integrand,X,Y)*Wy;
    
    % find L2 Error for p
    integrand =@(r,z) ((s(r,z) - approx_s(r,z)).^2).*r;
    
    integral_s = integral_s + Wx'*feval(integrand,X,Y)*Wy;    

    
end

err_u = sqrt(integral_u);
err_s = sqrt(integral_s);
% end