function [err_u,err_s] = errors_exact_weighted_HL_k_2_p1(p,t,ed,t_ed,basis_nodes,basis_NP1_edges,basis_RT_edges,basis_triangles,u_h,u_vec_r,u_vec_th,u_vec_z,s_h,s_vec_r,s_vec_th,s_vec_z,n)
% ERRORS_EXACT_WEIGHTED_HL_K_2_P1 - Calculate the errors of our solution
% (u_h, s_h) compared to the exact solution (u, s).
%   Hodge Laplacian k = 2 case, P1
%
% Syntax:
%     [err,grad_err,max_err] = 
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
%     basis_NP1_edges - a matrix representing piece-wise basis functions
%         for each edge in each triangle for the weighted fourier Nedelec
%         and P1 space. basis(i,:,T) represents the pieceiwise basis
%         function for the ith edge in triangle T.
%     basis_RT_edges - a matrix representing piece-wise basis functions
%         for each edge in each triangle for the weighted fourier Raviart 
%         Thomas space. basis(i,:,T) represents the pieceiwise basis 
%         function for the ith edge in triangle T.
%     basis_triangles - a vector representing piece-wise basis functions
%         for edge triangle. basis(1,T) represents the piecewise basis
%         function for the Tth triangle.
%     u_h - approximated solution of u
%     u_vec_r - exact solution u vector r component
%     u_vec_r - exact solution u vector theta component
%     u_vec_z - exact solution u vector z component
%     s_h - approximated solution of s
%     s_vec_r - exact solution s vector r component
%     s_vec_r - exact solution s vector theta component
%     s_vec_z - exact solution s vector z component
%     n - Hodge Laplacian on Axisymmetrix Domain and its discretization
%     weight
%
% Outputs:
%    err_u - L2 error for u approximation
%    err_v - L2 error for v approximation
%
% Author: Nicole Stock
% Date: Fall 2020

[~,triangles] = size(t);
[edges,~] = size(ed);
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

    % s_h (weighted fourier Nedelec & P1 space)

    In = basis_nodes(:,1,T);
    Jn = basis_nodes(:,2,T);
    Kn = basis_nodes(:,3,T);
    
    ai = In(1); bi = In(2); ci = In(3);
    aj = Jn(1); bj = Jn(2); cj = Jn(3);
    ak = Kn(1); bk = Kn(2); ck = Kn(3);
    
    ni = t(1,T) + edges;
    nj = t(2,T) + edges;
    nk = t(3,T) + edges;
 
    IeNP1 = basis_NP1_edges(:,1,T);
    JeNP1 = basis_NP1_edges(:,2,T);
    KeNP1 = basis_NP1_edges(:,3,T);

    Ai = IeNP1(1); Bi = IeNP1(2); Ci = IeNP1(3);
    Aj = JeNP1(1); Bj = JeNP1(2); Cj = JeNP1(3);
    Ak = KeNP1(1); Bk = KeNP1(2); Ck = KeNP1(3);

    ei = t_ed(1,T);
    ej = t_ed(2,T);
    ek = t_ed(3,T);
    
    % compute approximations for s_h
    approx_r_s =@(r,z) s_h(ni).*((-1./n).*(ci + ai.*r + bi.*z)) ...
        + s_h(nj).*((-1./n).*(cj + aj.*r + bj.*z)) ...
        + s_h(nk).*((-1./n).*(ck + ak.*r + bk.*z)) ...
        + s_h(ei).*((1./n).*(Bi.*r - Ai.*r.*z)) ...
        + s_h(ej).*((1./n).*(Bj.*r - Aj.*r.*z)) ...
        + s_h(ek).*((1./n).*(Bk.*r - Ak.*r.*z));
    
    approx_th_s =@(r,z) s_h(ni).*(ci + ai.*r + bi.*z) ...
        + s_h(nj).*(cj + aj.*r + bj.*z) ...
        + s_h(nk).*(ck + ak.*r + bk.*z);
    % theta component of edge basis functions is 0
    
    approx_z_s =@(r,z) s_h(ei).*((1./n).*(Ci.*r + Ai.*r.^2)) ...
        + s_h(ej).*((1./n).*(Cj.*r + Aj.*r.^2)) ...
        + s_h(ek).*((1./n).*(Ck.*r + Ak.*r.^2));
    % z component of nodal basis functions is 0
        
    % find L2 Error for s
    integrand =@(r,z) ((s_vec_r(r,z) - approx_r_s(r,z)).^2 ...
        + (s_vec_th(r,z) - approx_th_s(r,z)).^2 ...
        + (s_vec_z(r,z) - approx_z_s(r,z)).^2).*r;
    
    integral_s = integral_s + Wx'*feval(integrand,X,Y)*Wy;    
    
    % u_h (weighted fourier RT space)
    
    IeRT = basis_RT_edges(:,1,T);
    JeRT = basis_RT_edges(:,2,T);
    KeRT = basis_RT_edges(:,3,T);

    Ai = IeRT(1); Bi = IeRT(2); Ci = IeRT(3);
    Aj = JeRT(1); Bj = JeRT(2); Cj = JeRT(3);
    Ak = KeRT(1); Bk = KeRT(2); Ck = KeRT(3);

    ei = t_ed(1,T);
    ej = t_ed(2,T);
    ek = t_ed(3,T);
        
    Di = basis_triangles(1,T);
     
    ti = T + edges;
    
    % compute approximation for u_h
    approx_r_u =@(r,z) u_h(ei).*(Bi + Ai.*r) + u_h(ej).*(Bj + Aj.*r) ...
        + u_h(ek).*(Bk + Ak.*r);
    % r component of triangle basis functions is 0
    
    approx_th_u =@(r,z) u_h(ei).*(1./n).*(Bi + Ai.*r) ...
        + u_h(ej).*(1./n).*(Bj + Aj.*r) + u_h(ek).*(1./n).*(Bk + Ak.*r) ...
        + u_h(ti).*(1./n).*Di.*r;
    
    approx_z_u =@(r,z) u_h(ei).*(Ci + Ai.*z) + u_h(ej).*(Cj + Aj.*z) ...
        + u_h(ek).*(Ck + Ak.*z);
    % z component of triangle basis functions is 0

    % find L2 Error for u
    integrand =@(r,z) ((u_vec_r(r,z) - approx_r_u(r,z)).^2 ...
        + (u_vec_th(r,z) - approx_th_u(r,z)).^2 ...
        + (u_vec_z(r,z) - approx_z_u(r,z)).^2).*r;
    
    integral_u = integral_u + Wx'*feval(integrand,X,Y)*Wy;
end
err_s = sqrt(integral_s);
err_u = sqrt(integral_u);
% end