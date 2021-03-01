function mass_matrix = mass_matrix_weighted_HL_k_2_p1(p,t,ed,t_ed,basis_nodes,basis_edges,n)
% MASS_MATRIX_WEIGHTED_HL_K_2_P1 - Create mass matrix
%   Hodge Laplacian k = 2 case, lowest order
%   (zeta_i, zeta_j)_r where {zeta_k}k=1->(N+Ne) is the basis for Bh
%   (Bh is the weighted fourier modified Nedelec and P1 space)
%
% Syntax:
%     A = mass_matrix_weighted_HL_k_2_p1(p,t,ed,t_ed,basis_nodes,basis_edges,n)
%
% Inputs:
%     p - a 2xNumNodes matrix representing nodal coordinates.
%     t - a 4xNumTriangles matrix representing the element connectivity in
%         terms of node IDs. The end row of T represents the geometry face 
%         ID to which the element belongs.
%     basis_nodes - a matrix representing piece-wise basis functions for
%         each node in each triangle. basis(i,:,T) represents the
%         pieceiwise basis function for the ith node in triangle T.
%     basis_edges - a matrix representing piece-wise basis functions
%         for each edge in each triangle for the weighted fourier Nedelec 
%         and P1 space. basis(i,:,T) represents the pieceiwise basis 
%         function for the ith edge in triangle T.
%     n - Hodge Laplacian on Axisymmetrix Domain and its discretization
%     weight
%
% Outputs:
%     mass_matrix - mass matrix used to solve system of equations to
%         approximate solution
%
% Author: Nicole Stock
% Date: Fall 2020

[~,triangles] = size(t);
[~,nodes] = size(p);
[edges,~] = size(ed);
i_vec = zeros(1,triangles*36);
j_vec = zeros(1,triangles*36);
s_vec = zeros(1,triangles*36);
index = 1;

for T = 1:triangles
    
    % get coordinates of triangle T
    coordinates = zeros(3,2);
    for N = 1:3
        node = t(N,T);
        % get x,y coordinates of triangle
        coordinates(N,:) = p(:,node);
    end
        
    [R,Z,Wr,Wz] = triquad(7, coordinates);
    
    % integrate for each pair of vertices in the triangle
    for i = 1:6
        for j = i:6
            if i <= 3
                I = basis_nodes(:,i,T);
                ai = I(1);
                bi = I(2);
                ci = I(3);
                zeta_i_r =@(r,z) (-1./n).*(ci + ai.*r + bi.*z);
                zeta_i_th =@(r,z) ci + ai.*r + bi.*z;
                zeta_i_z =@(r,z) 0;
                global_i = t(i,T) + edges;
            else
                I = basis_edges(:,i-3,T);
                Ai = I(1);
                Bi = I(2);
                Ci = I(3);
                zeta_i_r =@(r,z) (1./n).*(Bi.*r - Ai.*r.*z);
                zeta_i_th =@(r,z) 0;
                zeta_i_z =@(r,z) (1./n).*(Ci.*r + Ai.*r.^2);
                global_i = t_ed(i-3,T);
            end
            if j <= 3
                J = basis_nodes(:,j,T);
                aj = J(1);
                bj = J(2);
                cj = J(3);
                zeta_j_r =@(r,z) (-1./n).*(cj + aj.*r + bj.*z);
                zeta_j_th =@(r,z) cj + aj.*r + bj.*z;
                zeta_j_z =@(r,z) 0;
                global_j = t(j,T) + edges;
            else
                J = basis_edges(:,j-3,T);
                Aj = J(1);
                Bj = J(2);
                Cj = J(3);
                zeta_j_r =@(r,z) (1./n).*(Bj.*r - Aj.*r.*z);
                zeta_j_th =@(r,z) 0;
                zeta_j_z =@(r,z) (1./n).*(Cj.*r + Aj.*r.^2);
                global_j = t_ed(j-3,T);
            end

            integrand =@(r,z) (zeta_i_r(r,z).*zeta_j_r(r,z) ...
                + zeta_i_th(r,z).*zeta_j_th(r,z) ...
                + zeta_i_z(r,z).*zeta_j_z(r,z)).*r;
            
            Q = Wr'*feval(integrand,R,Z)*Wz;               
          
            i_vec(index) = global_i;
            j_vec(index) = global_j;
            s_vec(index) = Q;
            index = index + 1;
            if global_j ~= global_i
                i_vec(index) = global_j;
                j_vec(index) = global_i;
                s_vec(index) = Q;
                index = index + 1;
            end
        end
    end
end

mass_matrix = sparse(i_vec,j_vec,s_vec,nodes+edges,nodes+edges);
% end
