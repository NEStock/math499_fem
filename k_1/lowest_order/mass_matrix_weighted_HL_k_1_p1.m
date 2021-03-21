function mass_matrix = mass_matrix_weighted_HL_k_1_p1(p,t,basis_nodes,n)
% MASS_MATRIX_WEIGHTED_HL_K_1_P1 - Create mass matrix
%   Hodge Laplacian k = 1 case, P1
%   (phi_i, phi_j)_r where {phi_k}k=1->N is the basis for Ah
%   (Ah is the weighted P1 space)
%
% Syntax:
%     A = mass_matrix_weighted_HL_k_1_p1(p,t,ed,t_ed,basis_nodes,basis_edges,n)
%
% Inputs:
%     p - a 2xNumNodes matrix representing nodal coordinates.
%     t - a 4xNumTriangles matrix representing the element connectivity in
%         terms of node IDs. The end row of T represents the geometry face 
%         ID to which the element belongs.
%     basis_nodes - a matrix representing piece-wise basis functions for
%         each node in each triangle. basis(i,:,T) represents the
%         pieceiwise basis function for the ith node in triangle T.
%     n - Fourier mode
%
% Outputs:
%     mass_matrix - mass matrix used to solve system of equations to
%         approximate solution
%
% Author: Nicole Stock
% Date: Fall 2020

[~,triangles] = size(t);
[~,nodes] = size(p);
i_vec = zeros(1,triangles*9);
j_vec = zeros(1,triangles*9);
s_vec = zeros(1,triangles*9);
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
    for i = 1:3
        for j = i:3
            I = basis_nodes(:,i,T);
            ai = I(1);
            bi = I(2);
            ci = I(3);
            phi_i =@(r,z) (1./n).*(ci.*r + ai.*r.^2 + bi.*r.*z);
            global_i = t(i,T);
            
            J = basis_nodes(:,j,T);
            aj = J(1);
            bj = J(2);
            cj = J(3);
            phi_j =@(r,z) (1./n).*(cj.*r + aj.*r.^2 + bj.*r.*z);
            global_j = t(j,T);
             
            integrand =@(r,z) phi_i(r,z).*phi_j(r,z).*r;
            
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

mass_matrix = sparse(i_vec,j_vec,s_vec,nodes,nodes);

% end