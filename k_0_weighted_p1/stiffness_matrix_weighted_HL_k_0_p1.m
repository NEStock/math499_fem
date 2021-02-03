function stiffness_matrix = stiffness_matrix_weighted_HL_k_0_p1(p,t,basis,n)
% STIFFNESS_MATRIX_WEIGHTED_HL_K_0 - Create stiffness matrix with weight n
%
% Syntax:
%     A = stiffness_matrix_weighted_HL_k_0(p,t,basis,k)
%
% Inputs:
%     p - a 2xNumNodes matrix representing nodal coordinates.
%     t - a 4xNumTriangles matrix representing the element connectivity in 
%         terms of node IDs. The end row of T represents the geometry face ID 
%         to which the element belongs.
%     basis - a 6x6xNumTriangles matrix representing piece-wise basis 
%         functions for each node in each triangle. basis(i,:,k) represents 
%         the pieceiwise basis function for the ith node in triangle k.
%     n - Hodge Laplacian on Axisymmetrix Domain and its discretization
%     weight
%
% Outputs:
%     stiffness_matrix - stiffness matrix
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
        % get x,y coordinates
        coordinates(N,:) = p(:,node);
    end
        
    [R,Z,Wr,Wz] = triquad(7, coordinates);
    
    % integrate for each pair of nodes in the triangle
    for i = 1:3
        for j = i:3
            I = basis(:,i,T);
            J = basis(:,j,T);
            
            ai = I(1);
            bi = I(2);
            ci = I(3);
                        
            aj = J(1);
            bj = J(2);
            cj = J(3);
            
            % integrate grad(I) * grad(J)
            % grad^k_rz(v) = [ partial_deriv_r(v)
            %                 (-k/r)*v
            %                 partial_deriv_z(v) ]
            % phi_k_i = (r/k)(ar + br + c)
            
            grad_i_r =@(r,z) (1./n).*(2.*ai.*r + bi.*z + ci);
            grad_i_th =@(r,z) ai.*r + bi.*z + ci;
            grad_i_z =@(r,z) (1./n).*(bi.*r);
            
            grad_j_r =@(r,z) (1./n).*(2.*aj.*r + bj.*z + cj);
            grad_j_th =@(r,z) aj.*r + bj.*z + cj;
            grad_j_z =@(r,z) (1./n).*(bj.*r);

            grad_integrand =@(r,z) (grad_i_r(r,z).*grad_j_r(r,z) ...
                + grad_i_th(r,z).*grad_j_th(r,z) ...
                + grad_i_z(r,z).*grad_j_z(r,z)).*r;
            Q = Wr'*feval(grad_integrand,R,Z)*Wz;               

            global_i = t(i,T);
            global_j = t(j,T);
                        
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

stiffness_matrix = sparse(i_vec,j_vec,s_vec,nodes,nodes);

% end