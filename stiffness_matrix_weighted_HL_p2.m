function stiffness_matrix = stiffness_matrix_weighted_p2(p,t,p2,t2,basis,n)
% STIFFNESS_MATRIX_WEIGHTED_P2 - Create stiffness matrix with weight n
%
% Syntax:
%     A = stiffness_matrix_weighted_p2(p,t,p2,t2,basis,k)
%
% Inputs:
%     p - a 2xNumNodes matrix representing nodal coordinates.
%     t - a 4xNumTriangles matrix representing the element connectivity in 
%         terms of node IDs. The end row of T represents the geometry face ID 
%         to which the element belongs
%     basis - a 3x3xNumTriangles matrix representing piece-wise basis 
%         functions for each node in each triangle. basis(i,:,k) represents 
%         the pieceiwise basis function for the ith node in triangle k.
%     n - given weight
%
% Outputs:
%     stiffness_matrix - stiffness matrix
%
% Author: Nicole Stock
% Date: Fall 2020

[~,triangles] = size(t);
[~,nodes] = size(p);
i_vec = zeros(1,triangles*36);
j_vec = zeros(1,triangles*36);
s_vec = zeros(1,triangles*36);
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
    for i = 1:6
        for j = i:6
            I = basis(:,i,T);
            J = basis(:,j,T);
            
            ai = I(1);
            bi = I(2);
            ci = I(3);
            di = I(4);
            ei = I(5);
            fi = I(6);
                        
            aj = J(1);
            bj = J(2);
            cj = J(3);
            dj = J(4);
            ej = J(5);
            fj = J(6);
            
            % integrate grad(I) * grad(J)
%%%% DOUBLE CHECK -- n > 0 for all of these problems!
%             if n == 0
%                 grad_i_r =@(r,z) 2.*I(1).*r + I(2).*z + I(4);
%                 grad_i_z =@(r,z) I(2).*r + 2.*I(3).*z + I(5);
%                 grad_j_r =@(r,z) 2.*J(1).*r + J(2).*z + J(4);
%                 grad_j_z =@(r,z) J(2).*r + 2.*J(3).*z + J(5);  
%                 
%                 grad_integrand =@(r,z) (grad_i_r(r,z).*grad_j_r(r,z) ...
%                 + grad_i_z(r,z).*grad_j_z(r,z)).*r;
%                 % integral grad(I) * grad(J)
%                 Q1 = Wr'*feval(grad_integrand,R,Z)*Wz;
%                 
%                 integrand =@(r,z) (I(1).*r.^2 + I(2).*r.*z + I(3).*z.^2 ...
%                     + I(4).*r + I(5).*z + I(6)).*(J(1).*r.^2 + J(2).*r.*z ...
%                     + J(3).*z.^2 + J(4).*r + J(5).*z + J(6)).*r;
%                 % integral I * J
%                 Q2 = Wr'*feval(integrand,R,Z)*Wz;
%             else

            % grad^k_rz(v) = [ partial_deriv_r(v)
            %                 (-k/r)*v
            %                 partial_deriv_z(v) ]
            % phi_k_i = (r/k)(ar^2 + brz + cz^2 + dr + ez + f)
            
            grad_i_r =@(r,z) (1./n).*(3.*ai.*r.^2 + 2.*bi.*r.*z + ci.*z.^2 ...
            + 2.*di.*r + ei.*z + fi);
            n_r_I =@(r,z) ai.*r.^2 + bi.*r.*z + ci.*z.^2 + di.*r ...
                + ei.*z + fi;
            grad_i_z =@(r,z) (1./n).*(bi.*r.^2 + 2.*ci.*r.*z + ei.*r);
            
            grad_j_r =@(r,z) (1./n).*(3.*aj.*r.^2 + 2.*bj.*r.*z + cj.*z.^2 ...
            + 2.*dj.*r + ej.*z + fj);
            n_r_J =@(r,z) aj.*r.^2 + bj.*r.*z + cj.*z.^2 + dj.*r ...
                + ej.*z + fj;
            grad_j_z =@(r,z) (1./n).*(bj.*r.^2 + 2.*cj.*r.*z + ej.*r);
            

            grad_integrand =@(r,z) (grad_i_r(r,z).*grad_j_r(r,z) ...
                + n_r_I(r,z).*n_r_J(r,z) ...
                + grad_i_z(r,z).*grad_j_z(r,z)).*r;
            Q = Wr'*feval(grad_integrand,R,Z)*Wz;               

            % if i or j == 4,5,6, they are midpoint nodes, so their global
            %   id is found in t2 added to the number of original nodes
            if i >= 4
                global_i = t2(i-3,T) + nodes;
            else
                global_i = t(i,T);
            end
            
            if j >=4
                global_j = t2(j-3,T) + nodes;
            else
                global_j = t(j,T);
            end
                        
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

[~,n2] = size(p2);
N = nodes + n2;
stiffness_matrix = sparse(i_vec,j_vec,s_vec,N,N);

% end