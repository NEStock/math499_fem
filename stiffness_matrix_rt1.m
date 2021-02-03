function stiffness_matrix = stiffness_matrix_rt1(p,t,ed,t_ed,basis,basis_div)
% STIFFNESS_MATRIX_RT1 - Create stiffness matrix
%
% Syntax:
%     A = stiffness_matrix_rt1(p,t,ed,t_ed,basis)
%
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
%     basis - a 8x2xNumTriangles matrix representing basis functions for
%         each node in each triangle. basis(i,:,T) represents the basis 
%         function for the ith node in triangle T.
%     basis_div - an 8xNumTriangles matrix representing divergence of the
%         basis functions for each node in each triangle. basis(i,:,T)  
%         represents the basis function for the ith node in triangle T.
%
% Outputs:
%     stiffness_matrix - stiffness matrix
%
% Author: Nicole Stock
% Date: Fall 2020

[~,triangles] = size(t);
[edges,~] = size(ed);
i_vec = zeros(1,triangles*64);
j_vec = zeros(1,triangles*64);
s_vec = zeros(1,triangles*64);
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
    
    % integrate for each pair of basis functions in the triangle
    for i = 1:8
        for j = i:8
            phi_dot =@(r,z) basis{i,1,T}(r,z).*basis{j,1,T}(r,z) ...
                + basis{i,2,T}(r,z).*basis{j,2,T}(r,z);
            
            div_dot =@(r,z) basis_div{i,T}(r,z).*basis_div{j,T}(r,z);

            integrand =@(r,z) (phi_dot(r,z) + div_dot(r,z));
            
            Q = Wr'*feval(integrand,R,Z)*Wz;

            if (i <= 3)
                % edge basis functions set 1
                global_i = t_ed(i,T);
            elseif (i <= 6)
                % edge basis functions set 2
                global_i = t_ed(i-3,T) + edges;
            elseif (i == 7)
                % non-normal basis functions set 1 (4s)
                global_i = edges*2 + T;
            else
                % non-normal basis functions set 2 (5s)
                global_i = edges*2 + triangles + T;
            end

            if (j <= 3)
                % edge basis functions set 1
                global_j = t_ed(j,T);
            elseif (j <= 6)
                % edge basis functions set 2
                global_j = t_ed(j-3,T) + edges;
            elseif (j == 7)
                % non-normal basis functions set 1 (4s)
                global_j = edges*2 + T;
            else
                % non-normal basis functions set 2 (5s)
                global_j = edges*2 + triangles + T;
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

N = edges*2 + triangles*2;
stiffness_matrix = sparse(i_vec,j_vec,s_vec,N,N);

% end