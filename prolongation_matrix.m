function prolongation_matrix = prolongation_matrix(nodes_km1,mid_nodes_km1,p,p2,t_km1,t2_km1,basis_km1,n)
% PROLONGATION_MATRIX - Create prolongation matrix
%
% Syntax:
%     prolongation_matrix = prolongation_matrix()
%
% Inputs:
%     nodes_km1 - Number of vertex nodes in mesh level - 1.
%     mid_nodes_km1 - Number of midpoint nodes in mesh level - 1.
%     p - a 2xNumNodes matrix representing nodal coordinates.
%     t_km1 - a 4xNumTriangles matrix representing the element connectivity
%         in terms of node IDs. The end row of T represents the geometry
%         face ID to which the element belongs. For current mesh level - 1.
%     p2 - a 2xNumNodes matrix representing midpoint nodal coordinates.
%     t2_km1 - a 4xNumTriangles matrix representing the element 
%         connectivity in terms of node IDs. The three node IDs in a column
%         are the three midpoints of the node IDS in corresponding column
%         in t. For current mesh level - 1.
%     basis_km1 - a 3x3xNumTriangles matrix representing piece-wise basis 
%         functions for each node in each triangle. basis(i,:,k) represents 
%         the pieceiwise basis function for the ith node in triangle k. For
%         current mesh level - 1.
%     n - Hodge Laplacian on Axisymmetrix Domain and its Discretization
%     weight
%    
% Outputs:
%     prolongation_matrix - Prolongation matrix. Used to extend mesh level
%     - 1 vectors.
%
% Author: Nicole Stock
% Date: Spring 2020

[~,nodes] = size(p);
[~,mid_nodes] = size(p2);
total_nodes = nodes + mid_nodes;
total_nodes_km1 = nodes_km1 + mid_nodes_km1;

% P is a total_nodes x total_nodes_km1 matrix
%   ij entry of P : basis_km1_j evaluated at node i in p / p2

[~,triangles_km1] = size(t_km1);
i_vec = zeros(1,triangles_km1*6*total_nodes);
j_vec = zeros(1,triangles_km1*6*total_nodes);
s_vec = zeros(1,triangles_km1*6*total_nodes);
index = 1;

for T = 1:triangles_km1
    for j = 1:6
        J = basis_km1(:,j,T);
        
        a = J(1);
        b = J(2);
        c = J(3);
        d = J(4);
        e = J(5);
        f = J(6);
        
        if n == 0
            basis_fn_j =@(r,z) a.*r.^2 + b.*r.*z + c.*z.^2 + d.*r ...
                        + e.*z + f;
        else
            basis_fn_j =@(r,z) (r./n).*(a.*r.^2 + b.*r.*z + c.*z.^2 ...
                        + d.*r + e.*z + f);
        end
        
        if j >= 4
            global_j = t2_km1(j-3,T) + nodes_km1;
        else
            global_j = t_km1(j,T);
        end
        
        for i = 1:nodes
            ij_val = feval(basis_fn_j,p(1,i),p(2,i));
            
            i_vec(index) = i;
            j_vec(index) = global_j;
            s_vec(index) = ij_val;
            index = index + 1;
        end  
            
        for i = 1:mid_nodes 
            ij_val = feval(basis_fn_j,p2(1,i),p2(2,i));

            i_vec(index) = i + nodes;
            j_vec(index) = global_j;
            s_vec(index) = ij_val;
            index = index + 1;
        end
    end
end

prolongation_matrix = sparse(i_vec,j_vec,s_vec,total_nodes,total_nodes_km1);

% end