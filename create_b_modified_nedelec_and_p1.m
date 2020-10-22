function b = create_b_modified_nedelec_and_p1(p,t,ed,t_ed,basis_nodes,basis_edges,n,f_vec_r,f_vec_th,f_vec_z)
% CREATE_B_MODIFIED_NEDELEC_AND_P1 - Create vector b such that 
%     stiffness_matrix * solution = b.
%
% Syntax:
%     b = create_b_modified_nedelec_and_p1(p,t,ed,t_ed,basis,f_vec_r,f_vec_th,f_vec_z)
%
% Inputs:
%     p - a 2xNumNodes matrix representing nodal coordinates.
%     t - a 4xNumTriangles matrix representing the element connectivity in 
%         terms of node IDs. The end row of T represents the geometry face ID 
%         to which the element belongs
%     ed - a 2xNumEdges matrix representing each edge as a row with
%         starting node in column 1 and the ending node in column 2.
%     t_ed - a 3xNumTriangles matrix representing the which edges
%         correspond to which triangles. t_ed(i,T) represents the ith edge
%         in triangle T.
%     basis_nodes - a matrix representing piece-wise basis functions for
%         each node in each triangle. basis(i,:,T) represents the
%         pieceiwise basis function for the ith node in triangle T.
%     basis_edges - a matrix representing piece-wise basis functions for 
%         each edge in each triangle. basis(i,:,T) represents the 
%         pieceiwise basis function for the ith edge in triangle T.
%     n - Hodge Laplacian on Axisymmetrix Domain and its discretization
%     weight
%     f_vec_r - given vector r component
%     f_vec_th - given vector theta component
%     f_vec_z - given vector z component
%
% Outputs:
%     b - vector such that stiffness_matrix * solution = b.
%
% Author: Nicole Stock
% Date: Fall 2020

[~,triangles] = size(t);
[~,nodes] = size(p);
[edges,~] = size(ed);
i_vec = zeros(1,triangles*6);
j_vec = ones(1,triangles*6);
s_vec = zeros(1,triangles*6);
index = 1;

for T = 1:triangles
    
    % get coordinates of triangle T (Tth col of t)
    coordinates = zeros(3,2);
    for N = 1:3
        node = t(N,T);
        % get x,y coordinates of triangle
        coordinates(N,:) = p(:,node);
    end
    
    [X,Y,Wx,Wy] = triquad(7, coordinates);
    
    % for each edge in the triangle
    for i = 1:6
        if i <= 3
            I = basis_nodes(:,i,T);
            ai = I(1);
            bi = I(2);
            ci = I(3);
            phi_i_r =@(r,z) (-1./n).*(ci + ai.*r + bi.*z);
            phi_i_th =@(r,z) ci + ai.*r + bi.*z;
            phi_i_z =@(r,z) 0;
            global_i = t(i,T);
        else
            I = basis_edges(:,i-3,T);
            Ai = I(1);
            Bi = I(2);
            Ci = I(3);
            phi_i_r =@(r,z) (1./n).*(Bi.*r - Ai.*r.*z);
            phi_i_th =@(r,z) 0;
            phi_i_z =@(r,z) (1./n).*(Ci.*r + Ai.*r.^2);
            global_i = t_ed(i-3,T) + nodes;
        end
    
        integrand =@(r,z) (f_vec_r(r,z).*phi_i_r(r,z) ...
            + f_vec_th(r,z).*phi_i_th(r,z) ...
            + f_vec_z(r,z).*phi_i_z(r,z)).*r;
        
        b_i = Wx'*feval(integrand,X,Y)*Wy;
                
        i_vec(index) = global_i;
        s_vec(index) = b_i;
        index = index + 1;
    end
end

b = sparse(i_vec,j_vec,s_vec,nodes+edges,1);