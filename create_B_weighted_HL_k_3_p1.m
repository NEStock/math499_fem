function B = create_B_weighted_HL_k_3_p1(p,t,ed,t_ed,basis_edges,basis_triangles)
% CREATE_B_WEIGHTED_HL_K_3_P1 - Create mass matrix
%   Hodge Laplacian k = 3 case, P1
%   (chi_i, div_rz^n(psi_j))_r where {psi_i}i=1->Ne+Nt is the basis for Ch
%     and {chi_j}j=1->Nt is the basis for Dh
%   (Ch is the weighted fourier Raviart Thomas space)
%   (Dh is the piecewise constant space)
%
% Syntax:
%     B = create_B_weighted_HL_k_3_p1(p,t,ed,t_ed,basis)
%
% Inputs:
%     p - a 2xNumNodes matrix representing nodal coordinates.
%     t - a 4xNumTriangles matrix representing the element connectivity in
%         terms of node IDs. The end row of T represents the geometry face 
%         ID to which the element belongs.
%     ed - a 2xNumEdges matrix representing each edge as a row with
%         starting node in column 1 and the ending node in column 2.
%     t_ed - a 3xNumTriangles matrix representing the which edges
%         correspond to which triangles. t_ed(i,T) represents the ith edge
%         in triangle T.
%     basis_edges - a matrix representing piece-wise basis functions for 
%         each edge in each triangle. basis(i,:,T) represents the 
%         pieceiwise basis function for the ith edge in triangle T.
%     basis_triangles - a vector representing piece-wise basis functions
%         for edge triangle. basis(1,T) represents the piecewise basis
%         function for the Tth triangle.
%
% Outputs:
%     B - B matrix used to solve system of equations to approximate
%         solution
%
% Author: Nicole Stock
% Date: Fall 2020

[~,triangles] = size(t);
[edges,~] = size(ed);
i_vec = zeros(1,triangles*4);
j_vec = zeros(1,triangles*4);
s_vec = zeros(1,triangles*4);
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
    
    % integrate for each pair of edges + 1 in the triangle
    for i = 1:4
        if i <= 3
            I = basis_edges(:,i,T);
            Ai = I(1);
            div_i =@(r,z) 2.*Ai;
            global_i = t_ed(i,T) + triangles;
        else
            Di = basis_triangles(1,T);
            div_i =@(r,z) -Di;
            global_i = T; % + edges;
        end

        % We assume n > 0 for all of these problems!

        integrand =@(r,z) (div_i(r,z)).*r;

        Q = Wr'*feval(integrand,R,Z)*Wz;               

        i_vec(index) = global_i;
        j_vec(index) = T;
        s_vec(index) = Q;
        index = index + 1;
    end
end

B = sparse(i_vec,j_vec,s_vec,edges+triangles,triangles);

% end