function F = create_F_weighted_HL_k_2_p1(p,t,ed,t_ed,basis_RT_edges,basis_triangles,f_r,f_th,f_z,n)
% CREATE_F_WEIGHTED_HL_K_2_P1 - Create F matrix
%   Hodge Laplacian k = 2 case, P1
%   (f, psi_i)_r where {psi_j}j=1->Ne+Nt is the basis for Ch
%   (Ch is the weighted fourier Raviart Thomas space)
%
% Syntax:
%     F = create_F_weighted_HL_k_2_p1(p,t,ed,t_ed,basis_nodes,basis_edges,n)
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
%     basis_RT_edges - a matrix representing piece-wise basis functions
%         for each edge in each triangle for the weighted fourier Raviart 
%         Thomas space. basis(i,:,T) represents the pieceiwise basis 
%         function for the ith edge in triangle T.
%     basis_triangles - a vector representing piece-wise basis functions
%         for edge triangle. basis(1,T) represents the piecewise basis
%         function for the Tth triangle.
%     f_vec_r - given function r component
%     f_vec_th - given function theta component
%     f_vec_z - given function z component
%
% Outputs:
%     F - F matrix used to solve system of equations to approximate
%         solution
%
% Author: Nicole Stock
% Date: Fall 2020

[~,triangles] = size(t);
[edges,~] = size(ed);
i_vec = zeros(1,triangles*4);
j_vec = ones(1,triangles*4);
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
            I = basis_RT_edges(:,i,T);
            Ai = I(1);
            Bi = I(2);
            Ci = I(3);
            phi_i_r =@(r,z) Bi + Ai.*r;
            phi_i_th =@(r,z) (1./n).*(Bi + Ai.*r);
            phi_i_z =@(r,z) Ci + Ai.*z;
            global_i = t_ed(i,T);
        else
            Di = basis_triangles(1,T);
            phi_i_r =@(r,z) 0;
            phi_i_th =@(r,z) (1./n).*Di.*r;
            phi_i_z =@(r,z) 0;
            global_i = T + edges;
        end
            
        integrand =@(r,z) (phi_i_r(r,z).*f_r(r,z) ...
            + phi_i_th(r,z).*f_th(r,z) ...
            + phi_i_z(r,z).*f_z(r,z)).*r;

        Q = Wr'*feval(integrand,R,Z)*Wz;               

        i_vec(index) = global_i;
        s_vec(index) = Q;
        index = index + 1;
    end
end

F = sparse(i_vec,j_vec,s_vec,edges+triangles,1);

% end
