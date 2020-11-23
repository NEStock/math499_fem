function mass_matrix = mass_matrix_weighted_HL_k_3_p1(p,t,ed,t_ed,basis_edges,basis_triangles,n)
% MASS_MATRIX_WEIGHTED_HL_K_3_P1 - Create mass matrix
%
% Syntax:
%     A = mass_matrix_weighted_HL_k_3_p1(p,t,ed,t_ed,basis_edges,basis_triangles,n)
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
[edges,~] = size(ed);
i_vec = zeros(1,triangles*16);
j_vec = zeros(1,triangles*16);
s_vec = zeros(1,triangles*16);
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
        for j = i:4
            if i <= 3
                I = basis_edges(:,i,T);
                Ai = I(1);
                Bi = I(2);
                Ci = I(3);
                phi_i_r =@(r,z) Bi + Ai.*r;
                phi_i_th =@(r,z) (1./n).*(Bi + Ai.*r);
                phi_i_z =@(r,z) Ci + Ai.*z;
                global_i = t_ed(i,T) + triangles;
            else
                Di = basis_triangles(1,T);
                phi_i_r =@(r,z) 0;
                phi_i_th =@(r,z) (1./n).*Di.*r;
                phi_i_z =@(r,z) 0;
                global_i = T; % + edges;
            end
            if j <= 3
                J = basis_edges(:,j,T);
                Aj = J(1);
                Bj = J(2);
                Cj = J(3);
                phi_j_r =@(r,z) Bj + Aj.*r;
                phi_j_th =@(r,z) (1./n).*(Bj + Aj.*r);
                phi_j_z =@(r,z) Cj + Aj.*z;
                global_j = t_ed(j,T) + triangles;
            else
                Dj = basis_triangles(1,T);
                phi_j_r =@(r,z) 0;
                phi_j_th =@(r,z) (1./n).*Dj.*r;
                phi_j_z =@(r,z) 0;
                global_j = T; % + edges;
            end

            % We assume n > 0 for all of these problems!

            integrand =@(r,z) (phi_i_r(r,z).*phi_j_r(r,z) ...
                + phi_i_th(r,z).*phi_j_th(r,z) ...
                + phi_i_z(r,z).*phi_j_z(r,z)).*r;
            
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

mass_matrix = sparse(i_vec,j_vec,s_vec,edges+triangles,edges+triangles);

% end