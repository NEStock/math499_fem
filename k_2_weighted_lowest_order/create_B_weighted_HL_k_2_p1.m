function B = create_B_weighted_HL_k_2_p1(p,t,ed,t_ed,basis_nodes,basis_NP1_edges,basis_RT_edges,basis_triangles,n)
% CREATE_B_WEIGHTED_HL_K_2_P1 - Create B matrix
%   Hodge Laplacian k = 2 case, P1
%   (curl_rz^n(zeta_i), psi_i)_r where {zeta_i}i=1->N+Ne is the basis for 
%     Bh and {psi_j}j=1->Ne+Nt is the basis for Ch
%   (Bh is the weighted fourier modified Nedelec and P1 space)
%   (Ch is the weighted fourier Raviart Thomas space)
%
% Syntax:
%     B = create_B_weighted_HL_k_2_p1(p,t,ed,t_ed,basis_nodes,basis_NP1_edges,basis_RT_edges,basis_triangles,n)
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
%     basis_nodes - a matrix representing piece-wise basis functions for
%         each node in each triangle. basis(i,:,T) represents the
%         pieceiwise basis function for the ith node in triangle T.
%     basis_NP1_edges - a matrix representing piece-wise basis functions
%         for each edge in each triangle for the weighted fourier Nedelec
%         and P1 space. basis(i,:,T) represents the pieceiwise basis
%         function for the ith edge in triangle T.
%     basis_RT_edges - a matrix representing piece-wise basis functions
%         for each edge in each triangle for the weighted fourier Raviart 
%         Thomas space. basis(i,:,T) represents the pieceiwise basis 
%         function for the ith edge in triangle T.
%     basis_triangles - a vector representing piece-wise basis functions
%         for edge triangle. basis(1,T) represents the piecewise basis
%         function for the Tth triangle.
%     n - Hodge Laplacian on Axisymmetrix Domain and its discretization
%     weight
%
% Outputs:
%     B - B matrix used to solve system of equations to approximate
%         solution
%
% Author: Nicole Stock
% Date: Fall 2020

[~,triangles] = size(t);
[~,nodes] = size(p);
[edges,~] = size(ed);
i_vec = zeros(1,triangles*24);
j_vec = zeros(1,triangles*24);
s_vec = zeros(1,triangles*24);
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
    
    % integrate for each pair of (edges+nodes) x (edges+1) in the triangle
    for i = 1:6
        for j = 1:4
            if i <= 3
                I = basis_NP1_edges(:,i,T);
                Ai = I(1);
                Bi = I(2);
                Ci = I(3);
                curl_i_r =@(r,z) -Ci - Ai.*r;
                curl_i_th =@(r,z) (-1./n).*(Ci + 3.*Ai.*r);
                curl_i_z =@(r,z) Bi - Ai.*z;
                global_i = t_ed(i,T);
            else
                I = basis_nodes(:,i-3,T);
                ai = I(1);
                bi = I(2);
                %ci = I(3);
                curl_i_r =@(r,z) -bi;
                curl_i_th =@(r,z) (-bi./n);
                curl_i_z =@(r,z) ai;
                global_i = t(i-3,T) + edges;                
            end
            if j <= 3
                J = basis_RT_edges(:,j,T);
                Aj = J(1);
                Bj = J(2);
                Cj = J(3);
                phi_j_r =@(r,z) Bj + Aj.*r;
                phi_j_th =@(r,z) (1./n).*(Bj + Aj.*r);
                phi_j_z =@(r,z) Cj + Aj.*z;
                global_j = t_ed(j,T);
            else
                Dj = basis_triangles(1,T);
                phi_j_r =@(r,z) 0;
                phi_j_th =@(r,z) (1./n).*Dj.*r;
                phi_j_z =@(r,z) 0;
                global_j = T + edges;
            end
            
            integrand =@(r,z) (curl_i_r(r,z).*phi_j_r(r,z) ...
                + curl_i_th(r,z).*phi_j_th(r,z) ...
                + curl_i_z(r,z).*phi_j_z(r,z)).*r;

            Q = Wr'*feval(integrand,R,Z)*Wz;               

            i_vec(index) = global_i;
            j_vec(index) = global_j;
            s_vec(index) = Q;
            index = index + 1;
        end
    end
end

B = sparse(i_vec,j_vec,s_vec,nodes+edges,edges+triangles);
% end