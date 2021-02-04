function stiffness_matrix = stiffness_matrix_weighted_modified_nedelec_and_p1(p,t,ed,t_ed,basis_nodes,basis_edges,n)
% STIFFNESS_MATRIX_WEIGHTED_MODIFIED_NEDELEC_AND_P1 - Create stiffness matrix
%
% Syntax:
%     A = stiffness_matrix_weighted_nedelec(p,t,ed,t_ed,basis)
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
%     basis_nodes - a matrix representing piece-wise basis functions for
%         each node in each triangle. basis(i,:,T) represents the
%         pieceiwise basis function for the ith node in triangle T.
%     basis_edges - a matrix representing piece-wise basis functions for 
%         each edge in each triangle. basis(i,:,T) represents the 
%         pieceiwise basis function for the ith edge in triangle T.
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
[edges,~] = size(ed);
i_vec = zeros(1,triangles*36);
j_vec = zeros(1,triangles*36);
s_vec = zeros(1,triangles*36);
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
    
    % integrate for each pair of edges in the triangle
    for i = 1:6
        for j = i:6
            if i <= 3
                I = basis_nodes(:,i,T);
                ai = I(1);
                bi = I(2);
                ci = I(3);
                phi_i_r =@(r,z) (-1./n).*(ci + ai.*r + bi.*z);
                phi_i_th =@(r,z) ci + ai.*r + bi.*z;
                phi_i_z =@(r,z) 0;
                curl_i_r =@(r,z) -bi;
                curl_i_th =@(r,z) (-bi./n);
                curl_i_z =@(r,z) ai;
                global_i = t(i,T);
            else
                I = basis_edges(:,i-3,T);
                Ai = I(1);
                Bi = I(2);
                Ci = I(3);
                phi_i_r =@(r,z) (1./n).*(Bi.*r - Ai.*r.*z);
                phi_i_th =@(r,z) 0;
                phi_i_z =@(r,z) (1./n).*(Ci.*r + Ai.*r.^2);
                curl_i_r =@(r,z) -Ci - Ai.*r;
                curl_i_th =@(r,z) (-1./n).*(Ci + 3.*Ai.*r);
                curl_i_z =@(r,z) Bi - Ai.*z;
                global_i = t_ed(i-3,T) + nodes;
            end
            if j <= 3
                J = basis_nodes(:,j,T);
                aj = J(1);
                bj = J(2);
                cj = J(3);
                phi_j_r =@(r,z) (-1./n).*(cj + aj.*r + bj.*z);
                phi_j_th =@(r,z) cj + aj.*r + bj.*z;
                phi_j_z =@(r,z) 0;
                curl_j_r =@(r,z) -bj;
                curl_j_th =@(r,z) (-bj./n);
                curl_j_z =@(r,z) aj;
                global_j = t(j,T);
            else
                J = basis_edges(:,j-3,T);
                Aj = J(1);
                Bj = J(2);
                Cj = J(3);
                phi_j_r =@(r,z) (1./n).*(Bj.*r - Aj.*r.*z);
                phi_j_th =@(r,z) 0;
                phi_j_z =@(r,z) (1./n).*(Cj.*r + Aj.*r.^2);
                curl_j_r =@(r,z) -Cj - Aj.*r;
                curl_j_th =@(r,z) (-1./n).*(Cj + 3.*Aj.*r);
                curl_j_z =@(r,z) Bj - Aj.*z;
                global_j = t_ed(j-3,T) + nodes;
            end

            % We assume n > 0 for all of these problems!
            % integrate (curl_rz^n(I) (dot) curl_rz^n(J) + I(dot)J)*r

            dot_ij =@(r,z) phi_i_r(r,z).*phi_j_r(r,z) ...
                + phi_i_th(r,z).*phi_j_th(r,z) ...
                + phi_i_z(r,z).*phi_j_z(r,z);
            dot_curl_ij =@(r,z) curl_i_r(r,z).*curl_j_r(r,z) ...
                + curl_i_th(r,z).*curl_j_th(r,z) ...
                + curl_i_z(r,z).*curl_j_z(r,z);
            
            integrand =@(r,z) (dot_curl_ij(r,z) + dot_ij(r,z)).*r;
            
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

stiffness_matrix = sparse(i_vec,j_vec,s_vec,nodes+edges,nodes+edges);

% end