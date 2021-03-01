function S = create_S_weighted_HL_k_1_p1(p,t,ed,t_ed,basis_nodes,basis_edges,n)
% CREATE_S_WEIGHTED_HL_K_1_P1 - Create S matrix
%   Hodge Laplacian k = 1 case, P1
%   (curl_rz^n(zeta_i), curl_rz^n(zeta_j))_r where {zeta_k}k=1->N+Ne is the
%     basis for Bh
%   (Bh is the weighted fourier modified Nedelec and P1 space)
%
% Syntax:
%     S = create_S_weighted_HL_k_1_p1(p,t,ed,t_ed,basis_nodes,basis_edges,n)
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
%     basis_edges - a matrix representing piece-wise basis functions for 
%         each edge in each triangle. basis(i,:,T) represents the 
%         pieceiwise basis function for the ith edge in triangle T.
%     n - Hodge Laplacian on Axisymmetrix Domain and its discretization
%     weight
%
% Outputs:
%     S - S matrix used to solve system of equations to approximate
%         solution
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
    
    % integrate for each pair of basis functions in the triangle
    for i = 1:6
        for j = i:6
            if i <=3
                I = basis_nodes(:,i,T);
                ai = I(1);
                bi = I(2);
                ci = I(3);
                curl_i_r =@(r,z) -bi;
                curl_i_th =@(r,z) (-bi./n);
                curl_i_z =@(r,z) ai;
                global_i = t(i,T);
            else
                I = basis_edges(:,i-3,T);
                Ai = I(1);
                Bi = I(2);
                Ci = I(3);
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
                curl_j_r =@(r,z) -bj;
                curl_j_th =@(r,z) (-bj./n);
                curl_j_z =@(r,z) aj;
                global_j = t(j,T);
            else
                J = basis_edges(:,j-3,T);
                Aj = J(1);
                Bj = J(2);
                Cj = J(3);
                curl_j_r =@(r,z) -Cj - Aj.*r;
                curl_j_th =@(r,z) (-1./n).*(Cj + 3.*Aj.*r);
                curl_j_z =@(r,z) Bj - Aj.*z;
                global_j = t_ed(j-3,T) + nodes;
            end

            integrand =@(r,z) (curl_i_r(r,z).*curl_j_r(r,z) ...
                + curl_i_th(r,z).*curl_j_th(r,z) ...
                + curl_i_z(r,z).*curl_j_z(r,z)).*r;

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

S = sparse(i_vec,j_vec,s_vec,nodes+edges,nodes+edges);

% end