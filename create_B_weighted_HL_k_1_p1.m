function B = create_B_weighted_HL_k_1_p1(p,t,ed,t_ed,basis_nodes,basis_edges,n)
% CREATE_B_WEIGHTED_HL_K_1_P1 - Create B matrix
%   Hodge Laplacian k = 1 case, P1
%   (grad_rz^n(phi_i), zeta_j)_r where {phi_k}k=1->N is the basis for Ah
%     and {zeta_l}l=1->N+Ne is the basis for Bh
%   (Ah is the weighted P1 space)
%   (Bh is the weighted fourier modified Nedelec and P1 space)
%
% Syntax:
%     B = create_B_weighted_HL_k_1_p1(p,t,ed,t_ed,basis_nodes,basis_edges,n)
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
%     B - B matrix used to solve system of equations to approximate
%         solution
%
% Author: Nicole Stock
% Date: Fall 2020

[~,triangles] = size(t);
[~,nodes] = size(p);
[edges,~] = size(ed);
i_vec = zeros(1,triangles*18);
j_vec = zeros(1,triangles*18);
s_vec = zeros(1,triangles*18);
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
    for i = 1:3
        for j = 1:6
            I = basis_nodes(:,i,T);
            ai = I(1);
            bi = I(2);
            ci = I(3);
            grad_i_r =@(r,z) (1./n).*(2.*ai.*r + bi.*z + ci);
            grad_i_th =@(r,z) -(ai.*r + bi.*z + ci);
            grad_i_z =@(r,z) (1./n).*bi.*r;
            global_i = t(i,T);

            if j <= 3
                J = basis_nodes(:,j,T);
                aj = J(1);
                bj = J(2);
                cj = J(3);
                phi_j_r =@(r,z) (-1./n).*(cj + aj.*r + bj.*z);
                phi_j_th =@(r,z) cj + aj.*r + bj.*z;
                phi_j_z =@(r,z) 0;
                global_j = t(j,T);
            else
                J = basis_edges(:,j-3,T);
                Aj = J(1);
                Bj = J(2);
                Cj = J(3);
                phi_j_r =@(r,z) (1./n).*(Bj.*r - Aj.*r.*z);
                phi_j_th =@(r,z) 0;
                phi_j_z =@(r,z) (1./n).*(Cj.*r + Aj.*r.^2);
                global_j = t_ed(j-3,T) + nodes;
            end

            integrand =@(r,z) (grad_i_r(r,z).*phi_j_r(r,z) ...
                + grad_i_th(r,z).*phi_j_th(r,z) ...
                + grad_i_z(r,z).*phi_j_z(r,z)).*r;

            Q = Wr'*feval(integrand,R,Z)*Wz;               

            i_vec(index) = global_i;
            j_vec(index) = global_j;
            s_vec(index) = Q;
            index = index + 1;
        end
    end
end

B = sparse(i_vec,j_vec,s_vec,nodes,nodes+edges);

% end