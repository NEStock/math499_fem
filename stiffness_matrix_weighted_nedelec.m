function stiffness_matrix = stiffness_matrix_weighted_nedelec(p,t,ed,t_ed,basis)
% STIFFNESS_MATRIX_WEIGHTED_NEDELEC - Create stiffness matrix with weight n
%
% Syntax:
%     A = stiffness_matrix_weighted_nedelec(p,t,p2,t2,basis,k)
%
% Inputs:
%     p - a 2xNumNodes matrix representing nodal coordinates.
%     t - a 4xNumTriangles matrix representing the element connectivity in 
%         terms of node IDs. The end row of T represents the geometry face ID 
%         to which the element belongs
%     basis - a 3x3xNumTriangles matrix representing piece-wise basis 
%         functions for each node in each triangle. basis(i,:,k) represents 
%         the pieceiwise basis function for the ith node in triangle k.
%     n - given weight
%
% Outputs:
%     stiffness_matrix - stiffness matrix
%
% Author: Nicole Stock
% Date: Fall 2020

[~,triangles] = size(t);
[edges,~] = size(ed);
i_vec = zeros(1,triangles*9);
j_vec = zeros(1,triangles*9);
s_vec = zeros(1,triangles*9);
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
    for i = 1:3
        for j = i:3
            I = basis(:,i,T);
            J = basis(:,j,T);
            
            ai = I(1);
            bi = I(2);
            ci = I(3);

            aj = J(1);
            bj = J(2);
            cj = J(3);
            
            % We assume n > 0 for all of these problems!
            
            % integrate (curl(I) * curl(J) + I(dot)J)*r
            % curl_rz(v) = curl_rz[ v_r ; v_z ]
            %            = partial_deriv_z(v_r) - partial_deriv_r(v_z)
            %
            % phi_i = [ B - Ay ; C + Ax ]
            %       = [ B - Az ; C + Ar ]
            %
            % curl_rz(phi_i) = partial_deriv_z(B - Az) - partial_deriv_r(C + Ar)
            %                = -A - A
            %                = -2A

            curl_i = @(r,z) -2.*(ai);            
            curl_j = @(r,z) -2.*(aj);
            
            phi_i_r = @(r,z) bi - ai.*z;
            phi_i_z = @(r,z) ci + ai.*r;
            
            phi_j_r = @(r,z) bj - aj.*z;
            phi_j_z = @(r,z) cj + aj.*r;
            
            phi_i_dot_phi_j =@(r,z) (phi_i_r(r,z).*phi_j_r(r,z) ...
                + phi_i_z(r,z).*phi_j_z(r,z));
                        
            integrand =@(r,z) (curl_i(r,z).*curl_j(r,z) ...
                + phi_i_dot_phi_j(r,z)).*r;
            
            Q = Wr'*feval(integrand,R,Z)*Wz;               

            global_i = t_ed(i,T);
            global_j = t_ed(j,T);
                        
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

stiffness_matrix = sparse(i_vec,j_vec,s_vec,edges,edges);

% end