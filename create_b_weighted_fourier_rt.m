function b = create_b_weighted_fourier_rt(p,t,ed,t_ed,basis_edges,basis_triangles,f_vec_r,f_vec_th,f_vec_z,n)
% CREATE_B_WEIGHTED_FOURIER_RT - Create vector b such that 
%     stiffness_matrix * solution = b.
%
% Syntax:
%     b = create_b_weighted_fourier_rt(p,t,ed,t_ed,basis,f_vec_r,f_vec_th,f_vec_z,n)
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
%     basis_edges - a matrix representing piece-wise basis functions for 
%         each edge in each triangle. basis(i,:,T) represents the 
%         pieceiwise basis function for the ith edge in triangle T.
%     basis_triangles - a vector representing piece-wise basis functions
%         for edge triangle. basis(1,T) represents the piecewise basis
%         function for the Tth triangle.
%     f_vec_r - given vector r component
%     f_vec_th - given vector theta component
%     f_vec_z - given vector z component
%     n - Hodge Laplacian on Axisymmetrix Domain and its discretization
%     weight
%
% Outputs:
%     b - vector such that stiffness_matrix * solution = b.
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
    
    % get coordinates of triangle T (Tth col of t)
    coordinates = zeros(3,2);
    for N = 1:3
        node = t(N,T);
        % get x,y coordinates of triangle
        coordinates(N,:) = p(:,node);
    end
    
    [X,Y,Wx,Wy] = triquad(7, coordinates);
    
    % for each edge + 1 in the triangle
    for i = 1:4
        if i <= 3
            I = basis_edges(:,i,T);
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
    
        integrand =@(r,z) (f_vec_r(r,z).*phi_i_r(r,z) ...
            + f_vec_th(r,z).*phi_i_th(r,z) ...
            + f_vec_z(r,z).*phi_i_z(r,z)).*r;
        
        b_i = Wx'*feval(integrand,X,Y)*Wy;
                
        i_vec(index) = global_i;
        s_vec(index) = b_i;
        index = index + 1;
    end
end

b = sparse(i_vec,j_vec,s_vec,edges+triangles,1);