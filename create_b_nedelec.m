function b = create_b_nedelec(p,t,ed,t_ed,basis,f_vec_r,f_vec_z)
% CREATE_B_NEDELEC - Create vector b
%
% Syntax:
%     b = create_b_nedelec(p,e,t,basis,f)
%
% Inputs:
%     p - a 2xNumNodes matrix representing nodal coordinates.
%     t - a 4xNumTriangles matrix representing the element connectivity in 
%         terms of node IDs. The end row of T represents the geometry face ID 
%         to which the element belongs
%     basis - a 3x3xNumTriangles matrix representing piece-wise basis 
%         functions for each node in each triangle. basis(i,:,k) represents 
%         the pieceiwise basis function for the ith node in triangle k. 
%     f - input function
%
% Outputs:
%     b - vector such that stiffness_matrix * b = solution.
%
% Author: Nicole Stock
% Date: Fall 2020

[~,triangles] = size(t);

i_vec = zeros(1,triangles*3);
j_vec = ones(1,triangles*3);
s_vec = zeros(1,triangles*3);
index = 1;

for T = 1:triangles
    
    % get coordinates of triangle T (Tth col of t)
    coordinates = zeros(3,2);
    for n = 1:3
        node = t(n,T);
        % get x,y coordinates of triangle
        coordinates(n,:) = p(:,node);
    end
    
    [X,Y,Wx,Wy] = triquad(7, coordinates);
    
    % for each edge in the triangle
    for i = 1:3
        I = basis(:,i,T);
            
        ai = I(1);
        bi = I(2);
        ci = I(3);
    
        phi_i_r = @(r,z) bi - ai.*z;
        phi_i_z = @(r,z) ci + ai.*r;
    
        integrand =@(r,z) (f_vec_r(r,z).*phi_i_r(r,z) ...
            + f_vec_z(r,z).*phi_i_z(r,z)).*r;
        
        b_i = Wx'*feval(integrand,X,Y)*Wy;
        
        global_i = t_ed(i,T);
        
        i_vec(index) = global_i;
        s_vec(index) = b_i;
        index = index + 1;
    end
end

[edges,~] = size(ed);
b = sparse(i_vec,j_vec,s_vec,edges,1);