function b = create_b_HL_k_0_p1(p,t,basis,f_fn,grad_f_r,grad_f_z,n)
% CREATE_B_HL_K_0_P1 - Create vector b such that
%     stiffness_matrix * solution = b.
%   Hodge Laplacian k = 0 case, P1
%
% Syntax:
%     b = create_b_HL_p1(p,t,p2,t2,basis,f_fn,grad_f_r,grad_f_z,n)
%
% Inputs:
%     p - a 2xNumNodes matrix representing nodal coordinates.
%     t - a 4xNumTriangles matrix representing the element connectivity in 
%         terms of node IDs. The end row of T represents the geometry face ID 
%         to which the element belongs.
%     basis - a 3x3xNumTriangles matrix representing piece-wise basis 
%         functions for each node in each triangle. basis(i,:,k) represents 
%         the pieceiwise basis function for the ith node in triangle k.
%     f - given function
%     grad_f_r - gradient(f) w.r.t r
%     grad_f_z - gradient(f) w.r.t z
%     n - Hodge Laplacian on Axisymmetrix Domain and its Discretization
%     weight
%
% Outputs:
%     b - vector such that stiffness_matrix * solution = b.
%
% Author: Nicole Stock
% Date: Fall 2020

[~,triangles] = size(t);
[~,nodes] = size(p);
i_vec = zeros(1,triangles*3);
j_vec = ones(1,triangles*3);
s_vec = zeros(1,triangles*3);
index = 1;

for T = 1:triangles
    
    % get coordinates of triangle T (Tth col of t)
    coordinates = zeros(3,2);
    for i = 1:3
        node = t(i,T);
        % get x,y coordinates
        coordinates(i,:) = p(:,node);
    end
    
    [X,Y,Wx,Wy] = triquad(7, coordinates);
    
    for i = 1:3
        
        I = basis(:,i,T);
        ai = I(1);
        bi = I(2);
        ci = I(3);
        
        % Calculate grad_rz^n* of basis function I
        grad_i_r =@(r,z) (1./n).*(2.*ai.*r + bi.*z + ci);
        grad_i_th =@(r,z) ai.*r + bi.*z + ci;
        grad_i_z =@(r,z) (1./n).*(bi.*r);
        
        integrand =@(r,z) (grad_f_r(r,z).*grad_i_r(r,z) ...
            + (n./r).*f_fn(r,z).*grad_i_th(r,z) ...
            + grad_f_z(r,z).*grad_i_z(r,z)).*r;
        
        % Integrate and solve
        b_i = Wx'*feval(integrand,X,Y)*Wy;
        
        global_i = t(i,T);
                
        i_vec(index) = global_i;
        s_vec(index) = b_i;
        index = index + 1;
    end
end

b = sparse(i_vec,j_vec,s_vec,nodes,1);