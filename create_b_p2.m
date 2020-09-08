function b = create_b_p2(p,t,p2,t2,basis,f)
% CREATE_B_P2 - Create vector b
%
% Syntax:
%     b = create_b_p2(p,e,t,basis,f)
%
% Inputs:
%     p - a 2xNumNodes matrix representing nodal coordinates.
%     t - a 4xNumTriangles matrix representing the element connectivity in 
%         terms of node IDs. The end row of T represents the geometry face ID 
%         to which the element belongs.
%     p2 - a 2xNumNodes matrix representing midpoint nodal coordinates.
%     t2 - a 4xNumTriangles matrix representing the element connectivity in 
%         terms of node IDs. The three node IDs in a column are the three
%         midpoints of the node IDS in corresponding column in t.
%     basis - a 3x3xNumTriangles matrix representing piece-wise basis 
%         functions for each node in each triangle. basis(i,:,k) represents 
%         the pieceiwise basis function for the ith node in triangle k. 
%
% Outputs:
%     b - vector such that stiffness_matrix * b = solution.
%
% Author: Nicole Stock
% Date: Spring 2020

[~,triangles] = size(t);
[~,nodes] = size(p);
i_vec = zeros(1,triangles*6);
j_vec = ones(1,triangles*6);
s_vec = zeros(1,triangles*6);
index = 1;

for T = 1:triangles
    
    % get coordinates of triangle T (Tth col of t)
    coordinates = zeros(3,2);
    for n = 1:3
        node = t(n,T);
        % get x,y coordinates
        coordinates(n,:) = p(:,node);
    end
    
    [X,Y,Wx,Wy] = triquad(7, coordinates);
    
    for i = 1:6
        
        I = basis(:,i,T);
    
        integrand =@(x,y) feval(f,x,y).*(I(1).*x.^2 + I(2).*x.*y ...
            + I(3).*y.^2 + I(4).*x + I(5).*y + I(6));
        b_i = Wx'*feval(integrand,X,Y)*Wy;
        
        if i >= 4
            global_i = t2(i-3,T) + nodes;
        else
            global_i = t(i,T);
        end
                
        i_vec(index) = global_i;
        s_vec(index) = b_i;
        index = index + 1;
    end
end

[~,n] = size(p);
[~,n2] = size(p2);
b = sparse(i_vec,j_vec,s_vec,n+n2,1);