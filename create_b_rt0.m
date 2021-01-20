function b = create_b_rt0(p,t,ed,t_ed,basis,f_vec_r,f_vec_z)
% CREATE_B_RT0 - Create vector b such that 
%     stiffness_matrix * solution = b.
%
% Syntax:
%     b = create_b_rt0(p,t,ed,t_ed,basis,f_vec_r,f_vec_z)
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
%     basis - an 3x2xNumTriangles matrix representing piece-wise basis 
%         functions for each node in each triangle. basis(i,:,T) represents 
%         the pieceiwise basis function for the ith node in triangle T. 
%     f_vec_r - given vector r component
%     f_vec_z - given vector z component
%
% Outputs:
%     b - vector such that stiffness_matrix * solution = b.
%
% Author: Nicole Stock
% Date: Fall 2020

[~,triangles] = size(t);
[edges,~] = size(ed);
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
    
    %f =@(r,z) [f_vec_r(r,z) ; f_vec_z(r,z)];
    
    % for each edge in the triangle
    for i = 1:3
        %integrand =@(r,z) dot(f(r,z), basis{i,T}(r,z));
        integrand =@(r,z) f_vec_r(r,z).*basis{i,1,T}(r,z) ...
            + f_vec_z(r,z).*basis{i,2,T}(r,z);
        
        b_i = Wx'*feval(integrand,X,Y)*Wy;
        %disp(b_i)
        
        global_i = t_ed(i,T);
        %disp(global_i)
        
        i_vec(index) = global_i;
        s_vec(index) = b_i;
        index = index + 1;
    end
end

N = edges;
b = sparse(i_vec,j_vec,s_vec,N,1);