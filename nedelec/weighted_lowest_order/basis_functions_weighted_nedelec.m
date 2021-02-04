function basis = basis_functions_weighted_nedelec(p,ed,t_ed)
% BASIS_FUNCTIONS_WEIGHTED_NEDELEC - Create a piecewise basis function for 
%   each node of a triangulation
%
% Syntax:
%     basis = basis_functions_weighted_nedelec(p,ed,t_ed)
% 
% Inputs:
%     p - a 2xNumNodes matrix representing nodal coordinates.
%     ed - a 2xNumEdges matrix representing each edge as a row with
%         starting node in column 1 and the ending node in column 2.
%     t_ed - a 3xNumTriangles matrix representing the which edges
%         correspond to which triangles. t_ed(i,T) represents the ith edge
%         in triangle T.
%
% Outputs:
%     basis - a matrix representing piece-wise basis functions for each
%         edge in each triangle. basis(i,:,T) represents the pieceiwise 
%         basis function for the ith edge in triangle T.
%
% Author: Nicole Stock
% Date: Fall 2020
    
[m,n] = size(t_ed);
basis = zeros(3,3,n);

% for each triangle (column in t_ed)
for T = 1:n
    Beta = zeros(3,3);
    % for each edge in triangle t (entry in in t_ed(:,T))
    for i = 1:m
        % get edge
        edge = t_ed(i,T);
        % get edge vertices
        v1 = ed(edge,1);
        v2 = ed(edge,2);
        % get x,y coordinates of each vertex
        x_v1 = p(1,v1);
        y_v1 = p(2,v1);
        x_v2 = p(1,v2);
        y_v2 = p(2,v2);

        Beta(i,1) = y_v2*x_v1 - y_v1*x_v2;
        Beta(i,2) = x_v2 - x_v1;
        Beta(i,3) = y_v2 - y_v1;
    end
    I = eye(3);
    basis(:,:,T) = Beta\I;
end

end


