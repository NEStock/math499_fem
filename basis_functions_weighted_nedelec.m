function basis = basis_functions_weighted_nedelec(p,t,ed,t_ed)
% BASIS_FUNCTIONS_WEIGHTED_NEDELEC - Create a piecewise basis function for 
%   each node of a triangulation
%
% Syntax:
%     basis = basis_functions_weighted_nedelec(p,e,t)
% 
% Inputs:
%     p - a 2-by-NumNodes matrix representing nodal coordinates.
%     t - a matrix representing the element connectivity in terms of node
%         IDs. The end row of T represents the geometry face ID to which the
%         element belongs
%
% Outputs:
%     basis - a matrix representing piece-wise basis functions for each node
%         in each triangle. basis(i,:,T) represents the pieceiwise basis 
%         function for the ith node in triangle T.
%
% Author: Nicole Stock
% Date: Fall 2020
    
[m,n] = size(t_ed);
basis = zeros(3,3,n);

% for each triangle (column in t_ed/t)
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


