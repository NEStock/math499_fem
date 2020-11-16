function [basis_edges, basis_triangles] = basis_functions_weighted_fourier_rt(p,t,ed,t_ed)
% BASIS_FUNCTIONS_WEIGHTED_FOURIER_RT - Create a piecewise basis function for 
%   each node of a triangulation
%
% Syntax:
%     basis = basis_functions_weighted_fourier_rt(p,ed,t_ed)
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
%
% Outputs:
%     basis_edges - a matrix representing piece-wise basis functions for 
%         each edge in each triangle. basis(i,:,T) represents the  
%         pieceiwise basis function for the ith edge in triangle T.
%     basis_triangles - a matrix representing piece-wise basis functions 
%         for each triangle. basis(1,T) represents the pieceiwise basis 
%         function for triangle T.
%
% Author: Nicole Stock
% Date: Fall 2020

[~,m] = size(t_ed);
basis_triangles = zeros(1,m);

% basis_edges(:,i,T) = basis for the ith node of triangle T
basis_edges = basis_functions_rt(p,t,ed,t_ed);

% basis_triangles(1,T) = basis for the Tth triangle

% for each triangle (column in t_ed)
for T = 1:m
    
    % get coordinates of triangle T
    rs = zeros(3,1);
    zs = zeros(3,1);
    for N = 1:3
        node = t(N,T);
        % get x,y coordinates of triangle
        rs(N,1) = p(1,node);
        zs(N,1) = p(2,node);
    end
    
    % get area of triangle
    Ta = (1/2)*(rs(1)*(zs(2) - zs(3)) + rs(2)*(zs(3) - zs(1)) ...
         + rs(3)*(zs(1) - zs(2)));
    
    basis_triangles(1,T) = 1/Ta;
end

end


