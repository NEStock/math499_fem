function [t_ed] = find_edges(t,ed)
% FIND_EDGES - Find the edges in each triangulation in the mesh
%
% Syntax:
%     [ed, t_ed] = find_edges(t)
% Inputs:
%     t - a 4xNumTriangles matrix representing the element connectivity in 
%         terms of node IDs. The end row of T represents the geometry face 
%         ID to which the element belongs.
%     ed - a 2xNumNodes matrix representing each edge as a row with
%         starting node in column 1 and the ending node in column 2.
%
% Outputs:
%     t_ed - a 3xNumTriangles matrix representing the which edges
%         correspond to which triangles. The three edge IDs in a column are
%         the three edges in a particular triangle.
%
% Author: Nicole Stock
% Date: Fall 2020

[m,~] = size(ed);
[mt,nt] = size(t);
t_ed = NaN(mt-1,nt);

% for each edge (rep as a row in edges_)
for i = 1:m
    % for each triangle
    for T = 1:nt
        % if each vertex of edge is a member of triangle T
        if ismember(ed(i,1),t(1:3,T)) && ismember(ed(i,2),t(1:3,T))
            if isnan(t_ed(1,T))
                t_ed(1,T) = i;
            elseif isnan(t_ed(2,T))
                t_ed(2,T) = i;
            else
                t_ed(3,T) = i;
            end                
        end
    end
end