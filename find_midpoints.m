function [p2,t2] = find_midpoints(p,t)
% FIND_MIDPOINTS - Find the midpoints between adjacent nodes in the mesh
%
% Syntax:
%     [p2,t2] = find_midpoints(p,t)
% Inputs:
%     p - a 2xNumNodes matrix representing nodal coordinates.
%     t - a 4xNumTriangles matrix representing the element connectivity in 
%         terms of node IDs. The end row of T represents the geometry face 
%         ID to which the element belongs.
%
% Outputs:
%     p2 - a 2xNumNodes matrix representing midpoint nodal coordinates.
%     t2 - a 3xNumTriangles matrix representing the element connectivity in 
%         terms of node IDs. The three node IDs in a column are the three
%         midpoints of the node IDS in corresponding column in t.
%
% Author: Nicole Stock
% Date: Spring 2020

[~,triangles] = size(t);
[~,nodes] = size(p);
p2 = NaN(2,nodes);
t2 = NaN(3,triangles);
counter_p2 = 1;

% for each triangle
for T = 1:triangles
    % get coordinate points of the 3 nodes in triangle T
    xs = NaN(3,1);
    ys = NaN(3,1);
    for i = 1:3
        node = t(i, T);
        % get x,y coordinates
        xs(i) = p(1,node);
        ys(i) = p(2,node);
    end

    % find midpoints for each side
    mid_xs = NaN(3,1);
    mid_ys = NaN(3,1);
    % x values
    mid_xs(1) = (xs(1) + xs(2))/2;
    mid_xs(2) = (xs(2) + xs(3))/2;
    mid_xs(3) = (xs(3) + xs(1))/2;
    % y values
    mid_ys(1) = (ys(1) + ys(2))/2;
    mid_ys(2) = (ys(2) + ys(3))/2;
    mid_ys(3) = (ys(3) + ys(1))/2;

    % check for duplicates
    [isfound,index] = ismember([mid_xs,mid_ys], p2.', 'rows');
    % add midpoints to p2 and t2
    for i = 1:3
        if ~isfound(i)
            p2(1,counter_p2) = mid_xs(i);
            p2(2,counter_p2) = mid_ys(i);
            t2(i,T) = counter_p2;
            counter_p2 = counter_p2 + 1;
        else
            t2(i,T) = index(i);
        end
    end
end
    
end