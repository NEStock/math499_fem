function basis = basis_functions_rt(p,t,ed,t_ed)
% BASIS_FUNCTIONS_RT - Create a piecewise basis function for 
%   each node of a triangulation
%
% Syntax:
%     basis = basis_functions_rt(p,ed,t_ed)
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

[~,n] = size(t_ed);
basis = zeros(3,3,n);

% for each triangle (column in t_ed)
for T = 1:n
    
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
    
    % calculate unit outward normals per edge
    % using the fact that vertices are listed counter-clockwise
    norm_rs = zeros(3,1);
    norm_zs = zeros(3,1);
    
    sq_3 = sqrt((rs(2) - rs(1))^2 + (zs(2) - zs(1))^2);
    norm_rs(3) = (zs(2) - zs(1))/sq_3;
    norm_zs(3) = (-rs(2) + rs(1))/sq_3;
    sq_1 = sqrt((rs(3) - rs(2))^2 + (zs(3) - zs(2))^2);
    norm_rs(1) = (zs(3) - zs(2))/sq_1;
    norm_zs(1) = (-rs(3) + rs(2))/sq_1;
    sq_2 = sqrt((rs(1) - rs(3))^2 + (zs(1) - zs(3))^2);
    norm_rs(2) = (zs(1) - zs(3))/sq_2;
    norm_zs(2) = (-rs(1) + rs(3))/sq_2;
    
    sigma = zeros(3,1);
    for i = 1:3
        if norm_rs(i) == 0
            % same or opp. direction as [0;1]
            if norm_zs(i) == 1
                sigma(i) = 1;
            else
                sigma(i) = -1;
            end
        elseif norm_zs(i) == 0
            % same or opp. direction as [1;0]
            if norm_rs(i) == 1
                sigma(i) = 1;
            else
                sigma(i) = -1;
            end
        else
            % same or opp. direction as [sqrt(2)/2; sqrt(2)/2]
            if norm_rs(i) > 0
                sigma(i) = 1;
            else
                sigma(i) = -1;
            end
        end
    end
        
    % for each edge in triangle t (entry in in t_ed(:,T))
    for i = 1:3
        % get edge
        edge = t_ed(i,T);
        % get edge vertices
        v1 = ed(edge,1);
        v2 = ed(edge,2);
                
        % get r,z coordinates of each vertex
        r_v1 = p(1,v1);
        z_v1 = p(2,v1);
        r_v2 = p(1,v2);
        z_v2 = p(2,v2);

        % get length of edge
        Ei = sqrt((r_v2 - r_v1)^2 + (z_v2 - z_v1)^2);
        
        % get point i (vertex in T, not part of edge i)
        s_i = 0;
        for j = 1:3
            if v1 ~= t(j,T) && v2 ~= t(j,T)
                Pi = t(j,T);
                % edge corresponds to sigma j value
                s_i = j;
            end
        end
        
        % get r,z coordinates of Pi
        r_Pi = p(1,Pi);
        z_Pi = p(2,Pi);

        % phi = sigmai*(|Ei|/2|T|)*([r;z] - Pi)
        s = sigma(s_i)*(Ei/(2*Ta));
        
        basis(1,i,T) = s; % A
        basis(2,i,T) = -s*r_Pi; %B
        basis(3,i,T) = -s*z_Pi; %C
    end
end

end


