function [basis, basis_div] = basis_functions_rt0(p,t,ed,t_ed)
% BASIS_FUNCTIONS_RT0 - Create a piecewise basis function for 
%   each node of a triangulation
%
% Syntax:
%     basis = basis_functions_rt0b(p,ed,t_ed)
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
%     basis - a 8xNumTriangles matrix representing piece-wise basis functions for each
%         edge in each triangle. basis(i,:,T) represents the pieceiwise 
%         basis function for the ith edge in triangle T.
%
% Author: Nicole Stock
% Date: Fall 2020

[~,triangles] = size(t_ed);
basis =  cell(3, 2, triangles);
basis_div =  cell(3, triangles);

% e-hat function
% e1r =@(r,z) sqrt(2).*r;
% e1z =@(r,z) sqrt(2).*z;
% 
% e2r =@(r,z) r - 1;
% e2z =@(r,z) z;
% 
% e3r =@(r,z) r;
% e3z =@(r,z) z - 1;

% for each triangle (column in t_ed)
for T = 1:triangles
    % get edges
    edge1 = t_ed(1,T);
    edge2 = t_ed(2,T);
    edge3 = t_ed(3,T);
    % get edge vertices
    vs = zeros(3,2);
    vs(1,1) = ed(edge1,1);
    vs(1,2) = ed(edge1,2);
    vs(2,1) = ed(edge2,1);
    vs(2,2) = ed(edge2,2);
    vs(3,1) = ed(edge3,1);
    vs(3,2) = ed(edge3,2);
    
    % get coordinates of triangle T
    rs = zeros(3,1);
    zs = zeros(3,1);
    for N = 1:3
        node = t(N,T);
        % get x,y coordinates of triangle
        rs(N,1) = p(1,node);
        zs(N,1) = p(2,node);
    end
    
    for M = 1:3
        node1 = t(1,T);
        node2 = t(2,T);
        
        if (node1 ~= vs(M,1) && node1 ~= vs(M,2))
            % node 1 is not in edge M
            vertex(1) = M;
        elseif (node2 ~= vs(M,1) && node2 ~= vs(M,2))
            % node 2 is not in edge M
            vertex(2) = M;
        else
            vertex(3) = M;
        end
    end
    
    % find the hypotenuse edge
    % set it equal to edge 1 (does not include vertex 1 as an endpoint)
    
    E1 = sqrt((rs(3) - rs(2))^2 + (zs(3) - zs(2))^2);
    E2 = sqrt((rs(1) - rs(3))^2 + (zs(1) - zs(3))^2);
    E3 = sqrt((rs(2) - rs(1))^2 + (zs(2) - zs(1))^2);
    
    norm_rs = zeros(3,1);
    norm_zs = zeros(3,1);
    
    norm_rs(3) = (zs(2) - zs(1))/E3;
    norm_zs(3) = (-rs(2) + rs(1))/E3;
    norm_rs(1) = (zs(3) - zs(2))/E1;
    norm_zs(1) = (-rs(3) + rs(2))/E1;
    norm_rs(2) = (zs(1) - zs(3))/E2;
    norm_zs(2) = (-rs(1) + rs(3))/E2;
    
    sigma = zeros(3,1);
    for i = 1:3
        if norm_rs(i) == 0
            % same or opp. direction as [0;1]
            if norm_zs(i) == 1
                sigma(i) = 1;
            else
                sigma(i) = -1;
            end
            % side opp vertex i is side 3
            % side 3
            ethree = vertex(i);
            three = i;
        elseif norm_zs(i) == 0
            % same or opp. direction as [1;0]
            if norm_rs(i) == 1
                sigma(i) = 1;
            else
                sigma(i) = -1;
            end
            % side 2
            etwo = vertex(i);
            two = i;
        else
            % same or opp. direction as [sqrt(2)/2; sqrt(2)/2]
            if norm_rs(i) > 0
                sigma(i) = 1;
            else
                sigma(i) = -1;
            end
            % side 1
            eone = vertex(i);
            one = i;
        end
    end
        
    Jta = rs(two) - rs(one);
    Jtb = rs(three) - rs(one);
    Jtc = zs(two) - zs(one);
    Jtd = zs(three) - zs(one);
    
    detJt = (rs(two) - rs(one)).*(zs(three) - zs(one)) - (rs(three) - rs(one)).*(zs(two) - zs(one));

    a1 = (zs(three) - zs(one))./detJt;
    a2 = (rs(one) - rs(three))./detJt;
    a3 = (- rs(one).*zs(three) + rs(three).*zs(one))./detJt;
    a4 = (zs(one) - zs(two))./detJt;
    a5 = (rs(two) - rs(one))./detJt;
    a6 = (- rs(two).*zs(one) + rs(one).*zs(two))./detJt;
    
    basis{eone,1,T} =@(r,z) sigma(one).*(1./detJt).*sqrt(2).*(Jta.*(a1.*r ...
        + a2.*z + a3) + Jtb.*(a4.*r + a5.*z + a6));
    basis{eone,2,T} =@(r,z) sigma(one).*(1./detJt).*sqrt(2).*(Jtc.*(a1.*r ...
        + a2.*z + a3) + Jtd.*(a4.*r + a5.*z + a6));

    basis{etwo,1,T} =@(r,z) sigma(two).*(1./detJt).*(Jta.*(a1.*r + a2.*z ...
        + a3 - 1) + Jtb.*(a4.*r + a5.*z + a6)); % r
    basis{etwo,2,T} =@(r,z) sigma(two).*(1./detJt).*(Jtc.*(a1.*r + a2.*z ...
        + a3 - 1) + Jtd.*(a4.*r + a5.*z + a6)); % z

    basis{ethree,1,T} =@(r,z) sigma(three).*(1./detJt).*(Jta.*(a1.*r + a2.*z ...
        + a3) + Jtb.*(a4.*r + a5.*z + a6 - 1));
    basis{ethree,2,T} =@(r,z) sigma(three).*(1./detJt).*(Jtc.*(a1.*r + a2.*z ...
        + a3) + Jtd.*(a4.*r + a5.*z + a6 - 1));

    basis_div{eone,T} = sigma(one).*(sqrt(2)./(detJt)).*(Jta.*a1 + Jtb.*a4) ...
        + sigma(one).*(sqrt(2)./(detJt)).*(Jtc.*a2 + Jtd.*a5);        
    basis_div{etwo,T} = sigma(two).*(1./detJt).*(Jta.*a1 + Jtb.*a4) ...
        +  sigma(two).*(1./detJt).*(Jtc.*a2 + Jtd.*a5);
    basis_div{ethree,T} = sigma(three).*(1./detJt).*(Jta.*a1 + Jtb.*a4) ...
        +  sigma(three).*(1./detJt).*(Jtc.*a2 + Jtd.*a5);
end

end
