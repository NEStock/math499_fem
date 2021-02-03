function [basis, basis_div] = basis_functions_rt1(p,t,ed,t_ed)
% BASIS_FUNCTIONS_RT1 - Create a piecewise basis function for 
%   each node of a triangulation
%
% Syntax:
%     basis = basis_functions_rt1(p,ed,t_ed)
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
%     basis - a 8x2xNumTriangles matrix representing basis functions for
%         each node in each triangle. basis(i,:,T) represents the basis 
%         function for the ith node in triangle T.
%     basis_div - an 8xNumTriangles matrix representing divergence of the
%         basis functions for each node in each triangle. basis(i,:,T)  
%         represents the basis function for the ith node in triangle T.
%
% Author: Nicole Stock
% Date: Fall 2020

[~,triangles] = size(t_ed);
basis = cell(8, 2, triangles);
basis_div = cell(8,triangles);

% Guassian quadature points
g1 = 1./2 - sqrt(3)./6;
g2 = 1./2 + sqrt(3)./6;

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
            esix = vertex(i) + 3;
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
            efive = vertex(i) + 3;
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
            efour = vertex(i) + 3;
            one = i;
        end
    end

    Jta = rs(two) - rs(one);
    Jtb = rs(three) - rs(one);
    Jtc = zs(two) - zs(one);
    Jtd = zs(three) - zs(one);
    
    detJt = (rs(two) - rs(one)).*(zs(three) - zs(one)) - (rs(three) - rs(one)).*(zs(two) - zs(one));
    abs_detJt = abs(detJt);

    a1 = (zs(three) - zs(one))./detJt;
    a2 = (rs(one) - rs(three))./detJt;
    a3 = (- rs(one).*zs(three) + rs(three).*zs(one))./detJt;
    a4 = (zs(one) - zs(two))./detJt;
    a5 = (rs(two) - rs(one))./detJt;
    a6 = (- rs(two).*zs(one) + rs(one).*zs(two))./detJt;
    
    b1 = a1./(g1 - g2);
    b2 = a2./(g1 - g2);
    b3 = (a3 - g2)./(g1 - g2);
    b4 = a4./(g1 - g2);
    b5 = a5./(g1 - g2);
    b6 = (a6 - g2)./(g1 - g2);
    
    c1 = a1./(g2 - g1);
    c2 = a2./(g2 - g1);
    c3 = (a3 - g1)./(g2 - g1);
    c4 = a4./(g2 - g1);
    c5 = a5./(g2 - g1);
    c6 = (a6 - g1)./(g2 - g1);
    
    % perform piola transformation on each basis function for the edge
    
    % edge basis functions set 1
    basis{eone,1,T} =@(r,z) sigma(one).*sqrt(2).*(1./abs_detJt).*((a4.*r ...
        + a5.*z + a6 - g2)./(g1 - g2)).*(Jta.*(a1.*r + a2.*z + a3) ...
        + Jtb.*(a4.*r + a5.*z + a6));
    
    basis{eone,2,T} =@(r,z) sigma(one).*sqrt(2).*(1./abs_detJt).*((a4.*r ...
        + a5.*z + a6 - g2)./(g1 - g2)).*(Jtc.*(a1.*r + a2.*z + a3) ...
        + Jtd.*(a4.*r + a5.*z + a6));
       
    basis{etwo,1,T} =@(r,z) sigma(two).*(1./abs_detJt).*(Jta.*((a4.*r ...
        + a5.*z + a6 - g1)./(g2 - g1)).*(a1.*r + a2.*z + a3 - 1) ...
        + Jtb.*((a4.*r + a5.*z + a6 - g1)./(g2 - g1)).*(a4.*r + a5.*z ...
        + a6));
    basis{etwo,2,T} =@(r,z) sigma(two).*(1./abs_detJt).*(Jtc.*((a4.*r ...
        + a5.*z + a6 - g1)./(g2 - g1)).*(a1.*r + a2.*z + a3 - 1) ...
        + Jtd.*((a4.*r + a5.*z + a6 - g1)./(g2 - g1)).*(a4.*r + a5.*z ...
        + a6));   

    basis{ethree,1,T} =@(r,z) sigma(three).*(1./abs_detJt).*(Jta.*((a1.*r ...
        + a2.*z + a3 - g2)./(g1 - g2)).*(a1.*r + a2.*z + a3) ...
        + (Jtb.*((a1.*r + a2.*z + a3 - g2)./(g1 - g2)).*(a4.*r + a5.*z ...
        + a6 - 1)));
    basis{ethree,2,T} =@(r,z) sigma(three).*(1./abs_detJt).*(Jtc.*((a1.*r ...
        + a2.*z + a3 - g2)./(g1 - g2)).*(a1.*r + a2.*z + a3) ...
        + (Jtd.*((a1.*r + a2.*z + a3 - g2)./(g1 - g2)).*(a4.*r + a5.*z ...
        + a6 - 1)));
   
    % edge basis functions set 2
    
    basis{efour,1,T} =@(r,z) sigma(one).*sqrt(2).*(1./abs_detJt).*((a4.*r ...
        + a5.*z + a6 - g1)./(g2 - g1)).*(Jta.*(a1.*r + a2.*z + a3) ...
        + Jtb.*(a4.*r + a5.*z + a6));
            
    basis{efour,2,T} =@(r,z) sigma(one).*sqrt(2).*(1./abs_detJt).*((a4.*r ...
        + a5.*z + a6 - g1)./(g2 - g1)).*(Jtc.*(a1.*r + a2.*z + a3) ...
        + Jtd.*(a4.*r + a5.*z + a6));
     
    basis{efive,1,T} =@(r,z) sigma(two).*(1./abs_detJt).*(Jta.*((a4.*r + a5.*z ...
        + a6 - g2)./(g1 - g2)).*(a1.*r + a2.*z + a3 - 1) + Jtb.*((a4.*r ...
        + a5.*z + a6 - g2)./(g1 - g2)).*(a4.*r + a5.*z + a6));

    basis{efive,2,T} =@(r,z) sigma(two).*(1./abs_detJt).*(Jtc.*((a4.*r + a5.*z ...
        + a6 - g2)./(g1 - g2)).*(a1.*r + a2.*z + a3 - 1) + Jtd.*((a4.*r ...
        + a5.*z + a6 - g2)./(g1 - g2)).*(a4.*r + a5.*z + a6));
    
    basis{esix,1,T} =@(r,z) sigma(three).*(1./abs_detJt).*(Jta.*((a1.*r ...
        + a2.*z + a3 - g1)./(g2 - g1)).*(a1.*r + a2.*z + a3) ...
        + (Jtb.*((a1.*r + a2.*z + a3 - g1)./(g2 - g1)).*(a4.*r + a5.*z ...
        + a6 - 1)));
    basis{esix,2,T} =@(r,z) sigma(three).*(1./abs_detJt).*(Jtc.*((a1.*r ...
        + a2.*z + a3 - g1)./(g2 - g1)).*(a1.*r + a2.*z + a3) ...
        + (Jtd.*((a1.*r + a2.*z + a3 - g1)./(g2 - g1)).*(a4.*r + a5.*z ...
        + a6 - 1)));
    
    % non-normal basis functions
    basis{7,1,T} =@(r,z) (1./abs_detJt).*(a4.*r + a5.*z + a6).*(Jta.*(a1.*r ...
        + a2.*z + a3)+ Jtb.*(a4.*r + a5.*z + a6 - 1));
    
    basis{7,2,T} =@(r,z) (1./abs_detJt).*(a4.*r + a5.*z + a6).*(Jtc.*(a1.*r ...
        + a2.*z + a3)+ Jtd.*(a4.*r + a5.*z + a6 - 1));
    
    basis{8,1,T} =@(r,z) (1./abs_detJt).*(Jta.*(a1.*r + a2.*z + a3).*(a1.*r ...
        + a2.*z + a3 - 1) + Jtb.*(a1.*r + a2.*z + a3).*(a4.*r + a5.*z ...
        + a6));
    basis{8,2,T} =@(r,z) (1./abs_detJt).*(Jtc.*(a1.*r + a2.*z + a3).*(a1.*r ...
        + a2.*z + a3 - 1) + Jtd.*(a1.*r + a2.*z + a3).*(a4.*r + a5.*z ...
        + a6));
    
    % basis divergence functions
    
    basis_div{eone,T} =@(r,z) sigma(one).*sqrt(2).*(1./abs_detJt).*(Jta.*(2.*b4.*a1.*r ...
        + (b4.*a2 + b5.*a1).*z + b4.*a3 + b6.*a1) + Jtb.*(2.*b4.*a4.*r ...
        + (b4.*a5 + b5.*a4).*z + b4.*a6 + b6.*a4)) ...
        + sigma(one).*sqrt(2).*(1./abs_detJt).*(Jtc.*(b4.*a2.*r + b5.*a1.*r ...
        + 2.*b5.*a2.*z + b5.*a3 + b6.*a2) + Jtd.*(b4.*a5.*r + b5.*a4.*r ...
        + 2.*b5.*a5.*z + b5.*a6 + b6.*a5));
    
    basis_div{etwo,T} =@(r,z) sigma(two).*(1./abs_detJt).*(Jta.*(2.*c4.*a1.*r ...
        + c4.*a2.*z + c4.*(a3 - 1) + c5.*a1.*z + c6.*a1) ...
        + Jtb.*(2.*c4.*a4.*r + c4.*a5.*z + c4.*a6 + c5.*a4.*z + c6.*a4)) ...
        + sigma(two).*(1./abs_detJt).*(Jtc.*(c4.*a2.*r + c5.*a1.*r ...
        + 2.*c5.*a2.*z + c5.*(a3 - 1) + c6.*a2) + Jtd.*(c4.*a5.*r ...
        + c5.*a4.*r + 2.*c5.*a5.*z + c5.*a6 + c6.*a5));
    
    basis_div{ethree,T} =@(r,z) sigma(three).*(1./abs_detJt).*(Jta.*(2.*b1.*a1.*r ...
        + b1.*a2.*z + b1.*a3 + b2.*a1.*z + b3.*a1) + Jtb.*(2.*b1.*a4.*r ...
        + b1.*a5.*z + b1.*(a6 - 1) + b2.*a4.*z + b3.*a4)) ...
        + sigma(three).*(1./abs_detJt).*(Jtc.*(b1.*a2.*r + b2.*a1.*r ...
        + 2.*b2.*a2.*z + b2.*a3 + b3.*a2) + Jtd.*(b1.*a5.*r + b2.*a4.*r ...
        + 2.*b2.*a5.*z + b2.*(a6 - 1) + b3.*a5));
    
    basis_div{efour,T} =@(r,z) sigma(one).*sqrt(2).*(1./abs_detJt).*(Jta.*( 2.*c4.*a1.*r ...
        + c4.*a2.*z + c4.*a3 + c5.*a1.*z + c6.*a1) ...
        + Jtb.*(2.*c4.*a4.*r + c4.*a5.*z + c4.*a6 + c5.*a4.*z ...
        + c6.*a4)) + sigma(one).*sqrt(2).*(1./abs_detJt).*(Jtc.*(c4.*a2.*r ...
        + c5.*a1.*r + 2.*c5.*a2.*z + c5.*a3 + c6.*a2) + Jtd.*(c4.*a5.*r ...
        + c5.*a4.*r + 2.*c5.*a5.*z + c5.*a6 + c6.*a5));

    basis_div{efive,T} =@(r,z) sigma(two).*(1./abs_detJt).*(Jta.*(2.*b4.*a1.*r ...
        + b4.*a2.*z + b4.*(a3 - 1) + b5.*a1.*z + b6.*a1) ...
        + Jtb.*(2.*b4.*a4.*r + b4.*a5.*z + b4.*a6 + b5.*a4.*z + b6.*a4)) ...
        + sigma(two).*(1./abs_detJt).*(Jtc.*(b4.*a2.*r + b5.*a1.*r ...
        + 2.*b5.*a2.*z + b5.*(a3 - 1) + b6.*a2) + Jtd.*(b4.*a5.*r ...
        + b5.*a4.*r + 2.*b5.*a5.*z + b5.*a6 + b6.*a5));

    basis_div{esix,T} =@(r,z) sigma(three).*(1./abs_detJt).*(Jta.*(2.*c1.*a1.*r ...
        + c1.*a2.*z + c1.*a3 + c2.*a1.*z + c3.*a1) + Jtb.*(2.*c1.*a4.*r ...
        + c1.*a5.*z + c1.*(a6 - 1) + c2.*a4.*z + c3.*a4)) ...
        + sigma(three).*(1./abs_detJt).*(Jtc.*(c1.*a2.*r + c2.*a1.*r ...
        + 2.*c2.*a2.*z + c2.*a3 + c3.*a2) + Jtd.*(c1.*a5.*r + c2.*a4.*r ...
        + 2.*c2.*a5.*z + c2.*(a6 - 1) + c3.*a5));

    basis_div{7,T} =@(r,z) (1./abs_detJt).*(Jta.*(2.*a1.*a4.*r + a1.*a5.*z ...
        + a1.*a6 + a2.*a4.*z + a3.*a4) + Jtb.*(2.*a4.*a4.*r + a4.*a5.*z ...
        + a4.*(a6 - 1) + a5.*a4.*z + a6.*a4)) ...
        + (1./abs_detJt).*(Jtc.*(a1.*a5.*r + a2.*a4.*r + 2.*a2.*a5.*z ...
        + a2.*a6 + a3.*a5) + Jtd.*(a4.*a5.*r + a5.*a4.*r + 2.*a5.*a5.*z ...
        + a5.*(a6 - 1) + a6.*a5));
    
    basis_div{8,T} =@(r,z) (1./abs_detJt).*(Jta.*(2.*a1.*a1.*r + a1.*a2.*z ...
        + a1.*(a3 - 1) + a2.*a1.*z + a3.*a1) + Jtb.*(2.*a1.*a4.*r ...
        + a1.*a5.*z + a1.*a6 + a2.*a4.*z + a3.*a4)) ...
        + (1./abs_detJt).*(Jtc.*(a1.*a2.*r + a2.*a1.*r + 2.*a2.*a2.*z ...
        + a2.*(a3 - 1) + a3.*a2) + Jtd.*(a1.*a5.*r + a2.*a4.*r ...
        + 2.*a2.*a5.*z + a2.*a6 + a3.*a5));
end

end
