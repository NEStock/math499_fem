function [err,grad_err,max_err] = errors_no_exact_weighted_p2(p,t,p2,t2,basis,u_h_km1,u_h_k,k)
% ERRORS_NO_EXACT_WEIGHTED_p2 - Calculate the errors of our solution u_h_km1
% compared to the approximate solution for the next mesh level (u_h_k).
%
% Syntax:
%     [err,grad_err,max_err] = errors_no_exact_weighted_p2(p,e,t,u_h_km1,u_h_k)
%
% Inputs:
%     p - a 2xNumNodes matrix representing nodal coordinates.
%     t - a 4xNumTriangles matrix representing the element connectivity in 
%         terms of node IDs. The end row of T represents the geometry face ID 
%         to which the element belongs
%     basis - a 3x3xNumTriangles matrix representing piece-wise basis 
%         functions for each node in each triangle. basis(i,:,k) represents 
%         the pieceiwise basis function for the ith node in triangle k. 
%     u_h_km1 - approximate solution for mesh level k-1
%     u_h_k - approximate solution for mesh level k
%     k - given weight
%
% Outputs:
%    err - L2 error
%    grad_err - L2 gradient error
%    max_err - max error
%
% Author: Nicole Stock
% Date: Spring 2020

[~,triangles] = size(t);
[~,nodes] = size(p);
[~,mid_nodes] = size(p2);

integral = 0;
grad_integral = 0;

for T = 1:triangles
    
    % get coordinates of triangle T
    coordinates = zeros(3,2);
    for n = 1:3
        node = t(n,T);
        % get x,y coordinates
        coordinates(n,:) = p(:,node);
    end
    
    [R,Z,Wr,Wz] = triquad(7, coordinates);

    b1 = basis(:,1,T);
    b2 = basis(:,2,T);
    b3 = basis(:,3,T);
    b4 = basis(:,4,T);
    b5 = basis(:,5,T);
    b6 = basis(:,6,T);
    
    n1 = t(1,T);
    n2 = t(2,T);
    n3 = t(3,T);
    n4 = t2(1,T);
    n5 = t2(2,T);
    n6 = t2(3,T);
    
    if k == 0
        % find L2 Error for u_h_km1 & u_h_k
        approx_km1 =@(r,z) u_h_km1(n1).*(b1(1).*r.^2 + b1(2).*r.*z ...
            + b1(3).*z.^2 + b1(4).*r + b1(5).*z + b1(6)) ...
            + u_h_km1(n2).*(b2(1).*r.^2 + b2(2).*r.*z ...
            + b2(3).*z.^2 + b2(4).*r + b2(5).*z + b2(6)) ...
            + u_h_km1(n3).*(b3(1).*r.^2 + b3(2).*r.*z ...
            + b3(3).*z.^2 + b3(4).*r + b3(5).*z + b3(6)) ...
            + u_h_km1(n4 + nodes).*(b4(1).*r.^2 + b4(2).*r.*z ...
            + b4(3).*z.^2 + b4(4).*r + b4(5).*z + b4(6)) ...
            + u_h_km1(n5 + nodes).*(b5(1).*r.^2 + b5(2).*r.*z ...
            + b5(3).*z.^2 + b5(4).*r + b5(5).*z + b5(6)) ...
            + u_h_km1(n6 + nodes).*(b6(1).*r.^2 + b6(2).*r.*z ...
            + b6(3).*z.^2 + b6(4).*r + b6(5).*z + b6(6));
        approx_k =@(r,z) u_h_k(n1).*(b1(1).*r.^2 + b1(2).*r.*z ...
            + b1(3).*z.^2 + b1(4).*r + b1(5).*z + b1(6)) ...
            + u_h_k(n2).*(b2(1).*r.^2 + b2(2).*r.*z ...
            + b2(3).*z.^2 + b2(4).*r + b2(5).*z + b2(6)) ...
            + u_h_k(n3).*(b3(1).*r.^2 + b3(2).*r.*z ...
            + b3(3).*z.^2 + b3(4).*r + b3(5).*z + b3(6)) ...
            + u_h_k(n4 + nodes).*(b4(1).*r.^2 + b4(2).*r.*z ...
            + b4(3).*z.^2 + b4(4).*r + b4(5).*z + b4(6)) ...
            + u_h_k(n5 + nodes).*(b5(1).*r.^2 + b5(2).*r.*z ...
            + b5(3).*z.^2 + b5(4).*r + b5(5).*z + b5(6)) ...
            + u_h_k(n6 + nodes).*(b6(1).*r.^2 + b6(2).*r.*z ...
            + b6(3).*z.^2 + b6(4).*r + b6(5).*z + b6(6));       
    else
        % find L2 Error for u_h_km1 & u_h_k
        approx_km1 =@(r,z) u_h_km1(n1).*(r./k).*(b1(1).*r.^2 + b1(2).*r.*z ...
            + b1(3).*z.^2 + b1(4).*r + b1(5).*z + b1(6)) ...
            + u_h_km1(n2).*(r./k).*(b2(1).*r.^2 + b2(2).*r.*z ...
            + b2(3).*z.^2 + b2(4).*r + b2(5).*z + b2(6)) ...
            + u_h_km1(n3).*(r./k).*(b3(1).*r.^2 + b3(2).*r.*z ...
            + b3(3).*z.^2 + b3(4).*r + b3(5).*z + b3(6)) ...
            + u_h_km1(n4 + nodes).*(r./k).*(b4(1).*r.^2 + b4(2).*r.*z ...
            + b4(3).*z.^2 + b4(4).*r + b4(5).*z + b4(6)) ...
            + u_h_km1(n5 + nodes).*(r./k).*(b5(1).*r.^2 + b5(2).*r.*z ...
            + b5(3).*z.^2 + b5(4).*r + b5(5).*z + b5(6)) ...
            + u_h_km1(n6 + nodes).*(r./k).*(b6(1).*r.^2 + b6(2).*r.*z ...
            + b6(3).*z.^2 + b6(4).*r + b6(5).*z + b6(6));
        approx_k =@(r,z) u_h_k(n1).*(r./k).*(b1(1).*r.^2 + b1(2).*r.*z ...
            + b1(3).*z.^2 + b1(4).*r + b1(5).*z + b1(6)) ...
            + u_h_k(n2).*(r./k).*(b2(1).*r.^2 + b2(2).*r.*z ...
            + b2(3).*z.^2 + b2(4).*r + b2(5).*z + b2(6)) ...
            + u_h_k(n3).*(r./k).*(b3(1).*r.^2 + b3(2).*r.*z ...
            + b3(3).*z.^2 + b3(4).*r + b3(5).*z + b3(6)) ...
            + u_h_k(n4 + nodes).*(r./k).*(b4(1).*r.^2 + b4(2).*r.*z ...
            + b4(3).*z.^2 + b4(4).*r + b4(5).*z + b4(6)) ...
            + u_h_k(n5 + nodes).*(r./k).*(b5(1).*r.^2 + b5(2).*r.*z ...
            + b5(3).*z.^2 + b5(4).*r + b5(5).*z + b5(6)) ...
            + u_h_k(n6 + nodes).*(r./k).*(b6(1).*r.^2 + b6(2).*r.*z ...
            + b6(3).*z.^2 + b6(4).*r + b6(5).*z + b6(6)); 
    end
    
    integrand =@(r,z) ((approx_km1(r,z) - approx_k(r,z)).^2).*r;
    
    integral = integral + Wr'*feval(integrand,R,Z)*Wz;
    
    % find Grad L2 Error for u_h_km1 & u_h_k
    grad_approx_r_km1 =@(r,z) u_h_km1(n1).*(2.*b1(1).*r + b1(2).*z ...
        + b1(4)) + u_h_km1(n2).*(2.*b2(1).*r + b2(2).*z + b2(4)) ...
        + u_h_km1(n3).*(2.*b3(1).*r + b3(2).*z + b3(4)) ...
        + u_h_km1(n4 + nodes).*(2.*b4(1).*r + b4(2).*z + b4(4)) ...
        + u_h_km1(n5 + nodes).*(2.*b5(1).*r + b5(2).*z + b5(4)) ...
        + u_h_km1(n6 + nodes).*(2.*b6(1).*r + b6(2).*z + b6(4));
    grad_approx_z_km1 =@(r,z) u_h_km1(n1).*(b1(2).*r + 2.*b1(3).*z ...
        + b1(5)) + u_h_km1(n2).*(b2(2).*r + 2.*b2(3).*z + b2(5)) ...
        + u_h_km1(n3).*(b3(2).*r + 2.*b3(3).*z + b3(5)) ...
        + u_h_km1(n4 + nodes).*(b4(2).*r + 2.*b4(3).*z + b4(5)) ...
        + u_h_km1(n5 + nodes).*(b5(2).*r + 2.*b5(3).*z + b5(5)) ...
        + u_h_km1(n6 + nodes).*(b6(2).*r + 2.*b6(3).*z + b6(5));
    
    grad_approx_r_k =@(r,z) u_h_k(n1).*(2.*b1(1).*r + b1(2).*z ...
        + b1(4)) + u_h_k(n2).*(2.*b2(1).*r + b2(2).*z + b2(4)) ...
        + u_h_k(n3).*(2.*b3(1).*r + b3(2).*z + b3(4)) ...
        + u_h_k(n4 + nodes).*(2.*b4(1).*r + b4(2).*z + b4(4)) ...
        + u_h_k(n5 + nodes).*(2.*b5(1).*r + b5(2).*z + b5(4)) ...
        + u_h_k(n6 + nodes).*(2.*b6(1).*r + b6(2).*z + b6(4));
    grad_approx_z_k =@(r,z) u_h_k(n1).*(b1(2).*r + 2.*b1(3).*z ...
        + b1(5)) + u_h_k(n2).*(b2(2).*r + 2.*b2(3).*z + b2(5)) ...
        + u_h_k(n3).*(b3(2).*r + 2.*b3(3).*z + b3(5)) ...
        + u_h_k(n4 + nodes).*(b4(2).*r + 2.*b4(3).*z + b4(5)) ...
        + u_h_k(n5 + nodes).*(b5(2).*r + 2.*b5(3).*z + b5(5)) ...
        + u_h_k(n6 + nodes).*(b6(2).*r + 2.*b6(3).*z + b6(5));
    
    grad_integrand =@(r,z) ((grad_approx_r_km1(r,z) ...
        - grad_approx_r_k(r,z)).^2 + (grad_approx_z_km1(r,z) ...
        - grad_approx_z_k(r,z)).^2).*r;
    
    grad_integral = grad_integral + Wr'*feval(grad_integrand,R,Z)*Wz;
end

% find Max Error
max_err = -Inf;
for n = 1:nodes+mid_nodes
    diff = abs(u_h_km1(n) - u_h_k(n));
    
    if diff > max_err
        max_err = diff;
    end
end

err = sqrt(integral);
grad_err = sqrt(grad_integral);

% end