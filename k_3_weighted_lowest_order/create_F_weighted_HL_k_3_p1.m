function F = create_F_weighted_HL_k_3_p1(p,t,f)
% CREATE_F_WEIGHTED_HL_K_3_P1 - Create mass matrix
%   Hodge Laplacian k = 3 case, P1
%   (f, chi_i)_r where {chi_j}j=1->Nt is the basis for Dh
%   (Dh is the weighted fourier Raviart Thomas space)
%
% Syntax:
%     B = create_F_weighted_HL_k_3_p1(p,t,ed,t_ed,basis)
%
% Inputs:
%     p - a 2xNumNodes matrix representing nodal coordinates.
%     t - a 4xNumTriangles matrix representing the element connectivity in
%         terms of node IDs. The end row of T represents the geometry face 
%         ID to which the element belongs.
%     f - given function
%
% Outputs:
%     F - F matrix used to solve system of equations to approximate
%         solution
%
% Author: Nicole Stock
% Date: Fall 2020

[~,triangles] = size(t);
i_vec = zeros(1,triangles);
j_vec = ones(1,triangles);
s_vec = zeros(1,triangles);
index = 1;

for T = 1:triangles
    
    % get coordinates of triangle T
    coordinates = zeros(3,2);
    for N = 1:3
        node = t(N,T);
        % get x,y coordinates of triangle
        coordinates(N,:) = p(:,node);
    end
        
    [R,Z,Wr,Wz] = triquad(7, coordinates);
   
    integrand =@(r,z) (f(r,z)).*r;

    Q = Wr'*feval(integrand,R,Z)*Wz;               

    i_vec(index) = T;
    s_vec(index) = Q;
    index = index + 1;
end

F = sparse(i_vec,j_vec,s_vec,triangles,1);

% end