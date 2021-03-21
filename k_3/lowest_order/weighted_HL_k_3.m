function [basis_edges,basis_triangles,z_h,p_h] = weighted_HL_k_3(f,gd,sf,ns,mesh,n)
%WEIGHTED_HL_K_3 Hodge Laplacian k = 3 lowest order Finite Element Method.
%   This program is set up to give approximations of the unknown solutions 
%   z and p.
%   Hodge Laplacian k = 3 case, lowest order
%   {psi_i}i=1->Ne+Nt is the basis for Ch0
%   {chi_j}j=1->Nt is the basis for Dh0
%   (Ch is the weighted fourier Raviart Thomas space)
%   (Dh is the piecewise constant space)
%   Solve for (z,p) in (Ch x Dh) s.t.
%       (z , w)_r - (p , div_rz^n(w))_r = 0
%       (div_rz^n(z) , s)_r = (f , s)_r
%           for all w in Ch, s in Dh
%
% Syntax:
%     [basis_edges,basis_triangles,z_h,p_h] 
%       = weighted_HL_k_3(f,gd,sf,ns,mesh,n)
%     f - given function
%     gd,sf,ns - outputs of pdepoly specifying domain
%     mesh - max mesh level
%     n - Fourier mode
%
% Outputs:
%     basis_edges - a matrix representing piece-wise basis functions for 
%         each edge in each triangle. basis(i,:,T) represents the 
%         pieceiwise basis function for the ith edge in triangle T.
%     basis_triangles - a vector representing piece-wise basis functions
%         for edge triangle. basis(1,T) represents the piecewise basis
%         function for the Tth triangle.
%     z_h - approximated solution for z
%     p_h - approximated solution for p
%
% Usage Exampled:
%    addpath ../../data ../data/
%    n = 1;
%    [z_vec_r,z_vec_th,z_vec_z,p_exact,f] = get_data_1(n);
%    mesh = 7;
%    pdepoly([0,1,1,0], [0,0,1,1]);
%       (OR) [gd,sf,ns] = get_gd_sf_ns([0,1,1,0],[0,0,1,1]);
%    [basis_edges,basis_triangles,z_h,p_h] = weighted_HL_k_3(f,gd,sf,ns,mesh,n)
% Dependencies:
%    find_edge_connectivity.m
%    basis_functions_weighted_fourier_rt.m
%    create_B_weighted_HL_k_3_p1.m
%    create_F_weighted_HL_k_3_p1.m
%    mass_matrix_weighted_HL_k_3_p1.m
%
% Author: Nicole Stock
% Date: Spring 2021

addpath('../../');
addpath('../../raviart_thomas/weighted_lowest_order/');
addpath('../../raviart_thomas/lowest_order/');

model=createpde(1);
g=decsg(gd,sf,ns);
geometryFromEdges(model,g);
[p,e,t]=initmesh(g,'hmax',inf);
%pdemesh(p,e,t, 'NodeLabels','on', 'ElementLabels','on');

% To ensure we refine every triangle the same
[~,num_node]=size(p);
it=zeros(1,num_node);
for i=1:num_node
    it(i)=i;
end   

for i = 2:mesh
    % Refine mesh to next level
    [p,e,t]=refinemesh(g,p,e,t,it,'regular');

    %pdemesh(p,e,t, 'NodeLabels','on', 'ElementLabels','on');
end
tr = triangulation(t(1:3,:)',p');
ed = edges(tr);
t_ed = find_edge_connectivity(t,ed,mesh);

% solve
[basis_edges,basis_triangles,z_h,p_h] = solve(p,e,t,ed,t_ed,f,n);

% end main
end

% subfunction
function [basis_edges,basis_triangles,z_h,p_h] = solve(p,e,t,ed,t_ed,f,n)
    [basis_edges, basis_triangles] = basis_functions_weighted_fourier_rt(p,t,ed,t_ed);
    M = mass_matrix_weighted_HL_k_3_p1(p,t,ed,t_ed,basis_edges,basis_triangles,n);
    B = create_B_weighted_HL_k_3_p1(p,t,ed,t_ed,basis_edges,basis_triangles);
    F = create_F_weighted_HL_k_3_p1(p,t,f);
    
    Bt = B.';
    MB = (M\B);
    p_h = (Bt*MB)\F;
    z_h = MB*p_h;
end
