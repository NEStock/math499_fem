function [basis_nodes,basis_NP1_edges,basis_RT_edges,basis_triangles,u_h,s_h] = weighted_HL_k_2(f_vec_r,f_vec_th,f_vec_z,gd,sf,ns,mesh,n)
%WEIGHTED_HL_K_2 Hodge Laplacian k = 2 lowest order Finite Element Method.
%   This program is set up to give approximations of the unknown solutions 
%   u and s.
%   Hodge Laplacian k = 2 case, lowest order
%   {zeta_j}j=1->N+Ne is the basis for Bh0
%   {psi_i}i=1->Ne+Nt is the basis for Ch0
%   (Bh is the weighted fourier modified Nedelec and P1 space)
%   (Ch is the weighted fourier Raviart Thomas space)
%   Solve for (s,u) in (Bh x Ch) s.t.
%       (s , w)_r - (curl_rz^n(w) , u)_r = 0
%       (curl_rz^n(s) , v)_r + (div_rz^n(u) , div_rz^n(v))_r = (f , v)_r
%           for all w in Bh, v in Ch
%
% Syntax:
%     [err] = weighted_HL_k_2(f_vec_r,f_vec_th,f_vec_z,gd,sf,ns,mesh,n)
%     f_vec_r - given function r component
%     f_vec_th - given function theta component
%     f_vec_z - given function z component
%     gd,sf,ns - outputs of pdepoly specifying domain
%     mesh - max mesh level
%     n - Fourier mode
%
% Outputs:
%     basis_nodes - a matrix representing piece-wise basis functions for
%         each node in each triangle. basis(i,:,T) represents the
%         pieceiwise basis function for the ith node in triangle T.
%     basis_NP1_edges - a matrix representing piece-wise basis functions
%         for each edge in each triangle for the weighted fourier Nedelec
%         and P1 space. basis(i,:,T) represents the pieceiwise basis
%         function for the ith edge in triangle T.
%     basis_RT_edges - a matrix representing piece-wise basis functions
%         for each edge in each triangle for the weighted fourier Raviart 
%         Thomas space. basis(i,:,T) represents the pieceiwise basis 
%         function for the ith edge in triangle T.
%     basis_triangles - a vector representing piece-wise basis functions
%         for edge triangle. basis(1,T) represents the piecewise basis
%         function for the Tth triangle.
%     u_h - approximated solution for u
%     s_h - approximated solution for s
%
% Usage Exampled:
%    addpath ../../data ../data/
%    n = 1;
%    [~,~,~,~,~,~,f_vec_r,f_vec_th,f_vec_z] = get_data_1(n);
%    mesh = 7;
%    pdepoly([0,1,1,0], [0,0,1,1]);
%       (OR) [gd,sf,ns] = get_gd_sf_ns([0,1,1,0],[0,0,1,1]);
%    [basis_nodes,basis_NP1_edges,basis_RT_edges,basis_triangles,u_h,s_h] 
%       = weighted_HL_k_2(f_vec_r,f_vec_th,f_vec_z,gd,sf,ns,mesh,n);
% Dependencies:
%    find_edge_connectivity.m
%    basis_functions_weighted_HL_k_2_p1.m
%    create_B_weighted_HL_k_2_p1.m
%    create_F_weighted_HL_k_2_p1.m
%    create_S_weighted_HL_k_2_p1.m
%    mass_matrix_weighted_HL_k_2_p1.m
%
% Author: Nicole Stock
% Date: Spring 2021

addpath('../../')
addpath('../../k_0/p1/')
addpath('../../nedelec/weighted_lowest_order/')
addpath('../../raviart_thomas/weighted_lowest_order/');
addpath('../../raviart_thomas/lowest_order/');
addpath('../../modified_weighted_nedelec_pk/p1/');

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
[basis_nodes,basis_NP1_edges,basis_RT_edges,basis_triangles,u_h,s_h]...
    = solve(p,e,t,ed,t_ed,f_vec_r,f_vec_th,f_vec_z,n);

% end main
end

% subfunction
function [basis_nodes,basis_NP1_edges,basis_RT_edges,basis_triangles,u_h,s_h] = solve(p,e,t,ed,t_ed,f_vec_r,f_vec_th,f_vec_z,n)
    [basis_nodes,basis_NP1_edges,basis_RT_edges,basis_triangles] = basis_functions_weighted_HL_k_2_p1(p,t,ed,t_ed);
    M = mass_matrix_weighted_HL_k_2_p1(p,t,ed,t_ed,basis_nodes,basis_NP1_edges,n);
    B = create_B_weighted_HL_k_2_p1(p,t,ed,t_ed,basis_nodes,basis_NP1_edges,basis_RT_edges,basis_triangles,n);
    S = create_S_weighted_HL_k_2_p1(p,t,ed,t_ed,basis_RT_edges,basis_triangles);
    F = create_F_weighted_HL_k_2_p1(p,t,ed,t_ed,basis_RT_edges,basis_triangles,f_vec_r,f_vec_th,f_vec_z,n);
    
    Bt = B.';
    MB = (M\B);
    u_h = (Bt*MB + S)\F;
    s_h = MB*u_h;
end
