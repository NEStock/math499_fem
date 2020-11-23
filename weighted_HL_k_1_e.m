function [err_u,err_s] = weighted_HL_k_1_e(f_vec_r,f_vec_th,f_vec_z,gd,sf,ns,mesh,u_vec_r,u_vec_th,u_vec_z,s,n)
%WEIGHTED_HL_K_1_E Hodge Laplacian k = 1 P1 Finite Element Method.
%   This program is set up to be given an exact solution.
%
% Syntax:
%     [err] = weighted_HL_k_1_e(f_vec_r,f_vec_th,f_vec_z,gd,sf,ns,mesh,u_vec_r,u_vec_th,u_vec_z,s,n)
%     f_vec_r - given function r component
%     f_vec_th - given function theta component
%     f_vec_z - given function z component
%     gd,sf,ns - outputs of pdepoly specifying domain
%     mesh - max mesh level
%     u_vec_r - exact solution z vector r component
%     u_vec_th - exact solution z vector theta component
%     u_vec_z - exact solution z vector z component
%     s - exact solution s
%
% Outputs:
%     err_u - array of L2 errors for mesh levels corresponding to indices
%     err_s - array of L2 errors for mesh levels corresponding to indices
%
% Usage Exampled:
%    n = 1;
%    [u_vec_r,u_vec_th,u_vec_z,s,f_vec_r,f_vec_th,f_vec_z] = get_data_1(n)
%    mesh = 7;
%    pdepoly([0,1,1,0], [0,0,1,1]);
%       (OR) [gd,sf,ns] = get_gd_sf_ns([0,1,1,0],[0,0,1,1]);
%    [err_u,err_s] = weighted_HL_k_1_e(f_vec_r,f_vec_th,f_vec_z,gd,sf,ns,mesh,u_vec_r,u_vec_th,u_vec_z,s,n)
% Dependencies:
%    find_edge_connectivity.m
%    basis_functions_weighted_modified_nedelec_and_p1.m
%    create_B_weighted_HL_k_1_p1.m
%    create_F_weighted_HL_k_1_p1.m
%    create_S_weighted_HL_k_1_p1.m
%    display_errors.m
%    errors_exact_weighted_HL_k_1_p1.m
%    mass_matrix_weighted_HL_k_1_p1.m
%
% Author: Nicole Stock
% Date: Fall 2020

addpath('data')

model=createpde(1);
g=decsg(gd,sf,ns);
geometryFromEdges(model,g);
[p,e,t]=initmesh(g,'hmax',inf);
tr = triangulation(t(1:3,:)',p');
ed = edges(tr);
t_ed = find_edge_connectivity(t,ed);
%pdemesh(p,e,t, 'NodeLabels','on', 'ElementLabels','on');

% Init error vector
err_u = zeros(1,mesh);
err_s = zeros(1,mesh);

if mesh > 1
    % To ensure we refine every triangle the same
    [~,num_node]=size(p);
    it=zeros(1,num_node);
    for i=1:num_node
        it(i)=i;
    end   

    [basis_nodes,basis_edges,u_h,s_h] = solve(p,e,t,ed,t_ed,f_vec_r,f_vec_th,f_vec_z,n);
    [err_u(1),err_s(1)] = errors_exact_weighted_HL_k_1_p1(p,t,t_ed,basis_nodes,basis_edges,u_h,u_vec_r,u_vec_th,u_vec_z,s_h,s,n);
    
    for i = 2:mesh
        % Refine mesh to next level
        [p,e,t]=refinemesh(g,p,e,t,it,'regular');
        tr = triangulation(t(1:3,:)',p');
        ed = edges(tr);
        t_ed = find_edge_connectivity(t,ed,i);
        %pdemesh(p,e,t, 'NodeLabels','on', 'ElementLabels','on');

        [basis_nodes,basis_edges,u_h,s_h] = solve(p,e,t,ed,t_ed,f_vec_r,f_vec_th,f_vec_z,n);
        [err_u(i),err_s(i)] = errors_exact_weighted_HL_k_1_p1(p,t,t_ed,basis_nodes,basis_edges,u_h,u_vec_r,u_vec_th,u_vec_z,s_h,s,n);
    end
    fprintf('u\n');
    display_errors(err_u);
    fprintf('s\n');
    display_errors(err_s);

end
% mesh level must be greater than 1

% end main
end

% subfunction
function [basis_nodes,basis_edges,u_h,s_h] = solve(p,e,t,ed,t_ed,f_vec_r,f_vec_th,f_vec_z,n)
    [basis_nodes,basis_edges] = basis_functions_weighted_modified_nedelec_and_p1(p,t,ed,t_ed);
    M = mass_matrix_weighted_HL_k_1_p1(p,t,basis_nodes,n);
    B = create_B_weighted_HL_k_1_p1(p,t,ed,t_ed,basis_nodes,basis_edges,n);
    S = create_S_weighted_HL_k_1_p1(p,t,ed,t_ed,basis_nodes,basis_edges,n);
    F = create_F_weighted_HL_k_1_p1(p,t,ed,t_ed,basis_nodes,basis_edges,f_vec_r,f_vec_th,f_vec_z,n);
    
    Bt = B.';
    MB = (M\B);
    u_h = (Bt*MB + S)\F;
    s_h = MB*u_h;
end
