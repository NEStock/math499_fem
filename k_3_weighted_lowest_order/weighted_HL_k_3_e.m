function [err_z,err_p] = weighted_HL_k_3_e(f,gd,sf,ns,mesh,z_vec_r,z_vec_th,z_vec_z,p_exact,n)
%WEIGHTED_HL_K_3_E Hodge Laplacian k = 3 P1 Finite Element Method.
%   This program is set up to be given an exact solution.
%   Hodge Laplacian k = 3 case, P1
%   {psi_i}i=1->Ne+Nt is the basis for Ch
%   {chi_j}j=1->Nt is the basis for Dh
%   (Ch is the weighted fourier Raviart Thomas space)
%   (Dh is the piecewise constant space)
%   Solve for (z,p) in (Ch x Dh) s.t.
%       (z , w)_r - (p , div_rz^n(w))_r = 0
%       (div_rz^n(z) , s)_r = (f , s)_r
%           for all w in Ch, s in Dh
%
% Syntax:
%     [err] = weighted_HL_k_3_e(f_vec_r,f_vec_th,f_vec_z,gd,sf,ns,mesh,z_vec_r,z_vec_th,z_vec_z,p_exact,n)
%     f - given function
%     gd,sf,ns - outputs of pdepoly specifying domain
%     mesh - max mesh level
%     z_vec_r - exact solution z vector r component
%     z_vec_th - exact solution z vector theta component
%     z_vec_z - exact solution z vector z component
%     p_vec_r - exact solution p vector r component
%     p_vec_th - exact solution p vector theta component
%     p_vec_z - exact solution p vector z component
%
% Outputs:
%     err_z - array of L2 errors for mesh levels corresponding to indices
%     err_p - array of L2 errors for mesh levels corresponding to indices
%
% Usage Exampled:
%    addpath ../data ../data/weighted_HL_k_3/
%    n = 1;
%    [z_vec_r,z_vec_th,z_vec_z,p_exact,f] = get_data_1(n);
%    mesh = 7;
%    pdepoly([0,1,1,0], [0,0,1,1]);
%       (OR) [gd,sf,ns] = get_gd_sf_ns([0,1,1,0],[0,0,1,1]);
%    [err_z,err_p] = weighted_HL_k_3_e(f,gd,sf,ns,mesh,z_vec_r,z_vec_th,z_vec_z,p_exact,n)
% Dependencies:
%    find_edge_connectivity.m
%    basis_functions_weighted_fourier_rt.m
%    create_B_weighted_HL_k_3_p1.m
%    create_F_weighted_HL_k_3_p1.m
%    display_errors.m
%    errors_exact_weighted_HL_k_3_p1.m
%    mass_matrix_weighted_HL_k_3_p1.m
%
% Author: Nicole Stock
% Date: Fall 2020

addpath('../')
addpath('../raviart_thomas_weighted_fourier_lowest_order/')

model=createpde(1);
g=decsg(gd,sf,ns);
geometryFromEdges(model,g);
[p,e,t]=initmesh(g,'hmax',inf);
tr = triangulation(t(1:3,:)',p');
ed = edges(tr);
t_ed = find_edge_connectivity(t,ed);
%pdemesh(p,e,t, 'NodeLabels','on', 'ElementLabels','on');

% Init error vector
err_z = zeros(1,mesh);
err_p = zeros(1,mesh);

if mesh > 1
    % To ensure we refine every triangle the same
    [~,num_node]=size(p);
    it=zeros(1,num_node);
    for i=1:num_node
        it(i)=i;
    end   

    [basis_edges,basis_triangles,z_h,p_h] = solve(p,e,t,ed,t_ed,f,n);
    [err_z(1),err_p(1)] = errors_exact_weighted_HL_k_3_p1(p,t,t_ed,basis_edges,basis_triangles,z_h,z_vec_r,z_vec_th,z_vec_z,p_h,p_exact,n);
    
    for i = 2:mesh
        % Refine mesh to next level
        [p,e,t]=refinemesh(g,p,e,t,it,'regular');
        tr = triangulation(t(1:3,:)',p');
        ed = edges(tr);
        t_ed = find_edge_connectivity(t,ed,i);
        %pdemesh(p,e,t, 'NodeLabels','on', 'ElementLabels','on');

        [basis_edges,basis_triangles,z_h,p_h] = solve(p,e,t,ed,t_ed,f,n);
        [err_z(i),err_p(i)] = errors_exact_weighted_HL_k_3_p1(p,t,t_ed,basis_edges,basis_triangles,z_h,z_vec_r,z_vec_th,z_vec_z,p_h,p_exact,n);
    end
    fprintf('z\n');
    display_errors(err_z);
    fprintf('p\n');
    display_errors(err_p);

end
% mesh level must be greater than 1

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

    %figure();
    %pdeplot(p,e,t, 'XYData',x, 'ZData', x, 'Mesh', 'on');
end
