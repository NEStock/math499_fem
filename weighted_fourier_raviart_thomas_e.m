function [err] = weighted_fourier_raviart_thomas_e(f_vec_r,f_vec_th,f_vec_z,gd,sf,ns,mesh,u_vec_r,u_vec_th,u_vec_z,n)
%WEIGHTED_FOURIER_RAVIART_THOMAS_E Fourier Raviart Thomas Space Finite 
%  Element Method. 
%   This program is set up to be given an exact solution.
%   Lowest Order Weighted Fourier Raviart Thomas Space
%
% Syntax:
%     [err] = weighted_fourier_raviart_thomas_e(f_vec_r,f_vec_th,f_vec_z,gd,sf,ns,mesh,u_vec_r,u_vec_th,u_vec_z,n)
%     f_vec_r - given vector r component
%     f_vec_z - given vector z component
%     gd,sf,ns - outputs of pdepoly specifying domain
%     mesh - max mesh level
%     u_vec_r - exact solution vector r component
%     u_vec_z - exact solution vector z component
%
% Outputs:
%     err - array of L2 errors for mesh levels corresponding to indices
%
% Usage Exampled:
%    n = 1;
%    [u_vec_r,u_vec_th,u_vec_z,f_vec_r,f_vec_th,f_vec_z] = get_data_1(n);
%    mesh = 7;
%    pdepoly([0,1,1,0], [0,0,1,1]);
%       (OR) [gd,sf,ns] = get_gd_sf_ns([0,1,1,0],[0,0,1,1]);
%    [err] = weighted_fourier_raviart_thomas_e(f_vec_r,f_vec_th,f_vec_z,gd,sf,ns,mesh,u_vec_r,u_vec_th,u_vec_z,n)
% Dependencies:
%    find_edge_connectivity.m
%    basis_functions_weighted_fourier_rt.m
%    create_b_fourier_rt.m
%    display_errors.m
%    errors_exact_weighted_fourier_rt.m
%    stiffness_matrix_weighted_fourier_rt.m
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
err = zeros(1,mesh);

if mesh > 1
    % To ensure we refine every triangle the same
    [~,num_node]=size(p);
    it=zeros(1,num_node);
    for i=1:num_node
        it(i)=i;
    end   

    [basis_edges,basis_triangles,x] = solve(p,e,t,ed,t_ed,f_vec_r,f_vec_th,f_vec_z,n);
    err(1) = errors_exact_weighted_fourier_rt(p,t,ed,t_ed,basis_edges,basis_triangles,x,u_vec_r,u_vec_th,u_vec_z,n);
    
    for i = 2:mesh
        % Refine mesh to next level
        [p,e,t]=refinemesh(g,p,e,t,it,'regular');
        tr = triangulation(t(1:3,:)',p');
        ed = edges(tr);
        t_ed = find_edge_connectivity(t,ed,i);
        %pdemesh(p,e,t, 'NodeLabels','on', 'ElementLabels','on');

        [basis_edges,basis_triangles,x] = solve(p,e,t,ed,t_ed,f_vec_r,f_vec_th,f_vec_z,n);
        err(i) = errors_exact_weighted_fourier_rt(p,t,ed,t_ed,basis_edges,basis_triangles,x,u_vec_r,u_vec_th,u_vec_z,n);
    end
    display_errors(err, nan(1,mesh), nan(1,mesh));

end
% mesh level must be greater than 1

% end main
end

% subfunction
function [basis_edges,basis_triangles,x] = solve(p,e,t,ed,t_ed,f_vec_r,f_vec_th,f_vec_z,n)
    [basis_edges, basis_triangles] = basis_functions_weighted_fourier_rt(p,t,ed,t_ed);
    S = stiffness_matrix_weighted_fourier_rt(p,t,ed,t_ed,basis_edges,basis_triangles,n);
    b = create_b_weighted_fourier_rt(p,t,ed,t_ed,basis_edges,basis_triangles,f_vec_r,f_vec_th,f_vec_z,n);
    x = S\b;

    %figure();
    %pdeplot(p,e,t, 'XYData',x, 'ZData', x, 'Mesh', 'on');
end
