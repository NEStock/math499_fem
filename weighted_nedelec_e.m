function [err] = weighted_nedelec_e(f_vec_r,f_vec_z,gd,sf,ns,mesh,u_vec_r,u_vec_z)
%WEIGHTED_NEDELEC_E Weighted . 
%   This program is set up to be given an exact solution.
%   Lowest Order Nedelec Space
%
% Syntax:
%     [err] = weighted_nedelec_e(f_vec_r,f_vec_z,gd,sf,ns,mesh,u_vec_r,u_vec_z)
%
% Inputs:
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
%    [u_vec_r,u_vec_z,f_vec_r,f_vec_z] = get_data_1();
%    mesh = 8;
%    pdepoly([0,1,1,0], [0,0,1,1]);
%       (OR) [gd,sf,ns] = get_gd_sf_ns([0,1,1,0],[0,0,1,1]);
%    [err] = weighted_nedelec_e(f_vec_r,f_vec_z,gd,sf,ns,mesh,u_vec_r,u_vec_z)
% Dependencies:
%    find_edges.m
%    basis_functions_weighted_nedelec.m
%    create_b_nedelec.m
%    display_errors.m
%    errors_exact_weighted_nedelec.m
%    stiffness_matrix_weighted_nedelec.m
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
t_ed = find_edges(t,ed);
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

    [basis,x] = solve(p,e,t,ed,t_ed,f_vec_r,f_vec_z);
    err(1) = errors_exact_weighted_nedelec(p,t,t_ed,basis,x,u_vec_r,u_vec_z);

    for i = 2:mesh
        % Refine mesh to next level
        [p,e,t]=refinemesh(g,p,e,t,it,'regular');
        tr = triangulation(t(1:3,:)',p');
        ed = edges(tr);
        t_ed = find_edges(t,ed);
        %pdemesh(p,e,t, 'NodeLabels','on', 'ElementLabels','on');

        [basis,x] = solve(p,e,t,ed,t_ed,f_vec_r,f_vec_z);
        err(i) = errors_exact_weighted_nedelec(p,t,t_ed,basis,x,u_vec_r,u_vec_z);
    end
    display_errors(err, nan(1,mesh), nan(1,mesh));

end
% mesh level must be greater than 1

% end main
end

% subfunction
function [basis,x] = solve(p,e,t,ed,t_ed,f_vec_r,f_vec_z)
    basis = basis_functions_weighted_nedelec(p,ed,t_ed);
    S = stiffness_matrix_weighted_nedelec(p,t,ed,t_ed,basis);
    b = create_b_nedelec(p,t,ed,t_ed,basis,f_vec_r,f_vec_z);
    x = S\b;

    figure();
    pdeplot(p,e,t, 'XYData',x, 'ZData', x, 'Mesh', 'on');
end
