function [err,grad_err,max_err] = weighted_nedelec_e(f_vec_r,f_vec_z,gd,sf,ns,mesh_level,u_vec_r,u_vec_z)
%WEIGHTED_NEDELEC_E Weighted . 
%   This program is set up to be given an exact solution.
%   Lowest Order Nedelec Space
%
% Syntax:
%     [err,grad_err,max_err] = 
%       weighted_nedelec_e(f_vec_r,f_vec_z,gd,sf,ns,mesh_level,u_vec_r,u_vec_z)
%
% Inputs:
%     f_vec - given vector
%     gd,sf,ns - outputs of pdepoly specifying domain
%     mesh_level - max mesh level
%     u - exact solution vector
%
% Outputs:
%     err - array of L2 errors for mesh levels corresponding to indices
%     grad_err - array of L2 gradient errors for mesh levels  
%         corresponding to indices
%     max_err - array of max errors for mesh levels corresponding to
%         indicies
%
% Usage Exampled:
%    [u_vec_r,u_vec_z,f_vec_r,f_vec_z] = get_data_1();
%    mesh = 8;
%    pdepoly([0,1,1,0], [0,0,1,1]);
%       (OR) [gd,sf,ns] = get_gd_sf_ns([0,1,1,0],[0,0,1,1]);
%    [err,grad_err,max_err] = weighted_nedelec_e(f_vec_r,f_vec_z,gd,sf,ns,mesh_level,u_vec_r,u_vec_z)
% Dependencies:
%    basis_functions_weighted_nedelec.m
%    create_b_nedelec.m
%    display_errors.m
%    errors_exact_weighted_p2.m
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

% Init error vectors
err = zeros(1,mesh_level);
grad_err = zeros(1,mesh_level);
max_err = zeros(1,mesh_level);

if mesh_level > 1
    % To ensure we refine every triangle the same
    [~,num_node]=size(p);
    it=zeros(1,num_node);
    for i=1:num_node
        it(i)=i;
    end   
        
    [basis,x] = solve(p,e,t,ed,t_ed,f_vec_r,f_vec_z);
    %[err(1),grad_err(1),max_err(1)] = errors_exact_weighted(p,t,basis,x,u_vec);

    for i = 2:mesh_level
        % Refine mesh to next level
        [p,e,t]=refinemesh(g,p,e,t,it,'regular');
        tr = triangulation(t(1:3,:)',p');
        ed = edges(tr);
        t_ed = find_edges(t,ed);
        %pdemesh(p,e,t, 'NodeLabels','on', 'ElementLabels','on');

        [basis,x] = solve(p,e,t,ed,t_ed,f_vec_r,f_vec_z);
        %[err(i),grad_err(i),max_err(i)] = errors_exact_weighted(p,t,basis,x,u_vec_r,u_vec_z);
    end
    display_errors(err,grad_err,max_err)

end
% mesh level must be greater than 1

% end main
end

% subfunction
function [basis,x] = solve(p,e,t,ed,t_ed,f_vec_r,f_vec_z)
    basis = basis_functions_weighted_nedelec(p,t,ed,t_ed);
    A = stiffness_matrix_weighted_nedelec(p,t,ed,t_ed,basis);
    b = create_b_nedelec(p,t,ed,t_ed,basis,f_vec_r,f_vec_z);
    x = A\b;
    disp(x)

    %figure();
    %pdeplot(p,e,t, 'XYData',x, 'ZData', x, 'Mesh', 'on');
end
