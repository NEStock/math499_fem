function [err] = raviart_thomas_0_e(f_vec_r,f_vec_z,gd,sf,ns,mesh,u_vec_r,u_vec_z)
%RAVIART_THOMAS_0_E Raviart Thomas-1 Space Finite Element Method. 
%   This program is set up to be given an exact solution.
%   Lowest Order Raviart Thomas Space
%
% Syntax:
%     [err] = raviart_thomas_0_e(f_vec_r,f_vec_z,gd,sf,ns,mesh,u_vec_r,u_vec_z)
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
%    addpath ../../data/ ../data/unweighted2D/
%    [u_vec_r,u_vec_z,f_vec_r,f_vec_z] = get_data_1();
%    mesh = 8;
%    pdepoly([0,1,1,0], [0,0,1,1]);
%       (OR) [gd,sf,ns] = get_gd_sf_ns([0,1,1,0],[0,0,1,1]);
%    [err] = raviart_thomas_0_e(f_vec_r,f_vec_z,gd,sf,ns,mesh,u_vec_r,u_vec_z)
% Dependencies:
%    find_edge_connectivity.m
%    basis_functions_rt1.m
%    create_b_rt1.m
%    display_errors.m
%    errors_exact_rt1.m
%    stiffness_matrix_rt1.m
%
% Author: Nicole Stock
% Date: Fall 2020

addpath('../../')

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

    [basis,x] = solve(p,e,t,ed,t_ed,f_vec_r,f_vec_z);
    err(1) = errors_exact_rt0(p,t,ed,t_ed,basis,x,u_vec_r,u_vec_z);

    for i = 2:mesh
        % Refine mesh to next level
        [p,e,t]=refinemesh(g,p,e,t,it,'regular');
        tr = triangulation(t(1:3,:)',p');
        ed = edges(tr);
        t_ed = find_edge_connectivity(t,ed,i);
        %pdemesh(p,e,t, 'NodeLabels','on', 'ElementLabels','on');

        [basis,x] = solve(p,e,t,ed,t_ed,f_vec_r,f_vec_z);
        err(i) = errors_exact_rt0(p,t,ed,t_ed,basis,x,u_vec_r,u_vec_z);
    end
    display_errors(err);

end
% mesh level must be greater than 1

% end main
end

% subfunction
function [basis,x] = solve(p,e,t,ed,t_ed,f_vec_r,f_vec_z)
    [basis, basis_div] = basis_functions_rt0(p,t,ed,t_ed);
    S = stiffness_matrix_rt0(p,t,ed,t_ed,basis,basis_div);
    b = create_b_rt0(p,t,ed,t_ed,basis,f_vec_r,f_vec_z);
    x = S\b;
    
    %disp(S)
    %disp(b)
    %figure();
    %pdeplot(p,e,t, 'XYData',x, 'ZData', x, 'Mesh', 'on');
end
