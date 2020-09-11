function [err,grad_err,max_err] = weighted_HL_k_0_e(u,grad_u_r,grad_u_z,gd,sf,ns,mesh_level,n)
%WEIGHTED_HL_K_0_e Weighted Hodge Laplacian with k = 0. 
%   This program is set up to be given an exact solution.
%
% Syntax:
%     [err,grad_err,max_err] = 
%       weighted_HL_k_0_e(u,grad_u_r,grad_u_z,gd,sf,ns,mesh_level,n)
%
% Inputs:
%     u - given function
%     gd,sf,ns - outputs of pdepoly specifying domain
%     mesh_level - max mesh level
%     n - Hodge Laplacian on Axisymmetrix Domain and its Discretization
%     weight
%
% Outputs:
%     err - array of L2 errors for mesh levels corresponding to indices
%     grad_err - array of L2 gradient errors for mesh levels  
%         corresponding to indices
%     max_err - array of max errors for mesh levels corresponding to
%         indicies
%
% Usage Exampled:
%    u = @(r,z) r^(1/2);
%    mesh = 8;
%    n = 1;
%    pdepoly([0,1,1,0], [0,0,1,1]);
%       (OR) [gd,sf,ns] = get_gd_sf_ns([0,1,1,0],[0,0,1,1])
%    [err,grad_err,max_err] = weighted_HL_k_0_e(u,grad_u_r,grad_u_z,gd,sf,ns,mesh,n)
% Dependencies:
%    basis_functions_weighted_p2.m
%    display_errors.m
%    errors_exact_weighted_p2.m
%    stiffness_matrix_weighted_p2.m
%
% Author: Nicole Stock
% Date: Fall 2020
addpath('data')

model=createpde(1);
g=decsg(gd,sf,ns);
geometryFromEdges(model,g);
[p,e,t]=initmesh(g,'hmax',inf);
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

    [p2,t2] = find_midpoints(p,t);

    [basis,Qh] = solve(p,p2,e,t,t2,u,grad_u_r,grad_u_z,n);
    [err(1),grad_err(1),max_err(1)] = errors_exact_weighted_p2(p,t,p2,t2,basis,Qh,n);

    for i = 2:mesh_level
        % Refine mesh to next level
        [p,e,t]=refinemesh(g,p,e,t,it,'regular');
        %pdemesh(p,e,t, 'NodeLabels','on', 'ElementLabels','on');

        % Find the midpoints for P2 nodal points
        [p2,t2] = find_midpoints(p,t);

        [basis,Qh] = solve(p,p2,e,t,t2,u,grad_u_r,grad_u_z,n);
        [err(i),grad_err(i),max_err(i)] = errors_exact_weighted_p2(p,t,p2,t2,basis,Qh,n);

    end
    display_errors(err,grad_err,max_err)

end
% mesh level must be greater than 1

% end main
end

% subfunction
function [basis,Qh] = solve(p,p2,e,t,t2,u,grad_u_r,grad_u_z,n)
    basis = basis_functions_weighted_p2(p,t,p2,t2);
    S = stiffness_matrix_weighted_p2(p,t,p2,t2,basis,n);
    b = create_b_p2(p,t,p2,t2,basis,u,grad_u_r,grad_u_z,n);
    Qh = S\b;

    figure();
    pdeplot([p,p2],e,t, 'XYData',Qh, 'ZData', Qh, 'Mesh', 'on');
end
