function [err] = weighted_HL_k_0_p1_e(f,grad_f_r,grad_f_z,gd,sf,ns,mesh,n,u,grad_u_r,grad_u_z)
%WEIGHTED_HL_K_0_P1_E Weighted Hodge Laplacian with k = 0. 
%   This program is set up to be given an exact solution.
%
% Syntax:
%     [err] = 
%       weighted_HL_k_0_p1_e(u,grad_u_r,grad_u_z,gd,sf,ns,mesh,n,u,grad_u_r,grad_u_z)
%
% Inputs:
%     f - given function
%     grad_f_r - gradient(f) w.r.t r
%     grad_f_z - gradient(f) w.r.t z
%     gd,sf,ns - outputs of pdepoly specifying domain
%     mesh - max mesh level
%     n - Hodge Laplacian on Axisymmetrix Domain and its discretization
%     weight
%     u - exact solution
%     grad_u_r - gradient(u) w.r.t r
%     grad_u_z - gradient(u) w.r.t z
%
% Outputs:
%     err - array of L2 errors for mesh levels corresponding to indices
%
% Usage Exampled:
%    addpath ../../data ../data/
%    [f,grad_f_r,grad_f_z,u,grad_u_r,grad_u_z] = get_data7();
%    mesh = 8;
%    n = 1;
%    pdepoly([0,1,1,0], [0,0,1,1]);
%       (OR) [gd,sf,ns] = get_gd_sf_ns([0,1,1,0],[0,0,1,1]);
%    [err] = weighted_HL_k_0_p1_e(f,grad_f_r,grad_f_z,gd,sf,ns,mesh,n,u,grad_u_r,grad_u_z);
% Dependencies:
%    basis_functions_weighted_HL_k_0_p1.m
%    create_b_HL_k_0_p1.m
%    display_errors.m
%    errors_exact_weighted_HL_k_0_p1.m
%    stiffness_matrix_weighted_HL_k_0_p1.m
%
% Author: Nicole Stock
% Date: Fall 2020

addpath('../../')

model=createpde(1);
g=decsg(gd,sf,ns);
geometryFromEdges(model,g);
[p,e,t]=initmesh(g,'hmax',inf);
%pdemesh(p,e,t, 'NodeLabels','on', 'ElementLabels','on');

% Init error vectors
err = zeros(1,mesh);

if mesh > 1
    % To ensure we refine every triangle the same
    [~,num_node]=size(p);
    it=zeros(1,num_node);
    for i=1:num_node
        it(i)=i;
    end

    [basis,Qh] = solve(p,e,t,f,grad_f_r,grad_f_z,n);
    [err(1)] = errors_exact_weighted_HL_k_0_p1(p,t,basis,Qh,n,u,grad_u_r,grad_u_z);

    for i = 2:mesh
        % Refine mesh to next level
        [p,e,t]=refinemesh(g,p,e,t,it,'regular');
        %pdemesh(p,e,t, 'NodeLabels','on', 'ElementLabels','on');

        [basis,Qh] = solve(p,e,t,f,grad_f_r,grad_f_z,n);
        [err(i)] = errors_exact_weighted_HL_k_0_p1(p,t,basis,Qh,n,u,grad_u_r,grad_u_z);

    end
    display_errors(err)

end
% mesh level must be greater than 1

% end main
end

% subfunction
function [basis,Qh] = solve(p,e,t,f,grad_f_r,grad_f_z,n)
    basis = basis_functions_weighted_HL_k_0_p1(p,t);
    S = stiffness_matrix_weighted_HL_k_0_p1(p,t,basis,n);
    b = create_b_HL_k_0_p1(p,t,basis,f,grad_f_r,grad_f_z,n);
    Qh = S\b;

   % figure();
   % pdeplot([p,p2],e,t, 'XYData',Qh, 'ZData', Qh, 'Mesh', 'on');
end
