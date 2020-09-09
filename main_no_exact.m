function [err,grad_err,max_err] = main_no_exact(u,gd,sf,ns,mesh_level,n)
%MAIN Summary of this function goes here
%   Detailed explanation goes here
%
% Syntax:
%     [err,grad_err,max_err] = 
%       main_no_exact(u,gd,sf,ns,mesh_level, n)
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
%    pdepoly([0,1,1,0], [-1,-1,1,1]);
%    [err,grad_err,max_err] = 
%       main_no_exact(u,gd,sf,ns,mesh_level, n)
% Dependencies:
%    basis_functions_weighted_p2.m
%    display_errors.m
%    errors_no_exact_weighted_p2.m
%    prolongation_matrix.m
%    stiffness_matrix_weighted_p2.m
%
% Note:
%    prolongation_matrix.m not working yet. Thus errors are wrong.
%
% Author: Nicole Stock
% Date: Fall 2020

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

    [basis,Qh] = solve(p,p2,e,t,t2,u,n);
    % errors_no_exact_weighted_p2(p,t,p2,t2,basis,u_h_km1,u_h_k,k)

    for i = 2:mesh_level
        basis_1 = basis;
        Qh_1 = Qh;
        t_1 = t;
        t2_1 = t2;
        [~,nodes_1] = size(p);
        [~,mid_nodes_1] = size(p2);
        
        % Refine mesh to next level
        [p,e,t]=refinemesh(g,p,e,t,it,'regular');

        % Find the midpoints for P2 nodal points
        [p2,t2] = find_midpoints(p,t);

        % Create prolongation matrix to extend mesh level - 1 vectors
        P = prolongation_matrix(nodes_1,mid_nodes_1,p,p2,t_1,t2_1,basis_1,n);

        % Extend Qh
        Qh_extended = P*Qh_1;
        
        % Solve
        [basis,Qh] = solve(p,p2,e,t,t2,u,n);
        [err(i),grad_err(i),max_err(i)] = errors_no_exact_weighted_p2(p,t,p2,t2,basis,Qh_extended,Qh,n);
    end
    
    display_errors(err,grad_err,max_err)
end
% mesh level must be greater than 1

% end main
end

% subfunction
function [basis,Qh] = solve(p,p2,e,t,t2,u,n)
    basis = basis_functions_weighted_p2(p,t,p2,t2);
    S = stiffness_matrix_weighted_p2(p,t,p2,t2,basis,n);
    b = create_b_p2(p,t,p2,t2,basis,u);
    Qh = S\b;

    figure();
    pdeplot([p,p2],e,t, 'XYData',Qh, 'ZData', Qh, 'Mesh', 'on');
end
