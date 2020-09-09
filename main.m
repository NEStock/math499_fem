function [err,grad_err,max_err] = main(u,gd,sf,ns,mesh_level,opt)
%MAIN Summary of this function goes here
%   Detailed explanation goes here
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

    [basis,Qh] = solve(p,p2,e,t,t2,f);
    [err(1),grad_err(1),max_err(1)] = errors_exact_p2(p,t,p2,t2,basis,Qh,u,grad_u_x,grad_u_y);
    % errors_exact_weighted_p2(...)

    for i = 2:mesh_level
        [p,e,t]=refinemesh(g,p,e,t,it,'regular');
        %pdemesh(p,e,t, 'NodeLabels','on', 'ElementLabels','on');

        [p2,t2] = find_midpoints(p,t);

        [basis,Qh] = solve(p,p2,e,t,t2,f);
        [err(i),grad_err(i),max_err(i)] = errors_exact_p2(p,t,p2,t2,basis,Qh,u,grad_u_x,grad_u_y);

    end
    display_errors(err,grad_err,max_err)

else
    [p2,t2] = find_midpoints(p,t);
    
    [basis,Qh] = solve(p,p2,e,t,t2,f);
    [err(1),grad_err(1),max_err(1)] = errors_exact_p2(p,t,p2,t2,basis,Qh,u,grad_u_x,grad_u_y);
end

% end main
end

% subfunction
function [basis,Qh] = solve(p,p2,e,t,t2,u)
    basis = basis_functions_weighted_p2(p,t,p2,t2);
    S = stiffness_matrix_weighted_p2(p,t,p2,t2,basis);
    b = create_b_p2(p,t,p2,t2,basis,u);
    Qh = S\b;

    figure();
    pdeplot([p,p2],e,t, 'XYData',Qh, 'ZData', Qh, 'Mesh', 'on');
end