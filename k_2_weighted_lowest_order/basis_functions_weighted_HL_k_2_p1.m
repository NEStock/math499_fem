function [basis_nodes,basis_NP1_edges,basis_RT_edges,basis_triangles] = basis_functions_weighted_HL_k_2_p1(p,t,ed,t_ed)
% BASIS_FUNCTIONS_WEIGHTED_HL_K_2_P1 - Create a piecewise 
%   basis function for each node of a triangulation
%
% Syntax:
%     [basis_nodes,basis_NP1_edges,basis_RT_edges,basis_triangles] = 
%           basis_functions_weighted_modified_nedelec_and_p1(p,t,ed,t_ed)
% 
% Inputs:
%     p - a 2xNumNodes matrix representing nodal coordinates.
%     t - a 4xNumTriangles matrix representing the element connectivity in 
%         terms of node IDs. The end row of T represents the geometry face ID 
%         to which the element belongs.
%     ed - a 2xNumEdges matrix representing each edge as a row with
%         starting node in column 1 and the ending node in column 2.
%     t_ed - a 3xNumTriangles matrix representing the which edges
%         correspond to which triangles. t_ed(i,T) represents the ith edge
%         in triangle T.
%
% Outputs:
%     basis_nodes - a matrix representing piece-wise basis functions for
%         each node in each triangle. basis(i,:,T) represents the
%         pieceiwise basis function for the ith node in triangle T.
%     basis_NP1_edges - a matrix representing piece-wise basis functions
%         for each edge in each triangle for the weighted fourier Nedelec
%         and P1 space. basis(i,:,T) represents the pieceiwise basis
%         function for the ith edge in triangle T.
%     basis_RT_edges - a matrix representing piece-wise basis functions
%         for each edge in each triangle for the weighted fourier Raviart 
%         Thomas space. basis(i,:,T) represents the pieceiwise basis 
%         function for the ith edge in triangle T.
%     basis_triangles - a vector representing piece-wise basis functions
%         for edge triangle. basis(1,T) represents the piecewise basis
%         function for the Tth triangle.
%
% Dependencies:
%     basis_functions_weighted_HL_p1.m
%     basis_functions_weighted_nedelec.m
%
% Author: Nicole Stock
% Date: Fall 2020

addpath('../raviart_thomas_weighted_fourier_lowest_order/');
addpath('../modified_nedelec_and_p1_weighted/');

  
% basis_RT_edges(:,i,T) = basis for the ith node of triangle T 
%       (fourier raviart thomas & P1)
% basis_triangles(:,T) = basis for the Tth triangle
[basis_RT_edges, basis_triangles] = basis_functions_weighted_fourier_rt(p,t,ed,t_ed);

% basis_nodes(:,i,T) = basis for the ith node of triangle T
% basis_NP1_edges(:,i,T) = basis for the ith node of triangle T 
%       (fourier nedelec & P1)
[basis_nodes,basis_NP1_edges] = basis_functions_weighted_modified_nedelec_and_p1(p,t,ed,t_ed);