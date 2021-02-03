function [basis_nodes,basis_edges] = basis_functions_weighted_modified_nedelec_and_p1(p,t,ed,t_ed)
% BASIS_FUNCTIONS_WEIGHTED_MODIFIED_NEDELEC_AND_P1 - Create a piecewise 
%   basis function for each node of a triangulation
%
% Syntax:
%     [basis_nodes,basis_edges] = basis_functions_weighted_modified_nedelec_and_p1(p,t,ed,t_ed)
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
%     basis_edges - a matrix representing piece-wise basis functions for 
%         each edge in each triangle. basis(i,:,T) represents the 
%         pieceiwise basis function for the ith edge in triangle T.
%
% Dependencies:
%     basis_functions_weighted_HL_p1.m
%     basis_functions_weighted_nedelec.m
%
% Author: Nicole Stock
% Date: Fall 2020
addpath('../k_0_weighted_p1/')
addpath('../nedelec_weighted_lowest_order/')
  
% basis_nodes(:,i,T) = basis for the ith node of triangle T
basis_nodes = basis_functions_weighted_HL_k_0_p1(p,t);

% basis_edges(:,i,T) = basis for the ith node of triangle T
basis_edges = basis_functions_weighted_nedelec(p,ed,t_ed);