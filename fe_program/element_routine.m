function [f,k,state] = element_routine(x, y, u, du, params, state_old, t)
%% ELEMENT_ROUTINE Returns element bar force and stiffness in 2D.
% Input:
% x - initial X coordinates of nodes.
% y - same but Y.
% u - element displacements.
% du - increment of element displacements.
% params - structure with material and cross-sectional parameters.
% state_old -- structure with element state variables from previous time step.
% t - element thickness, [mm].
%
% Output:
% f - element force vector.
% k - element stiffness matrix.
% state -- structure with updated state variables.

% Update node positions:
x_new = x + u(1:2:end) + du(1:2:end);
y_new = y + u(2:2:end) + du(2:2:end);
xy0 = [x', y'];
xy = [x_new', y_new'];

[f,state] = get_intern_eq_forces(xy0, xy, params, state_old, t);

[~,~,k] = get_tangent_matrices(xy0, xy, params, state_old, t);
end