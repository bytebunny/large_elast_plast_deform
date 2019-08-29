function [PK2_stress, state_var] = ...
    neo_hooke_plast(C, state_var_old, params)
%% NEO_HOOKE_PLAST Returns stress and new state variables.
% The function serves as a material routine.
% NOTE: material stiffness is expected to be computed numerically by the
% element routine!
% Inputs:
% C - right Cauchy-Green deformation tensor (in 9-vector representation).
% state_var_old - structure with state variables (plastic part of 
% deformation gradient, equivalent plastic strain) from previous time step.
% params - structure containing material parameters.
%
% Output:
% PK2_stress - 2nd Piola-Kirchhoff stress (in 9-vector representation).
% state_var - updated input structure.

max_iter = params.max_iter; % Max number of Newton iterations for local constitutive problem.
error_tol = params.error_tol; % Error tolerance for Newton's iteration.

delta_gamma = 0; % Initialise increment of equivalent plastic strain.
% Get plastic part of deformation gradient from previous time step:
Fp = state_var_old.Fp; 
inv_Fp_old = m_2_v9( inv( v9_2_m(Fp) ) );

[yield_func,~,~] = get_yield_func(C, inv_Fp_old, state_var_old, ...
                                        params, delta_gamma);

if yield_func <= 0
    %% Elastic response
    state_var = state_var_old;
    inv_Fp = inv_Fp_old;
    fprintf('/// Elastic step in constitutive driver.\n')
else
    %% Plastic response:     
    for iter=1:max_iter % Use Newton's method to solve local constitutive problem:
        [yield_func, yield_jac, inv_Fp] = get_yield_func(C, inv_Fp_old, state_var_old, ...
                                        params, delta_gamma);

        delta_gamma = delta_gamma - yield_func / yield_jac;

        if abs(yield_func) <= error_tol
           fprintf('/// Local constitutive problem converged in %d iterations.\n',iter)
           break
        end
    end
    if iter == max_iter
        error('/// Local constitutive problem did NOT converge in %d iterations.\n',max_iter)
    end
    % Update state_variables:
    state_var.Fp = m_2_v9( inv( v9_2_m(inv_Fp) ) );
    state_var.eq_plast_strain = state_var_old.eq_plast_strain + delta_gamma;
end
% Intermediate 2nd Piola-Kirchhoff stress:
[~, PK2_interm] = get_mandel_pk2interm_stress(C, inv_Fp, params);
% 2nd Piola-Kirchhoff:
inv_Fp_m = v9_2_m(inv_Fp);
PK2_stress = m_2_v9( inv_Fp_m * v9_2_m(PK2_interm) * inv_Fp_m' );
end