function [PK2_stress, state_var] = ...
    neo_hooke_plast_log_strain(C, state_var_old, params)
%% NEO_HOOKE_PLAST_LOG_STRAIN Returns stress and new state variables.
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

mu = params.mu; lambda = params.lambda; H = params.harden;
ident=[1 1 1 0 0 0 0 0 0]';

% Get plastic part of deformation gradient from previous time step:
Fp = state_var_old.Fp; 
inv_Fp_old = m_2_v9( inv( v9_2_m(Fp) ) );

yield_func = get_yield_func_log_strain(C,inv_Fp_old,state_var_old,params);

if yield_func <= 0
    %% Elastic response    
    state_var = state_var_old;
    inv_Fp = inv_Fp_old;
    
    % Trial Mandel stress:
    [M, ~] = get_mandel_trial_stress(C, inv_Fp, params);
else
    %% Plastic response:     
    % No Newton iteration, the local problem has closed-form solution.
    % Increment of equivalent plastic strain:
    delta_gamma = yield_func / (3*mu + H);

    [~,M_tr_dev] = get_mandel_trial_stress(C,inv_Fp_old,params);
    % Tensor nu from lecture notes:
    dyield_dM = sqrt(1.5) * M_tr_dev / (sqrt( M_tr_dev' * M_tr_dev )+1.e-16);
    A = expm( delta_gamma * v9_2_m(dyield_dM) );
    
    % Update plastic part of deformation gradient:
    Fp_m = A * v9_2_m(Fp);
    inv_Fp = m_2_v9( inv(Fp_m) );
    
    % Update Mandel stress:
    J=sqrt(det( v9_2_m(C) ));
    M = M_tr_dev - 2 * mu * delta_gamma * dyield_dM + ...
        (lambda + 2/3 * mu) * log(J) * ident;
    
    % Update state_variables:
    state_var.Fp = m_2_v9( inv( v9_2_m(inv_Fp) ) );
    state_var.eq_plast_strain = state_var_old.eq_plast_strain + delta_gamma;
end
% 2nd Piola-Kirchhoff:
inv_Fp_m = v9_2_m(inv_Fp);
inv_C_m = inv( v9_2_m(C) );
PK2_stress = m_2_v9( inv_C_m * v9_2_m(Fp)' * v9_2_m(M) * inv_Fp_m' );
end