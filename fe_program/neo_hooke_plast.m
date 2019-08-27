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
ident=[1 1 1 0 0 0 0 0 0]';
mu = params.mu; lambda = params.lambda; H = params.harden;

delta_gamma = 0; % Initialise increment of equivalent plastic strain.
C_m = v9_2_m(C);
% Get plastic part of deformation gradient from previous time step:
Fp = state_var_old.Fp; 
inv_Fp = m_2_v9( inv( v9_2_m(Fp) ) );

if get_yield_func(C, inv_Fp, state_var_old, params, 0) <= 0
    %% Elastic response
    state_var = state_var_old;
else
    %% Plastic response: 
    [M,~] = get_mandel_pk2interm_stress(C, inv_Fp, params);
    M_dev = M - (M'*ident)*ident;
    norm_M_dev = sqrt( m_2_v9( v9_2_m(M_dev)' )' * M_dev );
    dyield_dM = sqrt(1.5) * M_dev / (norm_M_dev + 1.e-16); % Avoid division by 0.
    
    for iter=1:max_iter % Use Newton's method to solve local constitutive problem:    
        inv_Fp_m_old = v9_2_m(inv_Fp); % Store for jacobian.
        dyield_dM_m_old = v9_2_m(dyield_dM);
             
        % Update inverse of plastic part of deformation gradient:
        inv_Fp = inv_Fp - delta_gamma * v9_2_m(inv_Fp) * v9_2_m(dyield_dM);
        inv_Fp_m = v9_2_m(inv_Fp);
        
        %% Jacobian for local problem:
        [M,PK2_interm] = get_mandel_pk2interm_stress(C, inv_Fp, params);
        M_dev = M - (M'*ident)*ident;
        norm_M_dev = sqrt( m_2_v9( v9_2_m(M_dev)' )' * M_dev );
        dyield_dM = sqrt(1.5) * M_dev / (norm_M_dev + 1.e-16); % Avoid division by 0.
        
        Ce_m = inv_Fp_m' * C_m * inv_Fp_m;
        Ce_inv = m_2_v9( inv(Ce_m) );
        Je=sqrt(det(Ce_m));
        dPK2_interm_dCe = 0.5*lambda * (Ce_inv * Ce_inv') + ...
                          (mu - lambda*log(Je)) * f9_open_u_9(Ce_inv,Ce_inv);
        
        dM_dCe = f9_open_u_9(ident,PK2_interm) + ...
                 f9_open_u_9(Ce,ident) * dPK2_interm_dCe;
             
        dCe_dinv_Fp = f9_open_l_9(ident, m_2_v9(inv_Fp_m'*C_m)) + ...
                      f9_open_l_9(m_2_v9(inv_Fp_m'*C_m),ident);
                  
        dinv_Fp_ddelta_gamma = - m_2_v9(inv_Fp_m_old * dyield_dM_m_old);
        
        yield_jac = dyield_dM * dM_dCe * dCe_dinv_Fp * dinv_Fp_ddelta_gamma ...
                    - H;
        %% Update guess of increment of equivalent plastic strain
        yield_func = get_yield_func(C, inv_Fp, state_var_old, params, delta_gamma);
        delta_gamma = delta_gamma - yield_jac * yield_func;
        
        if get_yield_func(C, inv_Fp, state_var_old, params, delta_gamma) <= error_tol
           fprintf('/// Local constitutive problem converged in %d iterations.\n',iter)
           break
        end
    end
    if iter == max_iter
        error('/// Local constitutive problem did NOT converge in %d iterations.\n',max_iter)
    end
    % Update plastic part of deformation gradient after convergence:
    inv_Fp = inv_Fp - delta_gamma * v9_2_m(inv_Fp) * v9_2_m(dyield_dM);
    % Update state_variables:
    state_var.Fp = m_2_v9( inv( v9_2_m(inv_Fp) ) );
    state_var.eq_plast_strain = state_var.eq_plast_strain + delta_gamma;
end
% Intermediate 2nd Piola-Kirchhoff stress:
[~, PK2_interm] = get_mandel_pk2interm_stress(C, inv_Fp, params);
% 2nd Piola-Kirchhoff:
inv_Fp_m = v9_2_m(inv_Fp);
PK2_stress = m_2_v9( inv_Fp_m * v9_2_m(PK2_interm) * inv_Fp_m' );
end