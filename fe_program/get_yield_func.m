function [yield_func, yield_jac, inv_Fp] = get_yield_func(C, inv_Fp_old, state_var, ...
                                                  params, delta_gamma)
%% GET_YIELD_FUNC returns value of the yield function.
% Input
% C - right Cauchy-Green deformation tensor (in 9-vector representation).
% inv_Fp_old - inverse of plastic part of deformation gradient from 
% previous time step.
% state_var - structure with state variables (plastic part of 
% deformation gradient, equivalent plastic strain).
% params - structure containing material parameters.
% delta_gamma - increment of equivalent plastic strain (scalar).
%
% Output:
% yield_func - scalar value of the yield function.
% yield_jac - Jacobian of yield function.
% inv_Fp - updated inverse of plastic part of deformation gradient.

yield_stress = params.init_yield;
mu = params.mu; lambda = params.lambda; H = params.harden;
ident=[1 1 1 0 0 0 0 0 0]';

inv_Fp_m_old = v9_2_m(inv_Fp_old);

[M_old,~] = get_mandel_pk2interm_stress(C, inv_Fp_old, params);
M_dev_old = M_old - 1/3*(M_old'*ident)*ident;
norm_M_dev_old = sqrt( M_dev_old' * M_dev_old );
dyield_dM_old = sqrt(1.5) * M_dev_old / (norm_M_dev_old + 1.e-16); % Avoid division by 0.
dyield_dM_m_old = v9_2_m(dyield_dM_old);

% Update inverse of plastic part of deformation gradient for next time
% step:
inv_Fp_m = inv_Fp_m_old - delta_gamma * inv_Fp_m_old * dyield_dM_m_old;
inv_Fp = m_2_v9(inv_Fp_m);

% Mandel stress:
[M,PK2_interm] = get_mandel_pk2interm_stress(C, inv_Fp, params);

M_dev = M - 1/3*(M'*ident)*ident;

norm_M_dev = sqrt( M_dev' * M_dev );

yield_func = sqrt(1.5) * norm_M_dev - ...
             (yield_stress + H * (state_var.eq_plast_strain + delta_gamma));
         
%% Jacobian for local problem:
dyield_dM = sqrt(1.5) * M_dev / (norm_M_dev + 1.e-16); % Avoid division by 0.
C_m = v9_2_m(C);
Ce_m = inv_Fp_m' * C_m * inv_Fp_m;
Ce = m_2_v9(Ce_m);
Ce_inv = m_2_v9( inv(Ce_m) );
Je=sqrt(det(Ce_m));
dPK2_interm_dCe = 0.5*lambda * (Ce_inv * Ce_inv') + ...
                  (mu - lambda*log(Je)) * f9_open_u_9(Ce_inv,Ce_inv);
        
dM_dCe = f9_open_u_9(ident,PK2_interm) + ...
         f9_open_u_9(Ce,ident) * dPK2_interm_dCe;
             
dCe_dinv_Fp = f9_open_l_9(ident, m_2_v9(inv_Fp_m'*C_m)) + ...
              f9_open_u_9(m_2_v9(inv_Fp_m'*C_m),ident);
                  
dinv_Fp_ddelta_gamma = - m_2_v9(inv_Fp_m_old * dyield_dM_m_old);
        
yield_jac = dyield_dM' * dM_dCe * dCe_dinv_Fp * dinv_Fp_ddelta_gamma - H;
end