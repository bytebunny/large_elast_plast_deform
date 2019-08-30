function [M_tr, M_tr_dev] = get_mandel_trial_stress(C, inv_Fp_old, params)
%% GET_MANDEL_TRIAL_STRESS returns Mandel trial stress and its deviatoric part.
% Input:
% C - right Cauchy-Green deformation tensor (in 9-vector representation).
% inv_Fp_old - inverse of the plastic part of deformation gradient from
% previous time step.
% params - structure containing material parameters.
%
% Output:
% M_tr - Mandel trial stress (in 9-vector representation).
% M_tr_dev - deviatoric part of M_tr.

mu = params.mu; lambda = params.lambda;
ident=[1 1 1 0 0 0 0 0 0]';

inv_Fp_old_m = v9_2_m(inv_Fp_old);
C_m = v9_2_m(C);
J=sqrt(det(C_m));

% Trial intermediate elastic deformation tensor:
Ce_tr = inv_Fp_old_m' * C_m * inv_Fp_old_m / J;

M_tr_dev = m_2_v9( mu * logm(Ce_tr) );

% Mandel trial stress:
M_tr = M_tr_dev + (lambda + 2/3 * mu) * log(J) * ident;
end