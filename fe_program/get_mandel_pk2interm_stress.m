function [M, PK2_interm] = get_mandel_pk2interm_stress(C, inv_Fp, params)
%% GET_MANDEL_PK2INTERM_STRESS returns Mandel and intermediate 2nd Piola-Kirchhoff stresses.
% Input:
% C - right Cauchy-Green deformation tensor (in 9-vector representation).
% inv_Fp - inverse of the plastic part of deformation gradient.
% params - structure containing material parameters.
%
% Output:
% M - Mandel stress (in 9-vector representation).
% PK2_interm - intermediate 2nd Piola-Kirchhoff stress (in 9-vector 
% representation).

mu = params.mu; lambda = params.lambda;
ident=[1 1 1 0 0 0 0 0 0]';

inv_Fp_m = v9_2_m(inv_Fp);
% Intermediate elastic deformation tensor:
Ce_m = inv_Fp_m' * v9_2_m(C) * inv_Fp_m;

inv_Ce = m_2_v9( inv(Ce_m) );

det_Ce=det(Ce_m);
J=sqrt(det_Ce);
% Intermediate 2nd Piola-Kirchhoff stress:
PK2_interm = mu * ident + (lambda * log(J) - mu) * inv_Ce;
% Mandel stress:
M = m_2_v9( Ce_m * v9_2_m(PK2_interm) );
end