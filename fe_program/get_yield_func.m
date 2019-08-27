function [yield_func] = get_yield_func(C, inv_Fp, state_var, ...
                                       params, delta_gamma)
%% GET_YIELD_FUNC returns value of the yield function.
% Input
% C - right Cauchy-Green deformation tensor (in 9-vector representation).
% state_var - structure with state variables (plastic part of 
% deformation gradient, equivalent plastic strain).
% inv_Fp - inverse of plastic part of deformation gradient.
% params - structure containing material parameters.
% delta_gamma - increment of equivalent plastic strain (scalar).
%
% Output:
% yield_func - scalar value of the yield function.

yield_stress = params.init_yield;
H = params.harden;
ident=[1 1 1 0 0 0 0 0 0]';

% Mandel stress:
[M,~] = get_mandel_pk2interm_stress(C, inv_Fp, params);

M_dev = M - (M'*ident)*ident;
norm_M_dev = sqrt( m_2_v9( v9_2_m(M_dev)' )' * M_dev );

yield_func = sqrt(1.5) * norm_M_dev - ...
             (yield_stress + H * (state_var.eq_plast_strain + delta_gamma));
end