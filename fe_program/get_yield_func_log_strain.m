function [yield_func] = get_yield_func_log_strain(C, inv_Fp_old, state_var, ...
                                                  params)
%% GET_YIELD_FUNC_LOG_STRAIN returns value of the yield function.
% NOTE: this corresponds to the trial yield function from ch.5 of Magnus's 
% lecture notes on large elasto-plastic deformations.
% Input
% C - right Cauchy-Green deformation tensor (in 9-vector representation).
% inv_Fp_old - inverse of plastic part of deformation gradient from 
% previous time step.
% state_var - structure with state variables (plastic part of 
% deformation gradient, equivalent plastic strain).
% params - structure containing material parameters.
%
% Output:
% yield_func - scalar value of the yield function.

yield_stress = params.init_yield;
H = params.harden;

[~,M_tr_dev] = get_mandel_trial_stress(C, inv_Fp_old, params);

yield_func = sqrt(1.5) * sqrt( M_tr_dev' * M_tr_dev ) - ...
             (yield_stress + H * state_var.eq_plast_strain);
end