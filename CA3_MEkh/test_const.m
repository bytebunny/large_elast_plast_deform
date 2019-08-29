clear all
close all
format compact
format short e
%%%

addpath( genpath('../fe_program/') ) % Add path to my routines. /Rostyslav

%definition of time history
tmax=100; no_of_timesteps=1.e2;
t=linspace(0,tmax,no_of_timesteps);

%max value of engineering strain
eps_max=0.2;


%definition of material parameters
Emod=2e5; % [MPa]
v=0.3;
mu=Emod/(2*(1+v)); lambda=Emod*v/((1-2*v)*(1+v));
para.mu = mu; para.lambda = lambda;
para.harden = Emod/20; para.init_yield = 400; % [MPa]
para.max_iter = 20; % Max number of Newton iterations.
para.error_tol = 1.e-3; % Error tolerance for Newton's iteration.

%small strain stiffness tensor (linear isotropic Hooke elasticity)
ident=[1 1 1 0 0 0 0 0 0]';
C_hooke=2*mu* f9_open_u_9(ident,ident)+lambda*ident*ident';

% allocate space and initialise internal variables
state_var = struct('eq_plast_strain',{0}, 'Fp',{ident});

for i=1:no_of_timesteps
   state_var_old = state_var;
   %deformation gradient  as a 3x3 matrix
   F12=eps_max*(i-1)/no_of_timesteps;
   Fm=... %[1 F12 0; 0 F11^(-v) 0; 0 0 F11^(-v)];
   [1 eps_max*(i-1)/no_of_timesteps 0; 0 1 0; 0 0 1];
   %translation to 9 vector representation
   F=m_2_v9(Fm);
   %right Cauchy Green deformation tensor as a 3x matrix
   Cm=Fm'*Fm;
   %translation to 9 vector representation
   C=m_2_v9(Cm);
   
   %call constitutive model
   %[S2,dS2_dE]=neo_hooke(C,para(1:2));
   [S2,state_var]=neo_hooke_plast(C,state_var_old,para);
   %numerical computation of stifness tensor to check that the
   %implementation
   for jj=1:9
      Cdiff=C; num_pert=1.e-7; Cdiff(jj)=Cdiff(jj)+num_pert;
      [S2diff,~]=neo_hooke_plast(Cdiff,state_var_old,para);
      dS2_dE(:,jj)=2*(S2diff-S2)/num_pert;
   end
   %sym dS2_dE, otherwise otherwise FE problem will not converge ... (sym
   %assumption is assumed in the derivations of the element stiffness)
   %
   %rows
   dS2_dE(4,:)=0.5*( dS2_dE(4,:)+dS2_dE(8,:) );
   dS2_dE(8,:)=dS2_dE(4,:);
   dS2_dE(5,:)=0.5*( dS2_dE(5,:)+dS2_dE(9,:) );
   dS2_dE(9,:)=dS2_dE(5,:);
   dS2_dE(6,:)=0.5*( dS2_dE(6,:)+dS2_dE(7,:) );
   dS2_dE(7,:)=dS2_dE(6,:);
   %columns
   dS2_dE(:,4)=0.5*( dS2_dE(:,4)+dS2_dE(:,8) );
   dS2_dE(:,8)=dS2_dE(:,4);
   dS2_dE(:,5)=0.5*( dS2_dE(:,5)+dS2_dE(:,9) );
   dS2_dE(:,9)=dS2_dE(:,5);
   dS2_dE(:,6)=0.5*( dS2_dE(:,6)+dS2_dE(:,7) );
   dS2_dE(:,7)=dS2_dE(:,6);

   
   %compute Cauchy stress (as a 3x3 matrix) as a push forward of S2
   sigma=1/det(Fm)*Fm*v9_2_m(S2)*Fm';
   
   %compute the sitffness c by push forward operation
   Ft=m_2_v9(Fm');
   c=1/det(Fm)*f9_open_u_9(F,F)*dS2_dE*f9_open_u_9(Ft,Ft);
   
   
   %save plot data
   sigma_out(i)=sigma(1,2); strain_out(i)=Fm(1,2); 
end

figure(1)
hold on
plot(strain_out,sigma_out,'b-')

f_name = ['../doc/data/sigma12_eps12.dat'];
f_id = fopen(f_name,'w');
header = '# eng strain 12 [-]   Cauchy stress 12, [MPa]';
fprintf(f_id, '%s\n', header);
fprintf(f_id, '%.4f             %.4f\n', [strain_out; sigma_out]);
fclose(f_id);






