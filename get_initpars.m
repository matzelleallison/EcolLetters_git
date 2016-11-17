function pars = get_initpars(X_0, par)
% created by Alli Matzelle 2016/5/7
  
  %% Syntax
  % pars = get_initpars(X_0, par)
  
  %% Description
  % sets life stage parameters (independently of other parameters so initial food
  % density accounted for) and inital conditions of state variables
  %
  % Input
  %
  % * par: structure with values of parameters
  % * X_0: initial food density to calculate quantities at birth
  %  
  % Output
  %
  % * lifePars : structure with values of parameters
  
%% create parameter structure...
vars_pull(par)  % unpack par structure
cPar = parscomp_st(par);    vars_pull(cPar);

f0 = X_0 / (K_hs + X_0);    % initial functional response
pars_U0 = [V_Hb, g, k_J, k_M, v];
[U_0, L_b, info] = initial_scaled_reserve(f0, pars_U0);

% birth
pars. L_b = L_b;  % length at birth at f

pars.E_0 = U_0 * p_Am;   % reserve
pars.V_0 = L_b^3;        % structural volume
pars.E_H0 = E_Hb;        % maturity
pars.E_R0 = 0;           % reproduction buffer

% metamorphosis
pars_lj = [g k l_T v_Hb v_Hj v_Hp];  % compose parameter vector 
[t_j, t_p, t_b, l_j, l_p, l_b, l_i, rho_j, rho_B, info] = get_tj(pars_lj, f0);% -, scaled times & lengths at f
pars.L_j = l_j * L_m;          % cm, structural length at metamorphosis

% pack initial conditions
% pars.EVHR_init = [E_0; V_0; E_H0; E_R0];  