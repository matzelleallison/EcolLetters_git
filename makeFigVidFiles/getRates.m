function [Lw, Ww, Wd, pX, pM, pLw, pWw, Tb_out] = getRates(DEBdata, smpl_interval)
% calculates Lw, Ww, Wd, pX, pM, pLw, pWw 
% called by: compare_noiseNonoise.m

% load parameters
load('metaData.mat')
load('data.mat')
estimData.metaData = metaData;  estimData.data = data;
% set par structure
[par, txtPar] = set_pars(smpl_interval, estimData);
cPar = parscomp_st(par);
vars_pull(cPar);

% interpolate before doing calculations
nt = 1 : (365.25*10);
f = interp1(DEBdata.t, real(DEBdata.f), nt);
Tb_out = interp1(DEBdata.t, real(DEBdata.Tb_out), nt);
X_out = interp1(DEBdata.t, real(DEBdata.X_out), nt);
Lw = interp1(DEBdata.t, real(DEBdata.Lw), nt);
Ww = interp1(DEBdata.t, real(DEBdata.W_w), nt);
Wd = interp1(DEBdata.t, real(DEBdata.W_d), nt);
V = interp1(DEBdata.t, real(DEBdata.V), nt);

% temperature and shape parameters
T_1 = par.T1;
Tpars = [par.T_A par.T_L par.T_H par.T_AL par.T_AH];
c_T = tempcorr(Tb_out + 273.15, T_1, Tpars);
initpars = get_initpars(X_out(1), par);
M = shapecorr(Lw, initpars.L_b, initpars.L_j); % shape correction factor

% correct for temperature and shape
p_AmT = c_T .* cPar.p_Am .* M; % specific maximum assimilation rate
p_XmT = p_AmT ./ par.kap_X;  % max specific ingestion rate
p_MT = c_T .* par.p_M;   % somatic maintenance flux

% calculate rate(t)
pX = f .* p_XmT .* V.^(2/3);
pM = p_MT .* V;             % somatic maintenance flux
pLw = diff(Lw); % growth rate (cm/day)
pWw = diff(Ww); % growth rate (g/day)
end