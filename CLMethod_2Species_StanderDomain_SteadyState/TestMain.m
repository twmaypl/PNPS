% PnP-Steric Test main fucntion
% Update on 2016/01/11
% Change log:
%        Considered the standard domain [-1,1]
%        Using Chebyshev-Legendre Method
%        This is the steady-state test case


function TestMain
clear all 
close all

global kB Temp e_unit epslion savesteps
global DiffIon ValIon ep StericG HOT
global AlphaBetaIon AlphaBetaPhi
global IonBc PhiBc dCdxBc 
global xmin xmax tstart tend Euler
global N CFL c_tau dt phi x C_ex w LGL_x T

[p]=parameters;
kB= p.kB; Temp = p.Temp; e_unit = p.e_unit; 
ValIon = p.ValIon; StericG = p.StericG; HOT = p.HOT;
DiffIon = p.DiffIon; epslion = p.epslion; 
IonBc = p.IonBc; PhiBc = p.PhiBc; dt=p.dt; ep=p.ep;
dCdxBc = p.dCdxBc; xmin = p.xmin; xmax = p.xmax; Euler=p.Euler;
tstart = p.tstart; tend = p.tend; 
AlphaBetaIon = p.AlphaBetaIon; AlphaBetaPhi = p.AlphaBetaPhi;
CFL = p.CFL; c_tau = p.c_tau; NN = p.NN; savesteps = p.savesteps;

errL_inf = zeros(length(NN));
errL_2 = zeros(length(NN));
convL_inf = zeros(length(NN));
convL_2 = zeros(length(NN));
for ii=1:length(NN)
    N=NN(ii);
 [Conc] = ModelPnP_PN_ODE15s;

 Cex_lgl = [C_ex(LGL_x) C_ex(LGL_x)];
 Conc_lgl = T * Conc;
 err_i = abs(Cex_lgl - Conc_lgl);
 errL_inf(ii) = max(max(abs(err_i)));
 disp([size(err_i) size(repmat(w,1,2))])
 errL_2(ii) = sqrt(sum(sum(err_i.^2 .* repmat(w,1,2))));
end

% output error and convergence order
for ii =2:length(NN)
convL_inf(ii)=log(errL_inf(ii)/errL_inf(ii-1)) / log(NN(ii-1)/NN(ii));
convL_2(ii)=log(errL_2(ii)/errL_2(ii-1)) / log(NN(ii-1)/NN(ii));
end

fprintf('\n')
fprintf('%4s %1s %12s %1s %5s %2s %12s %s %5s \n', 'N', '&', 'L_inf err', '&', 'Order', ...
                                                        '&', 'L_2 err', '&', 'Order')
for ii=1:length(NN)
    fprintf('%4d %s %12.4e %s %6.2f %2s %12.4e %s %6.2f %2s \n',NN(ii), '&', ...
            errL_inf(ii), '&', convL_inf(ii), '&', ...
            errL_2(ii), '&', convL_2(ii),'\\')
end
 


