% PnP-Steric Test main fucntion
% Update on 2015/09/13
% Change log:
%        output result every savesteps 
%        change the Bc input type
%
% Change log:
%        Add loop for N
%        Add loop for StericG
%        Plot error of left boundary of ion 1.


function TestMain
clear all 
close all

global kB Temp e_unit dielec savesteps 
global DiffIon ValIon StericG HOT phi
global AlphaBetaIon AreaStr dt
global IonBc PhiBc dCdxBc Euler
global xmin xmax tstart tend
global N CFL c_tau AlphaBetaPhi C_ex x w

% Input parameter
[p]=parameters;
kB= p.kB; Temp = p.Temp; e_unit = p.e_unit; dielec = p.dielec;
ValIon = p.ValIon; StericG = p.StericG(:,:,:); HOT = p.HOT(:,:,:);
DiffIon = p.DiffIon; AreaStr = p.AreaStr; %epslion = p.epslion; 
IonBc = p.IonBc; PhiBc = p.PhiBc;
dCdxBc = p.dCdxBc;  Euler=p.Euler; 
xmin = p.xmin; xmax = p.xmax; dt=p.dt;
tstart = p.tstart; tend = p.tend; 
AlphaBetaIon = p.AlphaBetaIon; AlphaBetaPhi = p.AlphaBetaPhi;
CFL = p.CFL; c_tau = p.c_tau; NN = p.NN; savesteps = p.savesteps;

runtime = zeros(length(NN));
errL_inf = zeros(length(NN));
errL_2 = zeros(length(NN));
convL_inf = zeros(length(NN));
convL_2 = zeros(length(NN));
for ii=1:length(NN)
    N = NN(ii);
    [Conc] = ModelPnP_PN_ODE15s;
    Cex = C_ex(x);
    err_i = abs(Cex - Conc);
    errL_inf(ii) = max(abs(err_i));
    errL_2(ii) = sqrt(sum(err_i.^2 .* w));
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
    fprintf('%4d %s %12.4e %s %6.2f %2s %12.4e %2s \n',NN(ii), '&', ...
            errL_inf(ii), '&', convL_inf(ii), '&', ...
            errL_2(ii), '\\')
end



