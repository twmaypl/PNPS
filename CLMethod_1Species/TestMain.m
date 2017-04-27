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

global kB Temp e_unit dielec epslion savesteps
global DiffIon ValIon  StericG HOT T
global AlphaBetaIon AreaStr 
global IonBc PhiBc dCdxBc 
global xmin xmax tstart tend Euler dt 
global N CFL c_tau AlphaBetaPhi w %phi C_ex LGL_x w

% input parameters
[p]=parameters;
kB= p.kB; Temp = p.Temp; e_unit = p.e_unit; dielec = p.dielec;
ValIon = p.ValIon; StericG = p.StericG; HOT = p.HOT;
DiffIon = p.DiffIon; AreaStr = p.AreaStr; epslion = p.epslion; 
IonBc = p.IonBc; PhiBc = p.PhiBc;% ep=p.ep;
dCdxBc = p.dCdxBc; xmin = p.xmin; xmax = p.xmax;
tstart = p.tstart; tend = p.tend; Euler = p.Euler; dt=p.dt;
AlphaBetaIon = p.AlphaBetaIon; AlphaBetaPhi = p.AlphaBetaPhi;
CFL = p.CFL; c_tau = p.c_tau; NN = p.NN; savesteps = p.savesteps;

for ii=1:length(NN)
N = NN(ii);
[Conc]=ModelPnP_PN_ODE15s;
end

