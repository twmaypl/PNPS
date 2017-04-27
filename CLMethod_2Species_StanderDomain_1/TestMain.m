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

global kB Temp e_unit epslion savesteps
global DiffIon ValIon ep StericG HOT
global AlphaBetaIon AlphaBetaPhi
global IonBc PhiBc dCdxBc 
global xmin xmax tstart tend Euler
global N CFL c_tau dt phi x

% % % % for id=2:2%1:11;
% % % % %filename = ['../Parameter/parameter',num2str(id),'.txt'];
% % % % filename = ['../Parameter/parameter',num2str(id),'.txt']
% % % % FID = fopen(filename,'r');
% % % % A = fscanf(FID, '%f');
% % % % fclose(FID);
% % % % td = A(1); kB = A(2); Temp = A(3); e_unit = A(4); dielec = A(5); xmin = A(6);
% % % % xmax = A(7); tstart = A(8); tend = A(9); CFL = A(10); Ns= A(11); 
% % % % DiffIon=[A(12) A(13)]; ValIon = [A(14) A(15)]; 
% % % % StericG = [A(16) A(17); A(18) A(19)] 
% % % % AlphaBetaIon(1:2,1:2,1) = [A(20) A(21); A(22) A(23)];
% % % % AlphaBetaIon(1:2,1:2,2) = [A(24) A(25); A(26) A(27)];
% % % % AreaStr='1.0*ones(size(x))';
% % % % NN=[8,12,16]*1;

[p]=parameters;
kB= p.kB; Temp = p.Temp; e_unit = p.e_unit; 
ValIon = p.ValIon; StericG = p.StericG; HOT = p.HOT;
DiffIon = p.DiffIon; epslion = p.epslion; 
IonBc = p.IonBc; PhiBc = p.PhiBc; dt=p.dt; ep=p.ep;
dCdxBc = p.dCdxBc; xmin = p.xmin; xmax = p.xmax; Euler=p.Euler;
tstart = p.tstart; tend = p.tend; 
AlphaBetaIon = p.AlphaBetaIon; AlphaBetaPhi = p.AlphaBetaPhi;
CFL = p.CFL; c_tau = p.c_tau; N = p.NN; savesteps = p.savesteps;

 [Conc] = ModelPnP_PN_ODE15s;
 disp ([Conc phi]);
% %  err1 = max(Conc(:,1) - (1+exp(-1) + cos(pi*x)))
% %  err2 = max(Conc(:,2) - (1+exp(-1) + cos(pi*x)))
%  err1 = max(Conc(:,1) - (ones(size(x)) * exp(-1)))
%  err2 = max(Conc(:,2) - (ones(size(x)) * exp(-1)))  


