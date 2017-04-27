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
global DiffIon ValIon ep StericG HOT
global AlphaBetaIon AreaStr 
global IonBc PhiBc dCdxBc ddxAdCdxBc ChemBc
global xmin xmax tstart tend Euler
global N CFL c_tau AlphaBetaPhi

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
kB= p.kB; Temp = p.Temp; e_unit = p.e_unit; dielec = p.dielec;
ValIon = p.ValIon; TotStericG = p.StericG(:,:,:); TotHOT = p.HOT(:,:,:);
DiffIon = p.DiffIon; AreaStr = p.AreaStr; epslion = p.epslion; 
IonBc = p.IonBc; PhiBc = p.PhiBc; dt=p.dt; ep=p.ep;
dCdxBc = p.dCdxBc; ddxAdCdxBc = p.ddxAdCdxBc; ChemBc = p.ChemBc; 
xmin = p.xmin; xmax = p.xmax; Euler=p.Euler;
tstart = p.tstart; tend = p.tend; 
AlphaBetaIon = p.AlphaBetaIon; AlphaBetaPhi = p.AlphaBetaPhi;
CFL = p.CFL; c_tau = p.c_tau; NN = p.NN; savesteps = p.savesteps;

err=zeros(length(NN),length(TotStericG(1,1,:)));

for IndxStericG = 2:2%1:length(TotStericG(1,1,:))
    StericG = TotStericG(:,:,IndxStericG)
    for IndxHOT = 5:5% 1:length(TotHOT)
        HOT = TotHOT(:,:,IndxHOT)
        for i = 1:length(NN)
            N = NN(i);
            ModelPnP_PN_ODE15s(dt);
%             [err(i,j)]= ModelPnP_PN(Bc);
        end
%         loglog(NN, err(:,j)), hold all
    end
end
% save('ERR_Boundary.mat','NN','err');



