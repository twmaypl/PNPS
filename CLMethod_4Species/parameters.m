% Input parameter for PNP-Steric simulation
% Update on 2015/09/28

function [input]=parameters
input.c_tau=2.0;                       %scaling parameter for adjjusting penalty parameters tau
input.NN=100;%[50 75 100 200 500];%50;      % Number of points  
input.CFL=0.35;%/10;%0.35/5;
input.ep = 0.001;                      % Parameter for initial condition
input.savesteps = 1;%1000;                % Output every # of steps
input.dt = 1e-5;

input.xmin     = -1;                    % Left bound
input.xmax     = 1; %1.5;              % Right bound
% SpaceDomain =[xmin xmax];

input.tstart = 0.0;                    % Start time
input.tend   = 20.0;%2.0;             % End time
% TimeInterval =[tstart tend];

input.kB=1.0;
input.Temp=1.0;
input.e_unit =  1.0;
input.dielec = 1.0;
input.epslion = 1E-13;
% PhyConst = [kB Temp e_unit dielec];

% Number of species
input.Ns = 4;       

% Diffusion constants
input.DiffIon = zeros(1,input.Ns);
%Dp = 2.0; Dn = 0.5;
D1 = 1; D2 = 1; D3 = 1; D4=1;
input.DiffIon =[D1 D2 D3 D4];%*10^(-3);

% Valence constants
input.ValIon = zeros(1,input.Ns);
%zp =  2.0; zn = -3.0;
z1 = 1.0; z2= -1.0; z3 = 2.0; z4=-2.0;
% zn = 0;
input.ValIon=[z1 z2 z3 z4];

% Steric effect parameters
input.StericG=zeros(input.Ns, input.Ns,6);
input.StericG(:,:,1) = [1 1.1 0 0; 
                        1.1 1 0 0;
                        0 0 1 1.1;
                        0 0 1.1 1];
input.StericG(:,:,2) = [1 1.2 0 0; 
                        1.2 1 0 0;
                        0 0 1 1.2;
                        0 0 1.2 1];
input.StericG(:,:,3) = [1 1.3 0 0; 
                        1.3 1 0 0;
                        0 0 1 1.3;
                        0 0 1.3 1];
input.StericG(:,:,4) = [1 1.4 0 0; 
                        1.4 1 0 0;
                        0 0 1 1.4;
                        0 0 1.4 1];
input.StericG(:,:,5) = [1 2 0 0; 
                        2 1 0 0;
                        0 0 1 4;
                        0 0 4 1];
input.StericG(:,:,6) = [0 0 0 0; 
                        0 0 0 0;
                        0 0 0 0;
                        0 0 0 0];
                    
% High Order Term parameters
input.Hot=zeros(input.Ns, input.Ns, 2);
input.HOT(:,:,1) = [0 0 0 0;
                    0 0 0 0;
                    0 0 0 0;
                    0 0 0 0];
input.HOT(:,:,2) = [1e-3   0    0    0;
                      0  1e-3   0    0;
                      0    0  1e-2   0;
                      0    0    0  1e-2];


% Boundary Condition
input.ChemBc = 1;   % if 1, BC is defined by chemical potential
                             % if 0, BC is defined by concentration

input.AlphaBetaIon = zeros(2,2,input.Ns);
input.AlphaBetaPhi = zeros(2,2);
input.AlphaBetaIon(1:2,1:2,1)=[1 0; 1 0];
input.AlphaBetaIon(1:2,1:2,2)=[1 0; 1 0];
input.AlphaBetaIon(1:2,1:2,3)=[1 0; 1 0];
input.AlphaBetaIon(1:2,1:2,4)=[1 0; 1 0];
input.AlphaBetaPhi(1:2,1:2) = [1 1e-1; 1 1e-1];

input.AreaStr='1.0*ones(size(x))';
% % % % Boundary Conditions
ion1bc_l = 0.001;%exp(input.xmin); 
ion1bc_r = 0.001;%exp(input.xmax);
ion2bc_l = 0.001;%exp(input.xmin); 
ion2bc_r = 0.001;%exp(input.xmax);
ion3bc_l = 3;
ion3bc_r = 3;
ion4bc_l = 3;%exp(input.xmin); 
ion4bc_r = 3;%exp(input.xmax);

dCdx1bc_l = 0;
dCdx1bc_r = 0;
dCdx2bc_l = 0;
dCdx2bc_r = 0;
dCdx3bc_l = 0;
dCdx3bc_r = 0;
dCdx4bc_l = 0; 
dCdx4bc_r = 0;

ddxAdCdx1bc_l = 0;
ddxAdCdx1bc_r = 0;
ddxAdCdx2bc_l = 0;
ddxAdCdx2bc_r = 0;
ddxAdCdx3bc_l = 0;
ddxAdCdx3bc_r = 0;
ddxAdCdx4bc_l = 0;
ddxAdCdx4bc_r = 0;

phi_l = -1e-3;%-1e-2;%input.xmin; 
phi_r = -1e-3;%1e-2;%input.xmax;
input.IonBc=[ion1bc_l ion2bc_l ion3bc_l ion4bc_l; 
             ion1bc_r ion2bc_r ion3bc_r ion4bc_r];
input.dCdxBc = [dCdx1bc_l dCdx2bc_l dCdx3bc_l dCdx4bc_l;
                dCdx1bc_r dCdx2bc_r dCdx3bc_r dCdx4bc_r];
input.ddxAdCdxBc = [ddxAdCdx1bc_l ddxAdCdx2bc_l ddxAdCdx3bc_l ddxAdCdx4bc_l;
                    ddxAdCdx1bc_r ddxAdCdx2bc_r ddxAdCdx3bc_r ddxAdCdx4bc_r];
 

input.PhiBc=[phi_l phi_r];
% disp(input.IonBc)
% disp([input.Pbc_l input.Nbc_l input.Pbc_r input.Nbc_r])

 
