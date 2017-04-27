% Input parameter for PNP-Steric simulation
% Update on 2015/09/03

function [input]=parameters
input.c_tau=5.0;                       %scaling parameter for adjjusting penalty parameters tau
input.NN=100;%[50 75 100 200 500];%50;      % Number of points  
input.CFL=0.3/55;%/10;%0.35/5;
input.ep = 0.001;                      % Parameter for initial condition
input.savesteps = 1;                % Output every # of steps
input.dt = 1e-3;%0.1;

input.xmin     = -1;                    % Left bound
input.xmax     = 1; %1.5;              % Right bound
% SpaceDomain =[xmin xmax];

input.tstart = 0.0;                    % Start time
input.tend   = 2e-2;%2.0;             % End time
% TimeInterval =[tstart tend];

input.kB=1.0;
input.Temp=1.0;
input.e_unit =  1.0;
input.dielec = 1.0;
input.epslion = 1E-13;
% PhyConst = [kB Temp e_unit dielec];

% Number of species
input.Ns = 2;       

% Diffusion constants
input.DiffIon = zeros(1,input.Ns);
%Dp = 2.0; Dn = 0.5;
Dp = 1; Dn = 0;
input.DiffIon =[Dp Dn];%*10^(-3);

% Valence constants
input.ValIon = zeros(1,input.Ns);
%zp =  2.0; zn = -3.0;
zp = 2.0; zn= -2.0;
% zn = 0;
input.ValIon=[zp zn];

% Steric effect parameters
input.StericG=zeros(input.Ns, input.Ns,4);
input.StericG(:,:,1)=[0 0; ...
                      0 0];
input.StericG(:,:,2)=[1 0.5; ...
                      0.5 1];
input.StericG(:,:,3)=[1 1.0; ...
                      1.0 1];
input.StericG(:,:,4)=[1 4; ...
                      4 1];


% High order term parameters
input.HOT=zeros(input.Ns, input.Ns,5);
input.HOT(:,:,1) = 1e-2 * [1 0;
                           0 1];
input.HOT(:,:,2) = [5e-3 0; ...
                    0 5e-3];
input.HOT(:,:,3) = [5e-1 0; ...
                    0 5e-1];
input.HOT(:,:,4) = [5e-4 0; ...
                    0 5e-4];
input.HOT(:,:,5) = [0 0; ...
                    0 0];                

% Boundary Condition
input.ChemBc = 1; 
input.AlphaBetaIon = zeros(2,2,input.Ns);
input.AlphaBetaPhi = zeros(2,2);
input.AlphaBetaIon(1:2,1:2,1)=[1 0; 1 0];
input.AlphaBetaIon(1:2,1:2,2)=[1 0; 1 0];
input.AlphaBetaPhi(1:2,1:2) = [1 1e-1; 1 1e-1];

input.AreaStr='1.0*ones(size(x))';
% % % % Boundary Conditions
input.ChemBc = 1;
ion1bc_l = 3;%exp(input.xmin); 
ion1bc_r = 3;%exp(input.xmax);
ion2bc_l = 3;%exp(input.xmin); 
ion2bc_r = 3;%exp(input.xmax);

dCdx1bc_l = 0;
dCdx1bc_r = 0;
dCdx2bc_l = 0;
dCdx2bc_r = 0; 

ddxAdCdx1bc_l = 0;
ddxAdCdx1bc_r = 0;
ddxAdCdx2bc_l = 0;
ddxAdCdx2bc_r = 0;

phi_l = -1e-4; %-1e-4;%input.xmin; 
phi_r = -1e-4; %1e-4;%input.xmax;

input.IonBc=[ion1bc_l ion2bc_l; ion1bc_r ion2bc_r];
input.dCdxBc = [dCdx1bc_l dCdx2bc_l; dCdx1bc_r dCdx2bc_r];
input.PhiBc=[phi_l phi_r]; 
input.ddxAdCdxBc = [ddxAdCdx1bc_l ddxAdCdx2bc_l; ddxAdCdx1bc_r ddxAdCdx2bc_r]; 
% disp([input.Pbc_l input.Nbc_l input.Pbc_r input.Nbc_r])

 
